#!/usr/bin/env Rscript
# =============================================================================
# qc_results.R: DTU result QC + manifest
# -----------------------------------------------------------------------------
# Sourced after lib/module1_prep.R, lib/qc_report.R and lib/saturn_dtu.R, whose
# helpers it reuses (logf, write_json, file_provenance, md_table, fmt_num,
# pkg_versions, assert, assertions_md, hard_failed, SATURN_DTU_VERSION).
# Defines functions only.
#
# Same shape as qc_report.R (compute_* -> assertions_* -> render_* ->
# finalize_* + build_manifest); the content is p-value and FDR QC of the satuRn
# result table rather than input-data QC. Significance is assessed on
# `regular_FDR`, the Benjamini-Hochberg adjustment of the raw p-value, at 0.05.
#
# Reported diagnostics:
#
#   * Genomic inflation factor
#         lambda = qchisq(median(p), 1, lower=F) / qchisq(0.5, 1, lower=F)
#     Under a mostly true null the chi^2_1 statistics follow their reference
#     distribution, so lambda ~ 1. Values above 1 indicate an over-dispersed
#     null or dense signal; lambda_dedup, computed on one minimum p-value per
#     gene, separates many-features-of-few-genes inflation from gene-level
#     signal.
#
#   * pi0 ~ 2 * mean(p > 0.5): Storey-style estimate of the null fraction. The
#     upper half of a Uniform(0,1) p-value histogram is nearly all null, so
#     doubling its mass estimates the null proportion; 1 - pi0 ~ signal fraction.
# =============================================================================

suppressPackageStartupMessages({ library(SummarizedExperiment) })

# ── genomic inflation factor from a vector of p-values ───────────────────────
gc_lambda <- function(p) {
    p <- p[is.finite(p) & p > 0 & p < 1]
    if (length(p) < 5) return(NA_real_)
    qchisq(median(p), df = 1, lower.tail = FALSE) / qchisq(0.5, df = 1, lower.tail = FALSE)
}

# =============================================================================
# METRICS: diagnostics + yield, computed from the result data.frame.
#   `res` must carry: pval, regular_FDR, estimates, t, feature_id, gene_id
#   (run_saturn_dtu_tier1.R attaches gene_id).
# =============================================================================
compute_qc_results <- function(res, fdr = 0.05, n_dropped = 0L) {
    gene <- res$gene_id
    p    <- res$pval
    fdrR <- res$regular_FDR

    hit_counts <- function(v) {                              # features + genes below `fdr`
        ok <- which(is.finite(v) & v < fdr)
        c(feat = length(ok), gene = length(unique(gene[ok])))
    }
    hR <- hit_counts(fdrR)

    list(
        n_features = nrow(res),
        n_genes    = length(unique(gene)),
        n_dropped  = as.integer(n_dropped),
        n_na_pval  = sum(is.na(p)),
        n_na_est   = sum(is.na(res$estimates)),
        n_zero_est = sum(res$estimates == 0, na.rm = TRUE),
        lambda_raw      = gc_lambda(p),
        lambda_dedup    = gc_lambda(suppressWarnings(                       # 1 min-p per gene
            tapply(p, gene, function(x) min(x, na.rm = TRUE)))),
        pi0             = min(1, 2 * mean(p > 0.5, na.rm = TRUE)),
        fdr             = fdr,
        hits_regFDR_feat = unname(hR['feat']), hits_regFDR_gene = unname(hR['gene'])
    )
}

# ── integrity assertions on the result table ─────────────────────────────────
RESULT_COLUMNS <- c('feature_id', 'gene_id', 'bucket', 'estimates', 'se', 'df',
                    't', 'pval', 'regular_FDR', 'q_BH', 'side_log2FC')

assertions_results <- function(res, m) {
    out <- list()
    out <- c(out, list(assert('results_nonempty', nrow(res) > 0, 'hard',
                              sprintf('%d features', nrow(res)))))
    out <- c(out, list(assert('pval_not_all_na', m$n_na_pval < nrow(res), 'hard',
                              sprintf('%d / %d NA pval', m$n_na_pval, nrow(res)))))
    out <- c(out, list(assert('expected_columns',
                              identical(colnames(res), RESULT_COLUMNS), 'hard',
                              paste(colnames(res), collapse = ', '))))
    out
}

# =============================================================================
# p-value / FDR panel: raw p and regular_FDR, the two distributions behind the
#   reported calls. Raw p is annotated with lambda; both with hits < the cutoff
#   (features | genes). Base graphics.
# =============================================================================
plot_dtu_pfdr <- function(res, png_path, title_prefix, fdr = 0.05) {
    gene <- res$gene_id
    panels <- list(
        list(nm = 'raw p',       v = res$pval,        col = '#6C8EBF', lam = TRUE),
        list(nm = 'regular_FDR', v = res$regular_FDR, col = '#9673A6', lam = FALSE))

    # cairo device works on headless compute nodes (no X11/DISPLAY needed).
    dev_type <- if (isTRUE(capabilities('cairo'))) 'cairo' else getOption('bitmapType', 'Xlib')
    grDevices::png(png_path, width = 1150, height = 460, res = 120, type = dev_type)
    on.exit(grDevices::dev.off(), add = TRUE)
    op <- par(mfrow = c(1, 2), mar = c(4.2, 4, 3.6, 1)); on.exit(par(op), add = TRUE)
    for (x in panels) {
        v  <- x$v; ok <- is.finite(v)
        nf <- sum(v < fdr, na.rm = TRUE)
        ng <- length(unique(gene[which(v < fdr)]))
        ttl <- if (x$lam)
            sprintf('%s %s (lambda=%.2f)\n<%.2g: %d feat / %d genes', title_prefix, x$nm, gc_lambda(v), fdr, nf, ng)
        else
            sprintf('%s %s\n<%.2g: %d feat / %d genes', title_prefix, x$nm, fdr, nf, ng)
        if (any(ok)) hist(v[ok], breaks = 40, col = x$col, border = 'white',
                          xlab = x$nm, cex.main = 0.95, main = ttl)
        else { plot.new(); title(main = paste(ttl, '\n(all NA)')) }
    }
    invisible(png_path)
}

# =============================================================================
# Markdown report: diagnostics, yield and totals.
# =============================================================================
render_results_qc_md <- function(meta, m, asserts, png_basename) {
    cal_df <- data.frame(
        metric = c('lambda_raw', 'lambda_dedup (per-gene min-p)', 'pi0 (~2*mean(p>0.5))'),
        value  = c(fmt_num(m$lambda_raw, 3), fmt_num(m$lambda_dedup, 3), fmt_num(m$pi0, 3)),
        stringsAsFactors = FALSE)

    yield_df <- data.frame(
        cutoff   = sprintf('regular_FDR < %.2f', m$fdr),
        features = m$hits_regFDR_feat,
        genes    = m$hits_regFDR_gene,
        stringsAsFactors = FALSE)

    totals_df <- data.frame(
        quantity = c('features tested', 'genes tested', 'degenerate fits dropped (pre-testDTU)',
                     'NA p-values', 'NA estimates', 'zero estimates'),
        value    = c(m$n_features, m$n_genes, m$n_dropped, m$n_na_pval, m$n_na_est, m$n_zero_est),
        stringsAsFactors = FALSE)


    paste(c(
        sprintf('# %s DTU result QC: %s (`%s`)', meta$tier_label, meta$cancer, meta$slug),
        '',
        sprintf('satuRn DTU testing · v%s · generated %s',
                SATURN_DTU_VERSION, format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
        '',
        sprintf('**Engine:** satuRn, paired block design `~ 0 + condition + patient_id`, contrast `conditionTumor - conditionNormal`.'),
        sprintf('**Test unit:** %s', meta$unit_desc),
        '',
        '## p-value diagnostics',
        md_table(cal_df),
        '',
        '- **lambda_raw**: genomic inflation of the raw p-values. Above 1 indicates an over-dispersed null or dense signal.',
        '- **lambda_dedup**: the same statistic on the per-gene *minimum* p-value. Taking the minimum over a gene\'s features biases it upward, more so for genes with many features.',
        '- **pi0**: null fraction estimated from the raw p-values.',
        '',
        '## Yield (FDR-corrected hits)',
        md_table(yield_df),
        '',
        '_`regular_FDR` is the Benjamini-Hochberg adjustment of the raw satuRn p-value._',
        '',
        '## Totals & fit health',
        md_table(totals_df),
        '',
        sprintf('## p-value / FDR panel'),
        sprintf('![p-value / FDR panel](%s)', png_basename),
        '',
        '## Acceptance assertions (machine-checkable)',
        assertions_md(asserts),
        '',
        sprintf('**QC verdict:** %s', if (hard_failed(asserts)) 'FAIL (hard assertion failed)' else 'PASS'),
        ''), collapse = '\n')
}

# =============================================================================
# Manifest + finalizer, called by run_saturn_dtu_tier1.R.
# =============================================================================
build_dtu_manifest <- function(tier, meta, m, asserts, inputs, outputs) {
    list(
        module = 'saturn_dtu', version = SATURN_DTU_VERSION, tier = tier,
        generated_at = format(Sys.time(), '%Y-%m-%dT%H:%M:%S'),
        cancer = meta$cancer, slug = meta$slug,
        engine = list(tool = 'satuRn', formula = '~ 0 + condition + patient_id',
                      contrast = 'conditionTumor - conditionNormal',
                      design = 'paired randomized complete block (block = patient_id)',
                      significance_scale = 'regular_FDR (Benjamini-Hochberg on the raw p-value)'),
        versions = pkg_versions(),
        thresholds = list(fdr = m$fdr),
        summary = list(n_features = m$n_features, n_genes = m$n_genes, n_dropped = m$n_dropped,
                       lambda_raw = m$lambda_raw, lambda_dedup = m$lambda_dedup, pi0 = m$pi0,
                       hits_regFDR_features = m$hits_regFDR_feat, hits_regFDR_genes = m$hits_regFDR_gene),
        inputs = inputs,
        outputs = outputs,
        assertions = lapply(asserts, function(a)
            list(name = a$name, level = a$level, pass = a$pass, detail = a$detail)),
        qc_pass = !hard_failed(asserts)
    )
}

#' Write results.{tsv,rds} + manifest.json + qc.md + the p-value/FDR PNG.
#' @param tier      'tier1'
#' @param meta      list(cancer, slug, tier_label, unit_desc)
#' @param res       finalized result data.frame (satuRn columns + gene_id +
#'                   bucket + side_log2FC).
#' @param out_dir   the existing output/<slug>/<tier>/ directory.
#' @param input_files named character vector of the source se.rds paths.
#' @param fdr,n_dropped  QC params.
#' @return invisibly TRUE if no hard assertion failed.
finalize_dtu_results <- function(tier, meta, res, out_dir, input_files, fdr = 0.05, n_dropped = 0L) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    slug <- meta$slug
    base <- sprintf('%s_%s_dtu', slug, tier)

    tsv_path <- file.path(out_dir, sprintf('%s_results.tsv', base))
    rds_path <- file.path(out_dir, sprintf('%s_results.rds', base))
    png_path <- file.path(out_dir, sprintf('%s_qc_pfdr.png', base))
    md_path  <- file.path(out_dir, sprintf('%s_qc.md',       base))
    json_path <- file.path(out_dir, sprintf('%s_manifest.json', base))

    data.table::fwrite(res, tsv_path, sep = '\t')
    saveRDS(res, rds_path)

    m       <- compute_qc_results(res, fdr = fdr, n_dropped = n_dropped)
    asserts <- assertions_results(res, m)
    plot_dtu_pfdr(res, png_path, title_prefix = meta$tier_label, fdr = fdr)
    writeLines(render_results_qc_md(meta, m, asserts, basename(png_path)), md_path)

    inputs <- lapply(names(input_files), function(k)
        c(list(name = k), file_provenance(input_files[[k]], md5 = TRUE)))
    out_files <- c(tsv_path, rds_path, png_path, md_path)
    outputs <- lapply(out_files, function(p)
        list(file = basename(p), md5 = unname(tools::md5sum(p))))
    write_json(build_dtu_manifest(tier, meta, m, asserts, inputs, outputs), json_path)

    logf('%s DTU written: %s (%d features / %d genes; dropped=%d; lambda_raw=%.2f lambda_dedup=%.2f) - QC %s',
         meta$tier_label, tsv_path, m$n_features, m$n_genes, m$n_dropped,
         m$lambda_raw, m$lambda_dedup, if (hard_failed(asserts)) 'FAIL' else 'PASS')
    invisible(!hard_failed(asserts))
}
