#!/usr/bin/env Rscript
# =============================================================================
# qc_report.R: input QC, assertions, manifest and markdown
# -----------------------------------------------------------------------------
# Sourced after lib/module1_prep.R, whose helpers it reuses (logf, write_json,
# file_provenance, m1_to_json, BROAD_CATS). For each cancer it:
#   - computes QC metrics (filter funnel, per-bucket counts, expression,
#     per-chromosome distribution, library sizes, sparsity, DeepLoc confidence),
#   - runs machine-checkable assertions on the output SE,
#   - writes <slug>_tier1_qc.md and <slug>_tier1_manifest.json,
#   - returns TRUE if no hard assertion failed (the entrypoint exits non-zero
#     otherwise).
# =============================================================================

suppressPackageStartupMessages({ library(SummarizedExperiment) })

# ── small helpers ────────────────────────────────────────────────────────────
pkg_versions <- function() {
    getv <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
    list(R = R.version.string,
         data.table           = getv('data.table'),
         SummarizedExperiment = getv('SummarizedExperiment'),
         satuRn               = getv('satuRn'),    # downstream engine (provenance)
         edgeR                = getv('edgeR'))
}

fmt_num <- function(x, d = 2) formatC(x, format = 'f', digits = d)

# markdown table from a data.frame (all cells coerced to character; pipes/newlines
# sanitized so free-text detail strings can't break the table layout)
md_table <- function(df) {
    df[] <- lapply(df, function(col) {
        v <- format(col, trim = TRUE)
        gsub('\n', ' ', gsub('|', '\\|', v, fixed = TRUE))
    })
    cn  <- colnames(df)
    hdr <- paste0('| ', paste(cn, collapse = ' | '), ' |')
    sep <- paste0('| ', paste(rep('---', length(cn)), collapse = ' | '), ' |')
    rows <- if (nrow(df))
        apply(df, 1, function(r) paste0('| ', paste(trimws(r), collapse = ' | '), ' |')) else character()
    paste(c(hdr, sep, rows), collapse = '\n')
}

assert <- function(name, pass, level = 'hard', detail = '')
    list(name = name, pass = isTRUE(pass), level = level, detail = as.character(detail))

assertions_md <- function(asserts) {
    df <- data.frame(
        check  = vapply(asserts, `[[`, character(1), 'name'),
        level  = vapply(asserts, `[[`, character(1), 'level'),
        status = vapply(asserts, function(a) if (a$pass) 'PASS' else 'FAIL', character(1)),
        detail = vapply(asserts, `[[`, character(1), 'detail'),
        stringsAsFactors = FALSE)
    md_table(df)
}

hard_failed <- function(asserts)
    any(vapply(asserts, function(a) a$level == 'hard' && !a$pass, logical(1)))

# gene -> chromosome lookup from the prep tx2gene table (if 'chr' present)
gene_chr_table <- function(prep, gene_ids) {
    if (!'chr' %in% colnames(prep$t2g)) return(NULL)
    g2c <- prep$t2g[!duplicated(prep$t2g$gene_id), c('gene_id', 'chr')]
    chr_of <- setNames(as.character(g2c$chr), g2c$gene_id)
    sort(table(chr_of[gene_ids]), decreasing = TRUE)
}

# =============================================================================
# Metrics + assertions + markdown
# =============================================================================
compute_qc_tier1 <- function(prep, se) {
    A   <- assay(se, 'counts')
    rd  <- as.data.frame(rowData(se))
    ktf <- prep$keep_tx_final
    gtx <- prep$g[ktf]; btx <- prep$tx_bucket[ktf]

    # buckets-per-gene over the DeepLoc-annotated set (pre multi-bucket drop)
    pr   <- paste(gtx, btx, sep = '|'); first <- !duplicated(pr)
    nbk  <- table(gtx[first])
    bpg  <- table(factor(as.integer(nbk), levels = 1:3))            # genes with 1 / 2 / 3 buckets

    # DeepLoc argmax confidence over the kept transcripts
    Pk  <- prep$P[ktf, , drop = FALSE]
    pmx <- Pk[cbind(seq_len(nrow(Pk)), max.col(Pk, ties.method = 'first'))]

    # library sizes by condition
    lib  <- prep$lib_sizes
    cnd  <- as.character(prep$col_meta[names(lib), 'condition'])
    dep  <- tapply(lib, cnd, median)

    # degenerate / one-condition features (note: avoid shadowing base::ncol)
    tcols <- prep$col_meta$sample[prep$col_meta$condition == 'Tumor']
    ncols <- prep$col_meta$sample[prep$col_meta$condition == 'Normal']
    zeroT <- rowSums(A[, tcols, drop = FALSE]) == 0
    zeroN <- rowSums(A[, ncols, drop = FALSE]) == 0

    list(
        n_features = nrow(se), n_genes = length(unique(rd$gene_id)), n_samples = ncol(se),
        bucket_counts = table(factor(rd$bucket, levels = BROAD_CATS)),
        buckets_per_gene = bpg,
        feat_mean_count = summary(rowMeans(A)),
        lib_median_M = c(Normal = unname(dep['Normal']) / 1e6, Tumor = unname(dep['Tumor']) / 1e6),
        tn_depth_ratio = unname(dep['Tumor'] / dep['Normal']),
        sparsity = mean(A == 0),
        deeploc_pmax_median = median(pmx), deeploc_pmax_frac_lt05 = mean(pmx < 0.5),
        chr_table = gene_chr_table(prep, rd$gene_id),
        flag_empty_buckets = setdiff(BROAD_CATS, rd$bucket),
        flag_degenerate_features = sum(zeroT & zeroN),
        flag_one_condition_features = sum(xor(zeroT, zeroN)),
        flag_depth_imbalance = abs(log2(unname(dep['Tumor'] / dep['Normal'])))
    )
}

assertions_tier1 <- function(prep, se) {
    A  <- assay(se, 'counts')
    rd <- as.data.frame(rowData(se))
    cd <- as.data.frame(colData(se))
    out <- list()

    out <- c(out, list(assert('dims_positive', nrow(se) > 0 && ncol(se) > 0, 'hard',
                              sprintf('%d features x %d samples', nrow(se), ncol(se)))))
    out <- c(out, list(assert('counts_finite', is.numeric(A) && all(is.finite(A)), 'hard')))
    out <- c(out, list(assert('counts_no_na', !anyNA(A), 'hard')))
    out <- c(out, list(assert('counts_nonneg', min(A) >= 0, 'hard', sprintf('min = %s', min(A)))))
    out <- c(out, list(assert('rowdata_schema',
                              all(c('isoform_id', 'gene_id', 'bucket') %in% colnames(rd)), 'hard')))
    tg <- table(rd$gene_id)
    out <- c(out, list(assert('genes_ge2_features', all(tg >= 2), 'hard',
                              sprintf('min features/gene = %d', min(tg)))))

    # mass conservation across ALL genes (disjoint bucket sums == source tx sums)
    se_g  <- rowsum(A, rd$gene_id)
    g_ktf <- prep$g[prep$keep_tx_final]
    src_g <- rowsum(prep$count_mat[prep$keep_tx_final, , drop = FALSE], g_ktf)
    src_g <- src_g[rownames(se_g), colnames(se_g), drop = FALSE]
    drift     <- max(abs(se_g - src_g))
    rel_drift <- max(abs(se_g - src_g) / pmax(abs(src_g), 1))   # scale-free; with fractional (double)
                                                                #   counts mass conserves only to float
                                                                #   precision (~1e-10), so test a tolerance,
                                                                #   not exact 0. A real mis-assignment is ~O(1).
    out <- c(out, list(assert('mass_conservation', rel_drift < 1e-6, 'hard',
                              sprintf('max |drift| = %.3g, max relative drift = %.2e (tol 1e-6)', drift, rel_drift))))

    out <- c(out, list(assert('coldata_schema',
                              all(c('condition', 'patient_id') %in% colnames(cd)), 'hard')))
    tab <- table(cd$patient_id, cd$condition)
    balanced <- all(tab[, 'Normal'] == 1) && all(tab[, 'Tumor'] == 1)
    out <- c(out, list(assert('balanced_pairs', balanced, 'hard',
                              sprintf('%d patients, each 1 Normal + 1 Tumor', nrow(tab)))))
    out
}

render_markdown_tier1 <- function(prep, se, m, asserts) {
    f <- prep$funnel
    funnel_df <- data.frame(
        step  = c('multi-iso genes (>=2 tx)', 'gene CI>10 in BOTH N&T (f2)',
                  '>=1 tx CI>10 in EITHER (f3)', 'ISAR-pass genes (f1&f2&f3)',
                  'ISAR-pass transcripts', 'DeepLoc-annotated genes',
                  'DeepLoc-annotated transcripts', 'Tier-1 testable genes (>=2 buckets)',
                  'Tier-1 features (gene x bucket)'),
        count = c(f$f1_ge2tx_genes, f$f2_bothNT_genes, f$f3_eitherTx_genes, f$isar_pass_genes,
                  f$isar_pass_tx, f$deeploc_annot_genes, f$deeploc_annot_tx,
                  m$n_genes, m$n_features),
        stringsAsFactors = FALSE)

    bucket_df <- data.frame(bucket = names(m$bucket_counts),
                            features = as.integer(m$bucket_counts), stringsAsFactors = FALSE)
    bpg_df <- data.frame(buckets_per_gene = names(m$buckets_per_gene),
                         genes = as.integer(m$buckets_per_gene), stringsAsFactors = FALSE)
    expr <- m$feat_mean_count
    expr_df <- data.frame(statistic = names(expr), value = fmt_num(as.numeric(expr), 1),
                          stringsAsFactors = FALSE)

    chr_md <- if (is.null(m$chr_table)) '_tx2gene has no `chr` column, per-chromosome view skipped._'
              else md_table(data.frame(chromosome = names(m$chr_table),
                                       features = as.integer(m$chr_table), stringsAsFactors = FALSE))

    flags <- c(
        sprintf('- Empty buckets (no feature): %s',
                if (length(m$flag_empty_buckets)) paste(m$flag_empty_buckets, collapse = ', ') else 'none'),
        sprintf('- Degenerate features (all-zero across all samples): %d', m$flag_degenerate_features),
        sprintf('- One-condition-only features (all-zero in N xor T): %d', m$flag_one_condition_features),
        sprintf('- Depth imbalance |log2(median T / median N)|: %s', fmt_num(m$flag_depth_imbalance, 2)))

    paste(c(
        sprintf('# Input QC: %s (`%s`)', prep$cancer, prep$slug),
        '',
        sprintf('Data preparation · v%s · generated %s',
                MODULE1_VERSION, format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
        '',
        '**Tier-1 unit:** `gene__bucket` features (disjoint DeepLoc-bucket sums), genes with >=2 occupied buckets.',
        '',
        '## Cohort',
        sprintf('- Samples: %d Tumor + %d Normal (%d matched pairs)',
                length(prep$tum_samples), length(prep$nor_samples), prep$n_pairs),
        sprintf('- Library depth (median M reads/sample): Normal=%s  Tumor=%s  (T/N=%s)',
                fmt_num(m$lib_median_M['Normal'], 1), fmt_num(m$lib_median_M['Tumor'], 1),
                fmt_num(m$tn_depth_ratio, 2)),
        '',
        '## Filter funnel',
        md_table(funnel_df),
        '',
        '## Features per localization bucket',
        md_table(bucket_df),
        '',
        '## Occupied buckets per gene (annotated set)',
        md_table(bpg_df),
        '',
        '## Feature expression (mean raw count per gene x bucket)',
        md_table(expr_df),
        sprintf('- Assay sparsity (fraction of zero cells): %s', fmt_num(m$sparsity, 3)),
        sprintf('- DeepLoc argmax confidence: median P_max=%s, %s%% below 0.5',
                fmt_num(m$deeploc_pmax_median, 2), fmt_num(100 * m$deeploc_pmax_frac_lt05, 0)),
        '',
        '## Per-chromosome feature distribution',
        chr_md,
        '',
        '## Data-quality flags',
        paste(flags, collapse = '\n'),
        '',
        '## Acceptance assertions (machine-checkable)',
        assertions_md(asserts),
        '',
        sprintf('**QC verdict:** %s', if (hard_failed(asserts)) 'FAIL (hard assertion failed)' else 'PASS'),
        ''), collapse = '\n')
}


# =============================================================================
# MANIFEST + finalizers (called by the entrypoints)
# =============================================================================
build_manifest <- function(tier, prep, outputs, asserts, full_md5, summary = NULL) {
    big <- c('tumor_log2', 'normal_log2')                           # ~1 GB, md5 only on demand
    inputs <- lapply(names(prep$files), function(k)
        c(list(name = k), file_provenance(prep$files[[k]], md5 = (full_md5 || !(k %in% big)))))
    list(
        module = 'module1_data_prep', version = MODULE1_VERSION, tier = tier,
        generated_at = format(Sys.time(), '%Y-%m-%dT%H:%M:%S'),
        cancer = prep$cancer, slug = prep$slug,
        summary = summary,                                         # headline numbers for cohort aggregation
        versions = pkg_versions(),
        thresholds = prep$thresholds,
        cohort = list(n_tumor = length(prep$tum_samples), n_normal = length(prep$nor_samples),
                      n_pairs = prep$n_pairs,
                      tumor_samples = prep$tum_samples, normal_samples = prep$nor_samples,
                      patient_ids = levels(prep$col_meta$patient_id)),
        filter_funnel = prep$funnel,
        inputs = inputs,
        outputs = outputs,
        assertions = lapply(asserts, function(a)
            list(name = a$name, level = a$level, pass = a$pass, detail = a$detail)),
        qc_pass = !hard_failed(asserts)
    )
}

finalize_tier1 <- function(prep, se, out_dir, full_md5 = FALSE) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    slug    <- prep$slug
    se_path <- file.path(out_dir, sprintf('%s_tier1_se.rds', slug))
    saveRDS(se, se_path)
    m       <- compute_qc_tier1(prep, se)
    asserts <- assertions_tier1(prep, se)
    writeLines(render_markdown_tier1(prep, se, m, asserts),
               file.path(out_dir, sprintf('%s_tier1_qc.md', slug)))
    n_genes <- length(unique(as.data.frame(rowData(se))$gene_id))
    outputs <- list(list(file = basename(se_path), md5 = unname(tools::md5sum(se_path)),
                         n_features = nrow(se), n_genes = n_genes, n_samples = ncol(se)))
    summ <- list(n_features = nrow(se), n_genes = n_genes, n_samples = ncol(se))
    write_json(build_manifest('tier1', prep, outputs, asserts, full_md5, summary = summ),
               file.path(out_dir, sprintf('%s_tier1_manifest.json', slug)))
    logf('Written: %s (%d features / %d genes), QC %s',
         se_path, nrow(se), length(unique(as.data.frame(rowData(se))$gene_id)),
         if (hard_failed(asserts)) 'FAIL' else 'PASS')
    invisible(!hard_failed(asserts))
}

