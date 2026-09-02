#!/usr/bin/env Rscript
# =============================================================================
# cohort_select.R: cohort selection + cancer-type stats
# -----------------------------------------------------------------------------
# Two modes:
#   (default)      Select TCGA cancers with >= --min-paired paired (Primary-Tumor
#                  + Solid-Tissue-Normal) patients. Writes, under <output-dir>:
#                    _cohort/cohort_manifest.tsv
#                    _cohort/cohort_summary.md
#                    _cohort/cohort_plots.pdf
#   --summarize    Re-read each passing cancer's manifest (after make_tier1_input.R
#                  has been run for every cancer) and append #testable
#                  genes/features to the summary table + plots.
#
# Usage:
#   Rscript cohort_select.R                    # select (metadata only, fast)
#   Rscript cohort_select.R --min-paired 10
#   Rscript cohort_select.R --summarize        # after make_tier1_input.R has run
# =============================================================================

.script_dir <- local({
    ca <- commandArgs(FALSE); m <- grepl('^--file=', ca)
    if (any(m)) normalizePath(dirname(sub('^--file=', '', ca[which(m)[1]]))) else normalizePath('.')
})
source(file.path(.script_dir, 'lib', 'module1_prep.R'))   # data.table, slugify, logf, write_json, m1_parse_cli

.repo_root <- dirname(.script_dir)
spec <- list(
    `input-dir`  = list(default = file.path(.repo_root, 'data', 'inputs'),  help = 'directory with metadata.tsv.gz'),
    `output-dir` = list(default = file.path(.repo_root, 'data', 'outputs'), help = 'root output directory'),
    `min-paired` = list(default = '10', help = 'minimum paired patients to qualify'),
    `summarize`  = list(default = FALSE, help = 'aggregate per-cancer tier manifests (post fan-out)', flag = TRUE)
)

# local markdown table (keeps this script independent of the SE stack)
md_table <- function(df) {
    df[] <- lapply(df, as.character); cn <- colnames(df)
    rows <- if (nrow(df)) apply(df, 1, function(r) paste0('| ', paste(trimws(r), collapse = ' | '), ' |')) else character()
    paste(c(paste0('| ', paste(cn, collapse = ' | '), ' |'),
            paste0('| ', paste(rep('---', length(cn)), collapse = ' | '), ' |'), rows), collapse = '\n')
}

# ── paired-patient counts per cancer (same pairing rule as the SE builder) ──
cohort_table <- function(metadata_path, min_paired) {
    meta <- fread(metadata_path)
    meta <- meta[`_sample_type` %in% c('Primary Tumor', 'Solid Tissue Normal')]
    meta[, patient_id := substr(sample, 1, 12)]
    pat <- meta[, .(nt = sum(`_sample_type` == 'Primary Tumor'),
                    nn = sum(`_sample_type` == 'Solid Tissue Normal')),
                by = .(detailed_category, patient_id)]
    pat[, paired := nt >= 1 & nn >= 1]
    tab <- pat[, .(n_paired = sum(paired),
                   n_patients_with_tumor  = sum(nt >= 1),
                   n_patients_with_normal = sum(nn >= 1)), by = detailed_category]
    smp <- meta[, .(n_tumor_samples  = sum(`_sample_type` == 'Primary Tumor'),
                    n_normal_samples = sum(`_sample_type` == 'Solid Tissue Normal')),
                by = detailed_category]
    tab <- merge(tab, smp, by = 'detailed_category')
    tab[, slug := vapply(detailed_category, slugify, character(1))]
    tab[, pass := n_paired >= min_paired]
    setorder(tab, -n_paired)
    tab[]
}

# ── read back the headline numbers a finalize_* manifest stores under "summary" ──
manifest_summary <- function(path) {
    keys <- c('n_features', 'n_genes', 'n_samples')
    out <- setNames(rep(NA_real_, length(keys)), keys)
    if (!file.exists(path)) return(out)
    if (requireNamespace('jsonlite', quietly = TRUE)) {
        s <- tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE)$summary, error = function(e) NULL)
        if (!is.null(s)) for (k in keys) if (!is.null(s[[k]])) out[k] <- as.numeric(s[[k]])
        return(out)
    }
    txt <- paste(readLines(path, warn = FALSE), collapse = '\n')      # fallback: targeted regex
    for (k in keys) {
        mm <- regmatches(txt, regexpr(sprintf('"%s"[[:space:]]*:[[:space:]]*[0-9.eE+-]+', k), txt))
        if (length(mm)) out[k] <- as.numeric(sub('.*:[[:space:]]*', '', mm))
    }
    out
}

main <- function() {
    opts <- m1_parse_cli(commandArgs(TRUE), spec, prog = 'cohort_select.R')
    min_paired <- as.integer(opts[['min-paired']])
    metadata_path <- file.path(opts[['input-dir']], 'metadata.tsv.gz')
    if (!file.exists(metadata_path)) stop(sprintf('metadata not found: %s', metadata_path))
    cohort_dir <- file.path(opts[['output-dir']], '_cohort')
    dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)

    tab <- cohort_table(metadata_path, min_paired)
    npass <- sum(tab$pass)
    logf('Cohort: %d cancers total, %d with >= %d paired patients', nrow(tab), npass, min_paired)

    if (isTRUE(opts$summarize)) {
        # append per-tier testable counts from each passing cancer's manifests
        addcols <- c('tier1_features', 'tier1_genes')
        vals <- matrix(NA_real_, nrow(tab), length(addcols), dimnames = list(NULL, addcols))
        for (i in which(tab$pass)) {
            slug <- tab$slug[i]
            s1 <- manifest_summary(file.path(opts[['output-dir']], slug, 'tier1',
                                             sprintf('%s_tier1_manifest.json', slug)))
            vals[i, 'tier1_features'] <- s1[['n_features']]
            vals[i, 'tier1_genes']    <- s1[['n_genes']]
        }
        for (k in addcols) tab[, (k) := vals[, k]]
        n_done <- sum(!is.na(vals[tab$pass, 'tier1_features']))
        logf('Summarize: found tier-1 manifests for %d / %d passing cancers', n_done, npass)
    }

    # cohort_manifest.tsv: one row per cancer, `pass` marks the analysed set
    manifest_path <- file.path(cohort_dir, 'cohort_manifest.tsv')
    fwrite(tab, manifest_path, sep = '\t')
    logf('Wrote %s', manifest_path)

    # cohort_summary.md
    show <- tab[, .(cancer = detailed_category, slug, n_paired,
                    tumor = n_tumor_samples, normal = n_normal_samples,
                    pass = ifelse(pass, 'yes', 'no'))]
    if (isTRUE(opts$summarize))
        show <- cbind(show, tab[, .(t1_features = tier1_features, t1_genes = tier1_genes)])
    md <- paste(c(
        '# Cohort summary',
        '',
        sprintf('Generated %s · min-paired = %d · qualifying cancers = %d / %d',
                format(Sys.time(), '%Y-%m-%d %H:%M:%S'), min_paired, npass, nrow(tab)),
        '',
        md_table(as.data.frame(show)),
        ''), collapse = '\n')
    writeLines(md, file.path(cohort_dir, 'cohort_summary.md'))
    logf('Wrote %s', file.path(cohort_dir, 'cohort_summary.md'))

    # cohort_plots.pdf (base graphics, no ggplot dependency)
    pdf(file.path(cohort_dir, 'cohort_plots.pdf'), width = 11, height = 7)
    on.exit(grDevices::dev.off(), add = TRUE)
    pass <- tab[pass == TRUE]
    par(mar = c(12, 4, 3, 1))
    cols <- ifelse(tab$pass, '#6C8EBF', '#cccccc')
    bp <- barplot(tab$n_paired, names.arg = tab$slug, las = 2, col = cols, border = NA, cex.names = 0.7,
                  main = sprintf('Paired patients per cancer (line = min %d)', min_paired),
                  ylab = 'paired patients')
    abline(h = min_paired, lty = 2, col = 'red'); text(bp, tab$n_paired, tab$n_paired, pos = 3, cex = 0.6)

    if (nrow(pass)) {
        m <- t(as.matrix(pass[, .(tumor = n_paired, normal = n_paired)]))   # 1 T + 1 N per paired patient
        barplot(m, names.arg = pass$slug, beside = TRUE, las = 2, cex.names = 0.7,
                col = c('#B85450', '#82B366'), border = NA,
                main = 'Tumor / Normal samples used (qualifying cancers)', ylab = 'samples')
        legend('topright', c('Tumor', 'Normal'), fill = c('#B85450', '#82B366'), border = NA, bty = 'n')
    }

    if (isTRUE(opts$summarize) && nrow(pass)) {
        if (any(!is.na(pass$tier1_features))) {
            barplot(pass$tier1_features, names.arg = pass$slug, las = 2, cex.names = 0.7,
                    col = '#9673A6', border = NA, main = 'Tier-1 testable features (gene x bucket)',
                    ylab = 'features')
            barplot(pass$tier1_genes, names.arg = pass$slug, las = 2, cex.names = 0.7,
                    col = '#6C8EBF', border = NA, main = 'Tier-1 testable genes (>=2 buckets)',
                    ylab = 'genes')
        }
    }
    logf('Wrote %s', file.path(cohort_dir, 'cohort_plots.pdf'))
    quit(save = 'no', status = 0)
}

tryCatch(main(), error = function(e) {
    logf('ERROR: %s', conditionMessage(e)); quit(save = 'no', status = 1)
})
