#!/usr/bin/env Rscript
# =============================================================================
# run_saturn_dtu_tier1.R: DTU testing entrypoint
# -----------------------------------------------------------------------------
# Standalone CLI. Reads the satuRn input (gene x localization-bucket SE) for one
# cancer, fits and tests it with satuRn under a paired block design, and writes
# the result table and p-value/FDR QC into the same directory. Switch calling is
# the next step, filter_switch_hits.R.
#
# Consumes:  data/outputs/<slug>/tier1/<slug>_tier1_se.rds
# Writes:
#   <slug>_tier1_dtu_results.{tsv,rds}   satuRn columns + gene_id + bucket +
#                                        q_BH + side_log2FC
#   <slug>_tier1_dtu_qc_pfdr.png         raw p | regular_FDR
#   <slug>_tier1_dtu_qc.md               diagnostics + yield metrics
#   <slug>_tier1_dtu_manifest.json       params, contrast, versions, lambda, ...
#
# Significance is reported on `regular_FDR`, the Benjamini-Hochberg adjustment
# of the raw satuRn p-value.
#
# Usage:
#   Rscript run_saturn_dtu_tier1.R --cancer-slug Breast_Invasive_Carcinoma
#   Rscript run_saturn_dtu_tier1.R --cancer "Stomach Adenocarcinoma"   # slugified
#   Rscript run_saturn_dtu_tier1.R --cancer-slug Kidney_Chromophobe --fdr 0.10
# =============================================================================

.script_dir <- local({
    ca <- commandArgs(FALSE); m <- grepl('^--file=', ca)
    if (any(m)) normalizePath(dirname(sub('^--file=', '', ca[which(m)[1]]))) else normalizePath('.')
})
source(file.path(.script_dir, 'lib', 'module1_prep.R'))   # data.table, slugify, logf, %||%, write_json
source(file.path(.script_dir, 'lib', 'qc_report.R'))       # md_table, fmt_num, pkg_versions, assert, hard_failed
source(file.path(.script_dir, 'lib', 'saturn_dtu.R'))      # saturn_std, compute_side_log2FC, SATURN_DTU_VERSION
source(file.path(.script_dir, 'lib', 'qc_results.R'))      # compute_qc_results, plot, finalize_dtu_results

.repo_root <- dirname(.script_dir)

# ── dependency-free CLI parser (prints SATURN_DTU_VERSION on --version) ─────
spec <- list(
    `cancer-slug` = list(default = NULL, help = 'filesystem slug = <output-dir>/<slug>/ (primary selector)'),
    `cancer`      = list(default = NULL, help = 'detailed_category string; slugified if --cancer-slug absent'),
    `output-dir`  = list(default = file.path(.repo_root, 'data', 'outputs'), help = 'root output directory'),
    `fdr`         = list(default = '0.05', help = 'FDR cutoff used for QC hit counts (NOT a gate)')
)
parse_cli <- function(argv) {
    opts <- lapply(spec, function(s) s$default)
    print_help <- function() {
        cat(sprintf('Usage: Rscript run_saturn_dtu_tier1.R [options]\n\nOptions:\n'))
        for (k in names(spec)) {
            d <- spec[[k]]$default; ds <- if (is.null(d)) '' else sprintf('  [default: %s]', as.character(d))
            cat(sprintf('  --%-12s %s%s\n', k, spec[[k]]$help, ds))
        }
        cat('  --version      Print version and exit\n  --help         Show this help and exit\n')
    }
    i <- 1L
    while (i <= length(argv)) {
        a <- argv[i]
        if (a == '--version') { cat(SATURN_DTU_VERSION, '\n'); quit(save = 'no', status = 0) }
        if (a %in% c('--help', '-h')) { print_help(); quit(save = 'no', status = 0) }
        if (!startsWith(a, '--')) stop(sprintf('Unexpected argument: %s', a))
        key <- sub('^--', '', a)
        if (!key %in% names(spec)) stop(sprintf('Unknown option: --%s', key))
        if (i + 1L > length(argv)) stop(sprintf('Option --%s requires a value', key))
        opts[[key]] <- argv[i + 1L]; i <- i + 2L
    }
    opts
}

# ── read a top-level string field from a manifest.json without jsonlite ─────
manifest_str_field <- function(path, key) {
    if (!file.exists(path)) return(NA_character_)
    txt <- paste(readLines(path, warn = FALSE), collapse = '\n')
    mm  <- regmatches(txt, regexpr(sprintf('"%s"[[:space:]]*:[[:space:]]*"[^"]*"', key), txt))
    if (!length(mm)) return(NA_character_)
    sub('.*:[[:space:]]*"([^"]*)"$', '\\1', mm)
}

# ── tumor / normal sample ids from an SE's colData ──────────────────────────
tn_samples <- function(se) {
    cd <- as.data.frame(colData(se))
    list(T = cd$sample[cd$condition == 'Tumor'], N = cd$sample[cd$condition == 'Normal'])
}

# ── column order for the result tables ──────────────────────────────────────
order_cols <- function(res, group_col) {
    stats  <- c('estimates', 'se', 'df', 't', 'pval', 'regular_FDR')
    front  <- c('feature_id', 'gene_id', group_col)
    res[, c(front, stats, 'q_BH', 'side_log2FC')]
}

# =============================================================================
# One gene x bucket SE
# =============================================================================
run_tier1 <- function(slug, cancer, out_root, fdr) {
    in_dir  <- file.path(out_root, slug, 'tier1')
    se_path <- file.path(in_dir, sprintf('%s_tier1_se.rds', slug))
    if (!file.exists(se_path)) stop(sprintf('SE not found: %s', se_path))
    logf('Tier-1: reading %s', se_path)
    se <- readRDS(se_path)

    fit <- saturn_std(se, verbose = TRUE)
    res <- fit$res
    rd  <- as.data.frame(rowData(se))                     # isoform_id, gene_id, bucket
    res$gene_id <- rd[res$feature_id, 'gene_id']
    res$bucket  <- rd[res$feature_id, 'bucket']
    tn  <- tn_samples(se)
    side <- compute_side_log2FC(as.matrix(assay(se, 'counts')), tn$T, tn$N)
    res$side_log2FC <- side[res$feature_id]
    res <- order_cols(res, 'bucket')

    meta <- list(cancer = cancer, slug = slug, tier_label = 'Tier-1',
                 unit_desc = 'gene x localization-bucket pseudo-features (`gene__bucket`); satuRn tests the bucket SHARE of the gene total.')
    finalize_dtu_results('tier1', meta, res, in_dir,
                         input_files = c(tier1_se = se_path), fdr = fdr, n_dropped = fit$n_dropped)
}

# =============================================================================
main <- function() {
    opts <- parse_cli(commandArgs(TRUE))
    if (is.null(opts[['cancer-slug']]) && is.null(opts$cancer))
        stop('Provide --cancer-slug <SLUG> (or --cancer "<detailed_category>")')
    slug <- if (!is.null(opts[['cancer-slug']])) opts[['cancer-slug']] else slugify(opts$cancer)
    out_root <- opts[['output-dir']]
    fdr  <- as.numeric(opts$fdr)

    cancer <- opts$cancer %||%
        manifest_str_field(file.path(out_root, slug, 'tier1', sprintf('%s_tier1_manifest.json', slug)), 'cancer')
    if (is.na(cancer) || is.null(cancer)) cancer <- slug

    logf('run_saturn_dtu_tier1.R v%s | slug="%s" cancer="%s" | fdr=%.3f | out=%s',
         SATURN_DTU_VERSION, slug, cancer, fdr, out_root)

    ok <- isTRUE(run_tier1(slug, cancer, out_root, fdr))
    if (!ok) { logf('Tier-1 DTU QC FAILED for %s, see *_tier1_dtu_qc.md', slug); quit(save = 'no', status = 2) }
    logf('Tier-1 DTU done for %s', slug)
    quit(save = 'no', status = 0)
}

tryCatch(main(), error = function(e) {
    logf('ERROR: %s', conditionMessage(e)); quit(save = 'no', status = 1)
})
