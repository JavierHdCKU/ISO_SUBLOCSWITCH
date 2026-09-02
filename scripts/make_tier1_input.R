#!/usr/bin/env Rscript
# =============================================================================
# make_tier1_input.R: build the gene x localization-bucket satuRn input
# -----------------------------------------------------------------------------
# Standalone CLI. Runs the shared data preparation, builds the SE (disjoint
# DeepLoc-bucket sums per gene; keeps genes with >=2 occupied buckets), writes
# <slug>_tier1_se.rds + manifest.json + qc.md under <output-dir>/<slug>/tier1/,
# and exits non-zero if a hard QC assertion fails.
#
# Usage:
#   Rscript make_tier1_input.R --cancer "Breast Invasive Carcinoma"
#   Rscript make_tier1_input.R --cancer "Lung Adenocarcinoma" \
#       --input-dir data/inputs --output-dir data/outputs
# =============================================================================

.script_dir <- local({
    ca <- commandArgs(FALSE); m <- grepl('^--file=', ca)
    if (any(m)) normalizePath(dirname(sub('^--file=', '', ca[which(m)[1]]))) else normalizePath('.')
})
source(file.path(.script_dir, 'lib', 'module1_prep.R'))
source(file.path(.script_dir, 'lib', 'qc_report.R'))

.repo_root <- dirname(.script_dir)
spec <- list(
    `cancer`       = list(default = NULL,  help = 'detailed_category string (required), e.g. "Lung Adenocarcinoma"'),
    `input-dir`    = list(default = file.path(.repo_root, 'data', 'inputs'),  help = 'directory with the 5 input files'),
    `output-dir`   = list(default = file.path(.repo_root, 'data', 'outputs'), help = 'root output directory'),
    `cancer-slug`  = list(default = NULL,  help = 'override the filesystem slug (default: slugified cancer)'),
    `ci-threshold` = list(default = '10',  help = 'ISAR lower-95%-CI-of-mean read threshold'),
    `seed`         = list(default = '42',  help = 'RNG seed (only used if --n-paired subsets)'),
    `n-paired`     = list(default = NULL,  help = 'subset to N paired patients (default: all)'),
    `full-md5`     = list(default = FALSE, help = 'also md5 the ~1 GB count matrices', flag = TRUE)
)

main <- function() {
    opts <- m1_parse_cli(commandArgs(TRUE), spec, prog = 'make_tier1_input.R')
    if (is.null(opts$cancer)) stop('--cancer is required (a detailed_category string)')

    seed     <- as.integer(opts$seed)
    n_paired <- if (is.null(opts[['n-paired']])) NULL else as.integer(opts[['n-paired']])
    set.seed(seed)

    logf('make_tier1_input.R v%s | cancer="%s" | input=%s | output=%s',
         MODULE1_VERSION, opts$cancer, opts[['input-dir']], opts[['output-dir']])

    prep <- prepare_module1_input(
        input_dir    = opts[['input-dir']],
        cancer       = opts$cancer,
        ci_threshold = as.numeric(opts[['ci-threshold']]),
        multi_iso    = TRUE,
        seed         = seed,
        n_paired     = n_paired,
        cancer_slug  = opts[['cancer-slug']],
        verbose      = TRUE)

    se <- build_tier1_long_SE(
        prep$count_mat[prep$keep_tx_final, , drop = FALSE],
        prep$tx_bucket[prep$keep_tx_final],
        prep$t2g, prep$col_meta, sample_col = 'sample', verbose = FALSE)

    out_dir <- file.path(opts[['output-dir']], prep$slug, 'tier1')
    ok <- finalize_tier1(prep, se, out_dir, full_md5 = isTRUE(opts[['full-md5']]))

    if (!ok) { logf('QC FAILED for %s - see qc.md', prep$slug); quit(save = 'no', status = 2) }
    logf('Tier-1 done for %s', prep$slug)
    quit(save = 'no', status = 0)
}

tryCatch(main(), error = function(e) {
    logf('ERROR: %s', conditionMessage(e)); quit(save = 'no', status = 1)
})
