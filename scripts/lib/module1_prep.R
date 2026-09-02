#!/usr/bin/env Rscript
# =============================================================================
# module1_prep.R: data preparation, shared library
# -----------------------------------------------------------------------------
# Sourced by make_tier1_input.R. Defines constants and functions only.
#
# Pipeline (per cancer):
#   recover raw counts  2^x - 1   (expected counts; fractional, no rounding)
#     -> pair samples    1 Primary Tumor + 1 Solid Tissue Normal per patient
#     -> expression filter (lower 95% CI of the mean > 10 reads)
#     -> DeepLoc top-ranked compartment -> bucket {TM, Sec, Int}
#     -> SummarizedExperiment (gene x bucket, disjoint sums, >=2 buckets/gene)
#
# Dependencies: data.table; SummarizedExperiment at SE-build time only
# (referenced as SummarizedExperiment:: so the metadata-only cohort step can
#  source this file without the Bioconductor stack); base `tools` for md5.
# =============================================================================

MODULE1_VERSION <- "1.2.0"

suppressPackageStartupMessages({
    library(data.table)
})

# ── Constants ────────────────────────────────────────────────────────────────
# DeepLoc's 10 compartment-probability columns. Order is immaterial to argmax,
# but is fixed here so the assignment is reproducible.
COMPS_ALL <- c('Cytoplasm','Nucleus','Extracellular','Cell membrane','Mitochondrion',
               'Plastid','Endoplasmic reticulum','Lysosome/Vacuole','Golgi apparatus','Peroxisome')

# Broad-bucket map: TM = the surface membrane, Sec = secreted, Int = the 8
# intracellular compartments lumped together.
COMP_TO_BUCKET <- c('Cell membrane'='TM', 'Extracellular'='Sec',
                    'Cytoplasm'='Int','Nucleus'='Int','Mitochondrion'='Int','Plastid'='Int',
                    'Lysosome/Vacuole'='Int','Peroxisome'='Int','Endoplasmic reticulum'='Int',
                    'Golgi apparatus'='Int')

BROAD_CATS <- c('TM', 'Sec', 'Int')   # bucket names / ordering used by the builders

# ── Tiny utilities (dependency-free; usable by entrypoints + qc_report.R) ─────

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# stderr logger with timestamp (per repo convention: logs to stderr, stdout is data)
logf <- function(fmt, ...) {
    cat(sprintf('[%s] %s\n', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), sprintf(fmt, ...)),
        file = stderr())
}

# Slug for filesystem paths: "Breast Invasive Carcinoma" -> "Breast_Invasive_Carcinoma"
slugify <- function(s) {
    s <- gsub('[^A-Za-z0-9]+', '_', s)
    gsub('^_+|_+$', '', s)
}

# Directory of the running Rscript (so --input-dir / --output-dir defaults are
# robust to the caller's cwd).
get_script_dir <- function() {
    ca <- commandArgs(FALSE)
    m  <- grepl('^--file=', ca)
    if (any(m)) return(normalizePath(dirname(sub('^--file=', '', ca[which(m)[1]]))))
    normalizePath('.')
}

# Minimal dependency-free CLI parser. `spec` is a named list of
# list(default=, help=, flag=FALSE). Handles --key value, --flag, --version, --help.
m1_parse_cli <- function(argv, spec, prog = 'script.R') {
    opts <- lapply(spec, function(s) s$default)
    print_help <- function() {
        cat(sprintf('Usage: Rscript %s [options]\n\nOptions:\n', prog))
        for (k in names(spec)) {
            d  <- spec[[k]]$default
            ds <- if (is.null(d)) '' else sprintf('  [default: %s]', as.character(d))
            cat(sprintf('  --%-14s %s%s\n', k, spec[[k]]$help %||% '', ds))
        }
        cat('  --version        Print version and exit\n')
        cat('  --help           Show this help and exit\n')
    }
    i <- 1L
    while (i <= length(argv)) {
        a <- argv[i]
        if (a == '--version') { cat(MODULE1_VERSION, '\n'); quit(save = 'no', status = 0) }
        if (a %in% c('--help', '-h')) { print_help(); quit(save = 'no', status = 0) }
        if (!startsWith(a, '--')) stop(sprintf('Unexpected argument: %s', a))
        key <- sub('^--', '', a)
        if (!key %in% names(spec)) stop(sprintf('Unknown option: --%s', key))
        if (isTRUE(spec[[key]]$flag)) { opts[[key]] <- TRUE; i <- i + 1L }
        else {
            if (i + 1L > length(argv)) stop(sprintf('Option --%s requires a value', key))
            opts[[key]] <- argv[i + 1L]; i <- i + 2L
        }
    }
    opts
}

# ── Minimal dependency-free JSON writer (manifests) ──────────────────────────
# Avoids a hard jsonlite dependency.
# Named lists -> JSON objects; unnamed lists / length>1 atomics -> JSON arrays.
m1_json_escape <- function(s) {
    s <- gsub('\\\\', '\\\\\\\\', s)
    s <- gsub('"',  '\\\\"',  s)
    s <- gsub('\n', '\\\\n',  s)
    s <- gsub('\t', '\\\\t',  s)
    s <- gsub('\r', '\\\\r',  s)
    s
}
m1_to_json <- function(x, indent = 0L) {
    pad  <- strrep('  ', indent)
    pad2 <- strrep('  ', indent + 1L)
    if (is.null(x)) return('null')
    if (is.list(x)) {
        if (length(x) == 0L) return(if (!is.null(names(x))) '{}' else '[]')
        nm <- names(x)
        if (!is.null(nm) && all(nzchar(nm))) {                          # object
            items <- vapply(seq_along(x), function(i)
                sprintf('%s"%s": %s', pad2, m1_json_escape(nm[i]),
                        m1_to_json(x[[i]], indent + 1L)), character(1))
            return(sprintf('{\n%s\n%s}', paste(items, collapse = ',\n'), pad))
        }
        items <- vapply(x, function(e) sprintf('%s%s', pad2, m1_to_json(e, indent + 1L)),
                        character(1))                                   # array of elements
        return(sprintf('[\n%s\n%s]', paste(items, collapse = ',\n'), pad))
    }
    if (length(x) == 0L) return('null')
    if (length(x) > 1L) {                                              # array of scalars
        items <- vapply(x, function(e) m1_to_json(e, indent + 1L), character(1))
        return(sprintf('[%s]', paste(items, collapse = ', ')))
    }
    if (is.logical(x)) return(if (is.na(x)) 'null' else tolower(as.character(x)))
    if (is.numeric(x)) return(if (is.na(x)) 'null' else trimws(formatC(x, format = 'g', digits = 15, width = 0)))
    if (is.na(x))      return('null')
    sprintf('"%s"', m1_json_escape(as.character(x)))
}
write_json <- function(x, path) writeLines(m1_to_json(x), path)

# Provenance record for one input file: the file name that belongs in the slot,
# plus size and md5 to identify which version of it was used. md5 is opt-in via
# --full-md5 for the ~1 GB count matrices.
file_provenance <- function(path, md5 = FALSE) {
    real <- normalizePath(path, mustWork = TRUE)
    info <- file.info(real)
    out  <- list(file = basename(path),
                 size_bytes = as.numeric(info$size))
    if (md5) out$md5 <- unname(tools::md5sum(real))
    out
}

# Resolve and validate the 5 input files under input_dir.
module1_input_files <- function(input_dir) {
    f <- list(
        tumor_log2  = file.path(input_dir, 'tumor_log2counts.tsv.gz'),
        normal_log2 = file.path(input_dir, 'normal_log2counts.tsv.gz'),
        metadata    = file.path(input_dir, 'metadata.tsv.gz'),
        tx2gene     = file.path(input_dir, 'tx2gene.tsv.gz'),
        deeploc     = file.path(input_dir, 'deeploc.csv')
    )
    miss <- !file.exists(unlist(f))
    if (any(miss))
        stop(sprintf('Missing input file(s) under %s:\n  %s',
                     input_dir, paste(unlist(f)[miss], collapse = '\n  ')))
    f
}

# =============================================================================
# PREP STEPS  (each is independently testable; prepare_module1_input strings
#              them together in order)
# =============================================================================

# (a) Paired samples: 1 Primary-Tumor + 1 Solid-Tissue-Normal per patient.
#     Deterministic: sort by sample id before the `.SD[1]` pick so the chosen
#     vial is reproducible across reruns and machines.
select_paired_samples <- function(metadata_path, cancer, seed = 42L, n_paired = NULL) {
    meta <- fread(metadata_path)
    if (!'detailed_category' %in% colnames(meta))
        stop("metadata is missing the 'detailed_category' column")
    mb <- meta[detailed_category == cancer]
    mb <- mb[`_sample_type` %in% c('Primary Tumor', 'Solid Tissue Normal')]   # T / N only (defensive)
    if (nrow(mb) == 0)
        stop(sprintf("No Primary Tumor / Solid Tissue Normal samples for cancer '%s'", cancer))
    mb[, patient_id := substr(sample, 1, 12)]
    paired <- mb[, .(nt = sum(`_sample_type` == 'Primary Tumor'),
                     nn = sum(`_sample_type` == 'Solid Tissue Normal')),
                 by = patient_id][nt >= 1 & nn >= 1, patient_id]
    if (length(paired) == 0)
        stop(sprintf("No paired (Tumor+Normal) patients for cancer '%s'", cancer))
    if (!is.null(n_paired) && length(paired) > n_paired) {           # optional subset (testing)
        set.seed(seed); paired <- sample(paired, n_paired)
    }
    mb <- mb[patient_id %in% paired]
    setorder(mb, patient_id, `_sample_type`, sample)                 # deterministic vial pick
    mb <- mb[, .SD[1], by = .(patient_id, `_sample_type`)]           # 1 T + 1 N / patient
    mb[, condition := factor(ifelse(`_sample_type` == 'Primary Tumor', 'Tumor', 'Normal'),
                             levels = c('Normal', 'Tumor'))]
    mb[]
}

# (b) Recover RAW counts for all multi-iso transcripts; assemble count_mat,
#     col_meta, and the (subset) tx2gene table (with chr if present, for QC).
load_raw_counts <- function(tumor_path, normal_path, tx2gene_path, mb, multi_iso = TRUE) {
    t2g_full <- fread(tx2gene_path)
    if (!all(c('transcript_id', 'gene_id') %in% colnames(t2g_full)))
        stop("tx2gene must contain 'transcript_id' and 'gene_id' columns")
    keep_cols <- c('transcript_id', 'gene_id', intersect('chr', colnames(t2g_full)))
    t2g <- t2g_full[, ..keep_cols]

    if (multi_iso) {                                                 # genes with >=2 transcripts
        multi_genes <- t2g[, .N, by = gene_id][N >= 2, gene_id]
        keep_tx0    <- t2g[gene_id %in% multi_genes, transcript_id]
    } else {
        keep_tx0 <- t2g$transcript_id
    }

    tum_samples <- mb[condition == 'Tumor',  sample]
    nor_samples <- mb[condition == 'Normal', sample]

    # Read the header alone, then pull ONLY this cohort's columns. The tumour
    # matrix holds every TCGA sample (~6.4k) while one cancer needs a few dozen,
    # so reading it whole costs ~4 GB of RAM and yields nothing extra. Column
    # order still follows the file, so the result is identical either way.
    read_rec <- function(path, samples) {
        hdr  <- names(fread(path, nrows = 0L))
        cols <- intersect(hdr[-1], samples)
        if (length(cols) == 0)
            stop(sprintf("No sample columns matched in %s.\n  file cols e.g.: %s\n  wanted e.g.:    %s",
                         basename(path),
                         paste(head(hdr[-1], 2), collapse = ', '),
                         paste(head(samples, 2), collapse = ', ')))
        dt <- fread(path, select = c(hdr[1], cols))
        setnames(dt, 1, 'tx'); dt[, tx := sub('\\.[0-9]+$', '', tx)]
        dt <- dt[tx %in% keep_tx0]
        x <- as.matrix(dt[, ..cols]); storage.mode(x) <- 'double'
        x <- 2^x - 1                                                # log2 -> RAW expected counts (fractional; non-negative by construction)
        matrix(x, nrow = nrow(dt), ncol = length(cols), dimnames = list(dt$tx, cols))
    }
    mt  <- read_rec(tumor_path,  tum_samples)
    mn  <- read_rec(normal_path, nor_samples)
    ctx <- intersect(rownames(mt), rownames(mn))
    if (length(ctx) == 0)
        stop("No transcripts shared between tumor and normal matrices after count recovery.")
    count_mat <- cbind(mt[ctx, , drop = FALSE], mn[ctx, , drop = FALSE])
    storage.mode(count_mat) <- 'double'   # fractional expected counts (satuRn quasi-likelihood accepts non-integer)

    col_meta <- as.data.frame(mb[match(colnames(count_mat), sample),
                                 .(sample, patient_id, condition)])
    col_meta$patient_id <- factor(col_meta$patient_id)
    rownames(col_meta) <- col_meta$sample
    stopifnot(all(colnames(count_mat) == col_meta$sample))

    t2g <- t2g[transcript_id %in% rownames(count_mat)]
    list(count_mat = count_mat, col_meta = col_meta, t2g = t2g)
}

# (c) Expression filter on raw counts. Lower 95% CI of the mean per condition:
#     CI_lo = mu - 1.96 * sd / sqrt(n), with
#     sd = sqrt( (Sum x^2 - n*mu^2) / (n-1) ).
isar_count_filter <- function(count_mat, g, cond, ci_threshold = 10) {
    lci <- function(m) sapply(c('Normal', 'Tumor'), function(L) {
        # double storage is required: x*x on integer counts overflows int32.
        x   <- m[, cond == L, drop = FALSE]; storage.mode(x) <- 'double'
        n <- ncol(x); mu <- rowMeans(x)
        sdr <- sqrt(pmax(0, (rowSums(x * x) - n * mu * mu) / (n - 1)))
        mu - 1.96 * sdr / sqrt(n)
    })
    g_lci <- lci(rowsum(count_mat, g))                              # gene-level CI
    t_lci <- lci(count_mat)                                         # transcript-level CI
    f1 <- names(table(g))[table(g) >= 2]                            # (1) gene has >=2 transcripts
    f2 <- rownames(g_lci)[ g_lci[, 'Normal'] > ci_threshold &
                           g_lci[, 'Tumor']  > ci_threshold ]       # (2) gene active in BOTH N AND T
    f3 <- unique(g[ rownames(t_lci)[ t_lci[, 'Normal'] > ci_threshold |
                                     t_lci[, 'Tumor']  > ci_threshold ] ])  # (3) >=1 tx active in EITHER
    keep_genes <- Reduce(intersect, list(f1, f2, f3))
    list(keep_genes = keep_genes, f1 = f1, f2 = f2, f3 = f3, g_lci = g_lci, t_lci = t_lci)
}

# (d) DeepLoc P (probabilities) -> argmax compartment -> broad bucket.
#     Uses the model's own ranking (no probability threshold, no tie machinery).
#     TM is restricted to cell-surface transcripts. A transcript is TM iff
#       (1) its argmax localization is 'Cell membrane'  AND
#       (2) DeepLoc's 'Membrane types' field reports a transmembrane helix.
#     Cell-membrane proteins that are Peripheral / Lipid-anchored / Soluble (no TM
#     helix: GPI/prenyl/peripheral) carry no transmembrane domain and are reassigned
#     to Int. Both conditions are needed: transmembrane proteins of intracellular
#     membranes (ER, Golgi, mitochondrion, lysosome) fall to Int through their
#     top-ranked compartment, and membrane type alone would pull them into TM.
assign_deeploc_buckets <- function(deeploc_path, tx_ids, tm_require_transmembrane = TRUE) {
    dl <- fread(deeploc_path)
    if (!'Protein_ID' %in% colnames(dl))
        stop("deeploc file is missing the 'Protein_ID' column")
    dl <- dl[Protein_ID %in% tx_ids]
    cols <- COMPS_ALL
    if (!all(cols %in% colnames(dl)))
        stop(sprintf('deeploc file missing compartment column(s): %s',
                     paste(setdiff(cols, colnames(dl)), collapse = ', ')))
    P <- as.matrix(dl[, ..cols]); rownames(P) <- dl$Protein_ID
    stopifnot(setequal(colnames(P), names(COMP_TO_BUCKET)))
    top_comp  <- setNames(colnames(P)[max.col(P, ties.method = 'first')], rownames(P))
    tx_bucket <- setNames(unname(COMP_TO_BUCKET[top_comp]), rownames(P))

    # ── strict TM: demote Cell-membrane tx with no predicted TM helix to Int ──
    memtype   <- setNames(rep(NA_character_, nrow(dl)), dl$Protein_ID)[names(tx_bucket)]
    n_demoted <- 0L
    if (isTRUE(tm_require_transmembrane)) {
        if (!'Membrane types' %in% colnames(dl))
            stop("deeploc file is missing the 'Membrane types' column (required for strict TM)")
        memtype <- setNames(as.character(dl[['Membrane types']]), dl$Protein_ID)[names(tx_bucket)]
        demote  <- tx_bucket == 'TM' & !grepl('Transmembrane', memtype, fixed = TRUE)
        demote[is.na(demote)] <- TRUE                  # no membrane-type call -> not a surface TM
        n_demoted <- sum(demote)
        tx_bucket[demote] <- 'Int'
    }
    stopifnot(!anyNA(tx_bucket))
    list(P = P, top_comp = top_comp, tx_bucket = tx_bucket,
         memtype = memtype, n_tm_demoted = n_demoted,
         tm_require_transmembrane = isTRUE(tm_require_transmembrane))
}

# ── Orchestrator: shared prep bundle, steps (a)-(e) end-to-end ───────────────
prepare_module1_input <- function(input_dir, cancer,
                                   ci_threshold = 10, multi_iso = TRUE,
                                   tm_require_transmembrane = TRUE,
                                   seed = 42L, n_paired = NULL,
                                   cancer_slug = NULL, verbose = TRUE) {
    files <- module1_input_files(input_dir)
    if (is.null(cancer_slug)) cancer_slug <- slugify(cancer)

    ## (a) paired samples
    mb <- select_paired_samples(files$metadata, cancer, seed = seed, n_paired = n_paired)
    tum_samples <- mb[condition == 'Tumor',  sample]
    nor_samples <- mb[condition == 'Normal', sample]
    if (verbose) logf('Cohort %s [%s]: %d paired patients (%d Tumor + %d Normal samples)',
                      cancer, cancer_slug, length(unique(mb$patient_id)),
                      length(tum_samples), length(nor_samples))

    ## (b) raw counts + tx2gene + col_meta
    rc <- load_raw_counts(files$tumor_log2, files$normal_log2, files$tx2gene, mb, multi_iso = multi_iso)
    count_mat <- rc$count_mat; t2g <- rc$t2g; col_meta <- rc$col_meta
    cond <- as.character(col_meta$condition)
    g <- setNames(t2g$gene_id, t2g$transcript_id)[rownames(count_mat)]
    if (verbose) logf('Recovered raw counts: %d transcripts x %d samples',
                      nrow(count_mat), ncol(count_mat))

    ## (c) expression filter on raw counts
    isar <- isar_count_filter(count_mat, g, cond, ci_threshold = ci_threshold)
    keep_genes <- isar$keep_genes
    kept_tx    <- names(g)[g %in% keep_genes]
    if (verbose) logf('Expression filter: f1=%d f2=%d f3=%d -> kept %d genes / %d transcripts',
                      length(isar$f1), length(isar$f2), length(isar$f3),
                      length(keep_genes), length(kept_tx))

    ## (d) DeepLoc compartment -> bucket
    dl <- assign_deeploc_buckets(files$deeploc, rownames(count_mat),
                                 tm_require_transmembrane = tm_require_transmembrane)
    tx_bucket <- dl$tx_bucket; top_comp <- dl$top_comp; P <- dl$P
    if (verbose) logf('%d Cell-membrane transcripts without a TM helix -> Int',
                      dl$n_tm_demoted)

    ## (e) testable set: filter-passing and DeepLoc-annotated
    keep_tx_final <- intersect(kept_tx, names(tx_bucket))
    lib_sizes <- colSums(count_mat)
    if (verbose) logf('Testable: %d transcripts / %d genes',
                      length(keep_tx_final), length(unique(g[keep_tx_final])))

    funnel <- list(
        multi_iso_genes     = length(isar$f1),                      # genes with >=2 tx in count_mat
        f1_ge2tx_genes      = length(isar$f1),
        f2_bothNT_genes     = length(isar$f2),
        f3_eitherTx_genes   = length(isar$f3),
        isar_pass_genes     = length(keep_genes),
        isar_pass_tx        = length(kept_tx),
        deeploc_annot_genes = length(unique(g[keep_tx_final])),
        deeploc_annot_tx    = length(keep_tx_final),
        tm_demoted_to_int   = dl$n_tm_demoted               # Cell-membrane tx w/o TM helix -> Int
    )

    list(
        cancer = cancer, slug = cancer_slug, files = files,
        count_mat = count_mat, col_meta = col_meta, t2g = as.data.frame(t2g),
        g = g, lib_sizes = lib_sizes, cond = cond,
        tum_samples = tum_samples, nor_samples = nor_samples,
        tx_bucket = tx_bucket, top_comp = top_comp, P = P,
        keep_genes = keep_genes, kept_tx = kept_tx, keep_tx_final = keep_tx_final,
        n_pairs = length(intersect(
            col_meta$patient_id[col_meta$condition == 'Tumor'],
            col_meta$patient_id[col_meta$condition == 'Normal'])),
        funnel = funnel,
        thresholds = list(ci_threshold = ci_threshold, multi_iso = multi_iso,
                          tm_require_transmembrane = tm_require_transmembrane,
                          tm_definition = if (isTRUE(tm_require_transmembrane))
                              'cell_membrane_AND_transmembrane_helix' else 'cell_membrane_argmax',
                          seed = seed,
                          n_paired = if (is.null(n_paired)) 'all' else n_paired)
    )
}

# =============================================================================
# SE builder: gene x localization-bucket SummarizedExperiment.
# =============================================================================

# ── Builder: gene x bucket SE, keeps genes with >=2 buckets ───────
build_tier1_long_SE <- function(count_mat, tx_bucket, t2g, col_data,
                                 sample_col = 'sample', verbose = TRUE) {

    stopifnot(all(c('transcript_id', 'gene_id') %in% colnames(t2g)))
    stopifnot(sample_col %in% colnames(col_data))

    # Align tx_bucket to count_mat rows (drop tx with no DeepLoc / no bucket)
    tx_bucket <- tx_bucket[rownames(count_mat)]
    has_bucket <- !is.na(tx_bucket)
    if (verbose) cat(sprintf('Tier 1 build: %d / %d tx with a bucket assignment\n',
                             sum(has_bucket), length(tx_bucket)))

    # tx -> gene mapping for tx with a bucket
    tx2g <- setNames(t2g$gene_id, t2g$transcript_id)
    tx2g <- tx2g[names(tx_bucket)[has_bucket]]
    genes <- unique(tx2g)

    if (verbose) cat(sprintf('Building Tier 1 SE: %d genes x up to %d buckets x %d samples\n',
                             length(genes), length(BROAD_CATS), ncol(count_mat)))

    out_rows  <- list()
    out_names <- character()
    out_gid   <- character()
    out_buck  <- character()
    progress_step <- max(1, length(genes) %/% 20)
    for (gi in seq_along(genes)) {
        g <- genes[gi]
        tx_g <- names(tx2g)[tx2g == g]
        for (B in BROAD_CATS) {
            tx_gB <- intersect(tx_g, names(tx_bucket)[tx_bucket == B & has_bucket])
            if (length(tx_gB) == 0) next        # empty bucket for this gene -> skip
            key <- sprintf('%s__%s', g, B)
            out_rows[[key]] <- colSums(count_mat[tx_gB, , drop = FALSE])
            out_names <- c(out_names, key)
            out_gid   <- c(out_gid,   g)
            out_buck  <- c(out_buck,  B)
        }
        if (verbose && gi %% progress_step == 0)
            cat(sprintf('  ... %d / %d genes aggregated\n', gi, length(genes)))
    }

    tier1_mat <- do.call(rbind, out_rows)
    rownames(tier1_mat) <- out_names

    row_meta <- data.frame(
        isoform_id = out_names,
        gene_id    = out_gid,
        bucket     = out_buck,
        stringsAsFactors = FALSE
    )
    rownames(row_meta) <- out_names

    # Drop genes with < 2 non-empty buckets (no DTU testable)
    n_per_gene <- table(row_meta$gene_id)
    keep_genes <- names(n_per_gene[n_per_gene >= 2])
    keep <- row_meta$gene_id %in% keep_genes
    tier1_mat <- tier1_mat[keep, , drop = FALSE]
    row_meta  <- row_meta[keep, , drop = FALSE]
    if (verbose) cat(sprintf('After multi-bucket filter: %d rows across %d genes\n',
                             nrow(tier1_mat), length(keep_genes)))

    col_df <- as.data.frame(col_data)
    col_df <- col_df[match(colnames(tier1_mat), col_df[[sample_col]]), , drop = FALSE]
    rownames(col_df) <- col_df[[sample_col]]

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(counts = tier1_mat),
        colData = col_df,
        rowData = row_meta
    )

    # Mass-conservation sanity (disjoint sums should be EXACT modulo dropped tx)
    g0 <- keep_genes[1]
    tx0 <- intersect(names(tx2g)[tx2g == g0], names(tx_bucket)[has_bucket])
    expected <- colSums(count_mat[tx0, , drop = FALSE])
    bucket_rows <- paste(g0, BROAD_CATS, sep = '__')
    bucket_rows <- bucket_rows[bucket_rows %in% rownames(tier1_mat)]
    got <- colSums(tier1_mat[bucket_rows, , drop = FALSE])
    drift <- max(abs(expected - got))
    # counts are fractional doubles, so drift is a float: %g, never %d.
    if (verbose) cat(sprintf('Mass-conservation check on %s: max |drift| = %.3g (expect ~0)\n',
                             g0, drift))

    se
}

