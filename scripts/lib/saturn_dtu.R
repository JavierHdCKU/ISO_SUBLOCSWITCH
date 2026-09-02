#!/usr/bin/env Rscript
# =============================================================================
# saturn_dtu.R: DTU testing engine, shared library
# -----------------------------------------------------------------------------
# Sourced by run_saturn_dtu_tier1.R. Defines functions only. Fits and tests
# gene x localization-bucket features with satuRn under a paired block design,
# and adds a model-free side log2FC descriptor. Switch calling is a separate
# step (filter_switch_hits.R).
#
# Statistical model (per feature i of gene g, sample s). satuRn's quasi-binomial
# GLM on the SHARE of the gene total (Gilis et al. 2022, F1000Research):
#
#     logit( p_{i,s} ) = log L_s + beta_0 + beta_{patient[s]} + beta_cond * 1{T}
#
# where p_{i,s} is feature i's share of its gene's reads in sample s and
# beta_cond the Tumor-vs-Normal effect on the logit-share scale. The design
# `~ 0 + condition + patient_id` is a randomized complete block design:
# `condition` is the main effect and `patient_id` the block absorbing
# patient-to-patient baseline variability. satuRn is driven through its own GLM
# interface rather than through IsoformSwitchAnalyzeR's importRdata, which
# injects SVA surrogate variables that on a paired tumour/normal design
# correlate with `condition` and cancel the signal.
#
# Dependencies: satuRn, SummarizedExperiment, BiocParallel.
# =============================================================================

SATURN_DTU_VERSION <- "1.0.0"

# Statistics carried into the result table. Significance is assessed on
# regular_FDR, the Benjamini-Hochberg adjustment of the raw p-value.
REPORTED_STATS <- c('estimates', 'se', 'df', 't', 'pval', 'regular_FDR')

suppressPackageStartupMessages({
    library(satuRn)
    library(SummarizedExperiment)
    library(BiocParallel)
})

# -----------------------------------------------------------------------------
# saturn_std(): fit and test one SummarizedExperiment with the paired
#                satuRn machinery. Returns native satuRn columns plus manual BH.
# -----------------------------------------------------------------------------
#' @param se      SummarizedExperiment. Requires colData columns `condition`
#'                (factor, levels c('Normal','Tumor')) and `patient_id`. rowData
#'                must carry `gene_id` so satuRn knows the binomial denominator
#'                (the per-gene total). Counts assay may be fractional expected
#'                counts (satuRn's quasi-binomial GLM does not require integers).
#' @param verbose passed to fitDTU; also reports the degenerate-fit drop count.
#' @return list(res, se, n_dropped):
#'   * res: data.frame keyed by feature (rownames = feature_id) with satuRn's
#'     estimates, se, df, t, pval and regular_FDR, plus feature_id and q_BH.
#'     Significance is assessed on `regular_FDR`, the Benjamini-Hochberg
#'     adjustment of the raw p-value. `q_BH` recomputes that adjustment with
#'     p.adjust(pval, 'BH') as a cross-check that the two agree.
#'   * se: the fitted and tested SE (rowData[["fitDTUResult_TvsN"]] in sync,
#'     full length), so callers can read per-feature dispersions.
#'   * n_dropped: number of degenerate per-feature fits dropped before testDTU.
saturn_std <- function(se, verbose = FALSE) {
    stopifnot('condition' %in% colnames(colData(se)), 'patient_id' %in% colnames(colData(se)))
    grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)            # satuRn draws during testDTU
    se <- satuRn::fitDTU(object = se, formula = ~ 0 + condition + patient_id,
                         parallel = FALSE, BPPARAM = SerialParam(), verbose = verbose)
    design <- model.matrix(~ 0 + condition + patient_id, data = as.data.frame(colData(se)))
    # Contrast = conditionTumor - conditionNormal on the logit-share scale.
    L <- matrix(0, ncol(design), 1, dimnames = list(colnames(design), "TvsN"))
    L["conditionTumor", 1] <- 1; L["conditionNormal", 1] <- -1

    # Drop degenerate per-feature fits before testDTU: a model with a NULL
    # vcovUnsc makes satuRn's varContrast test a zero-length condition and abort.
    # Keep models with a conformable vcov (the rest would return NA anyway) and
    # re-pad the result to full length.
    models <- rowData(se)[["fitDTUModels"]]
    ncoef  <- nrow(L)
    good   <- if (is.null(models)) rep(TRUE, nrow(se)) else
              vapply(as.list(models), function(m) {
                  vc <- tryCatch(m@params$vcovUnsc, error = function(e) NULL)
                  is.matrix(vc) && nrow(vc) == ncoef
              }, logical(1))
    n_drop <- sum(!good)
    if (verbose && n_drop) cat(sprintf('  dropped %d/%d degenerate fits before testDTU\n', n_drop, length(good)))

    if (!any(good)) {                                              # all fits degenerate -> all-NA result
        res <- data.frame(estimates = NA_real_, se = NA_real_, df = NA_real_, t = NA_real_,
                          pval = NA_real_, regular_FDR = NA_real_,
                          row.names = rownames(se))
    } else {
        se_ok  <- satuRn::testDTU(object = se[good, ], contrasts = L,
                                  sort = FALSE, diagplot1 = FALSE, diagplot2 = FALSE)
        res_ok <- as.data.frame(rowData(se_ok)[["fitDTUResult_TvsN"]])
        res    <- res_ok[match(rownames(se), rownames(res_ok)), , drop = FALSE]   # NA-pad dropped rows
        rownames(res) <- rownames(se)
    }
    rowData(se)[["fitDTUResult_TvsN"]] <- res                      # keep se in sync (full length)
    res <- res[, intersect(colnames(res), REPORTED_STATS), drop = FALSE]
    res$feature_id <- rownames(res)
    res$q_BH       <- p.adjust(res$pval, "BH")
    list(res = res, se = se, n_dropped = n_drop)
}

# -----------------------------------------------------------------------------
# compute_side_log2FC(): model-free direction and magnitude on raw counts
# -----------------------------------------------------------------------------
# satuRn's beta_cond is on the logit-share scale: it says the share moved, not
# which way the absolute count mass went. The side log2FC of the raw feature
# counts supplies that direction and magnitude:
#
#     side_log2FC(i) = log2( ( mean_T count(i) + 1 ) / ( mean_N count(i) + 1 ) )
#
# The +1 pseudocount keeps it finite when a feature is silent in one condition.
# count_mat is the gene x bucket count matrix.
#'
#' @param count_mat numeric matrix, features x samples (raw counts).
#' @param samples_T,samples_N character sample-id vectors (Tumor / Normal).
#' @return named numeric vector keyed by rownames(count_mat).
compute_side_log2FC <- function(count_mat, samples_T, samples_N) {
    sT <- intersect(samples_T, colnames(count_mat))
    sN <- intersect(samples_N, colnames(count_mat))
    stopifnot(length(sT) > 0, length(sN) > 0)
    val <- log2((rowMeans(count_mat[, sT, drop = FALSE]) + 1) /
                (rowMeans(count_mat[, sN, drop = FALSE]) + 1))
    setNames(as.numeric(val), rownames(count_mat))
}
