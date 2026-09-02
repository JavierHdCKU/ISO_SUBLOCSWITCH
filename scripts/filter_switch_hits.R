#!/usr/bin/env Rscript
# =============================================================================
# filter_switch_hits.R: localization-switch calling across all cancers
# -----------------------------------------------------------------------------
# Applies the switch filter to every cancer with a satuRn DTU result and
# emits, per cancer, the gene-level src->dst switch calls; plus a combined table
# and a markdown report quantifying how many hits are lost at each filter step.
#
# Filter funnel (per gene x localization bucket = TM / Sec / Int):
#   0. testable        every bucket satuRn tested
#   1. significant     regular_FDR < 0.05  (Benjamini-Hochberg on the raw
#                      satuRn p-value)
#   2. valid effect    is.finite(side_log2FC) & side_log2FC != 0
#   3. expressed       mean_count >= 10 in both Normal and Tumor, applied per
#                      bucket (a near-empty compartment inflates side_log2FC
#                      through the +1 pseudocount)
#   4. >=2 buckets     a switch is a reallocation, so the gene needs >=2 buckets
#                      surviving steps 1-3
#   5. directional     dst = bucket with max side_log2FC (gained most share),
#                      src = min (lost most share); keep genes with src != dst
# switch_type records the expression background: cross-compartment (buckets move
# both ways at ~constant total), upregulated or downregulated and switch (all
# buckets move the same way). Direction: share_i = FC_i / FC_gene, so ranking
# buckets by count log2FC ranks them by share gained or lost.
#
# Usage:  Rscript filter_switch_hits.R [--fdr 0.05] [--min-count 10]
#                                      [--root DIR] [--tx2gene FILE]
# Outputs: <root>/<cancer>/tier1/<cancer>_tier1_filtered_hits.csv   (per cancer)
#          <root>/_cohort/tier1_filtered_hits_ALL.csv               (combined)
#          <root>/_cohort/tier1_filtering_report.md                 (funnel+why)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(SummarizedExperiment) })

.script_dir <- local({ ca <- commandArgs(FALSE); m <- grepl('^--file=', ca)
  if (any(m)) normalizePath(dirname(sub('^--file=', '', ca[which(m)[1]]))) else normalizePath('.') })
.repo <- dirname(.script_dir)
args <- commandArgs(TRUE)
opt  <- function(f, d) { i <- which(args == f); if (length(i)) args[i + 1] else d }
FDR   <- as.numeric(opt('--fdr', '0.05'))
MINC  <- as.numeric(opt('--min-count', '10'))
ROOT  <- normalizePath(opt('--root', file.path(.repo, 'data', 'outputs')))
TX2G  <- opt('--tx2gene', file.path(.repo, 'data', 'inputs', 'tx2gene.tsv.gz'))
COH   <- file.path(ROOT, '_cohort'); dir.create(COH, showWarnings = FALSE, recursive = TRUE)
logf  <- function(...) cat(sprintf(...), '\n')

# ---- cancers with a Tier-1 DTU result ---------------------------------------
slugs <- setdiff(list.dirs(ROOT, recursive = FALSE, full.names = FALSE),
                 c('_cohort', 'figures', 'hits_topology'))
slugs <- slugs[file.exists(file.path(ROOT, slugs, 'tier1', sprintf('%s_tier1_dtu_results.rds', slugs)))]
logf('%d cancers found', length(slugs))

# ---- load DTU results (gene x bucket) ---------------------------------------
t1 <- rbindlist(lapply(slugs, function(s)
  as.data.table(readRDS(file.path(ROOT, s, 'tier1', sprintf('%s_tier1_dtu_results.rds', s))))[, cancer := s]),
  use.names = TRUE, fill = TRUE)

# ---- per-bucket expression (mean count per condition) + dIF, one SE pass -----
bucket_stats <- function(se) {
  A <- as.matrix(assay(se, 'counts')); g <- as.data.frame(rowData(se))$gene_id
  cond <- as.character(colData(se)$condition)
  share <- A / rowsum(A, g)[g, , drop = FALSE]; share[!is.finite(share)] <- NA
  data.table(feature_id   = rownames(A),
             mean_count_N = rowMeans(A[, cond == 'Normal', drop = FALSE]),
             mean_count_T = rowMeans(A[, cond == 'Tumor',  drop = FALSE]),
             dIF          = rowMeans(share[, cond == 'Tumor',  drop = FALSE], na.rm = TRUE) -
                            rowMeans(share[, cond == 'Normal', drop = FALSE], na.rm = TRUE))
}
bs <- rbindlist(lapply(slugs, function(s)
  bucket_stats(readRDS(file.path(ROOT, s, 'tier1', sprintf('%s_tier1_se.rds', s))))[, cancer := s]))
t1 <- merge(t1, bs, by = c('cancer', 'feature_id'), all.x = TRUE)

# ---- staged funnel ----------------------------------------------------------
s1 <- t1[regular_FDR < FDR]                                            # significant
s2 <- s1[is.finite(side_log2FC) & side_log2FC != 0]                    # valid effect
s3 <- s2[mean_count_N >= MINC & mean_count_T >= MINC]                  # expressed (= hit buckets)
g2 <- s3[, .(nb = uniqueN(bucket)), by = .(cancer, gene_id)][nb >= 2]  # >=2 buckets
hits <- s3[g2[, .(cancer, gene_id)], on = c('cancer', 'gene_id')]      # buckets of switch-eligible genes

cnt <- function(dt) dt[, .(features = .N, genes = uniqueN(gene_id)), by = cancer]
F <- Reduce(function(a, b) merge(a, b, by = 'cancer', all = TRUE),
  list(cnt(t1)[, .(cancer, feat_tested = features, gene_tested = genes)],
       cnt(s1)[, .(cancer, feat_sig = features, gene_sig = genes)],
       cnt(s2)[, .(cancer, feat_fc = features, gene_fc = genes)],
       cnt(s3)[, .(cancer, feat_expr = features, gene_expr = genes)],
       g2[, .(gene_ge2 = .N), by = cancer]))
for (j in names(F)) F[is.na(get(j)), (j) := 0L]

# ---- direction + switch_type per gene ---------------------------------------
trans <- hits[, {
  up <- side_log2FC > 0
  dst <- bucket[which.max(side_log2FC)]; src <- bucket[which.min(side_log2FC)]; nb <- uniqueN(bucket)
  list(switch_type = if (any(up) && any(!up)) 'cross-compartment'
                     else if (all(up)) 'upregulated and switch' else 'downregulated and switch',
       transition = if (nb >= 2 && src != dst) paste0(src, '->', dst) else NA_character_,
       src = src, dst = dst, n_sig_buckets = nb)
}, by = .(cancer, gene_id)]
trans <- trans[!is.na(transition)]
F <- merge(F, trans[, .(gene_switch = .N), by = cancer], by = 'cancer', all.x = TRUE)
F[is.na(gene_switch), gene_switch := 0L]

reading <- fread(sep = '|', header = TRUE, text = '
transition|reading
Int->TM|Neo-surface target: product surfaces; new antibody/CAR target appears
Sec->TM|Target gain (anchored): soluble form becomes membrane-bound
TM->Int|Target loss / escape: antigen internalized, leaves the surface
TM->Sec|Decoy formation / shedding: surface antigen becomes a soluble decoy
Int->Sec|Neo-secreted form: possible decoy/biomarker
Sec->Int|Loss of secreted form')
trans <- merge(trans, reading, by = 'transition', all.x = TRUE)

# ---- attach src/dst per-bucket stats (one row per gene) ---------------------
key <- hits[, .(cancer, gene_id, bucket, side_log2FC, dIF, regular_FDR, mean_count_N, mean_count_T)]
trans <- merge(trans, key[, .(cancer, gene_id, dst = bucket,
              dst_side_log2FC = round(side_log2FC, 3), dst_dIF = round(dIF, 3),
              dst_regFDR = signif(regular_FDR, 2), dst_meanT = round(mean_count_T, 1))],
              by = c('cancer', 'gene_id', 'dst'), all.x = TRUE)
trans <- merge(trans, key[, .(cancer, gene_id, src = bucket,
              src_side_log2FC = round(side_log2FC, 3), src_dIF = round(dIF, 3),
              src_regFDR = signif(regular_FDR, 2), src_meanN = round(mean_count_N, 1))],
              by = c('cancer', 'gene_id', 'src'), all.x = TRUE)

# ---- gene symbol ------------------------------------------------------------
strip_v <- function(x) sub('\\.[0-9]+$', '', x)
if (!file.exists(TX2G)) stop(sprintf('tx2gene not found: %s (override with --tx2gene)', TX2G))
g2n <- unique(fread(TX2G, select = c('gene_id', 'gene_name')))
g2n[, gene_id := strip_v(gene_id)]
name_of <- setNames(g2n$gene_name, g2n$gene_id)[!duplicated(g2n$gene_id)]
trans[, gene_name := name_of[strip_v(gene_id)]]
trans[is.na(gene_name) | gene_name == '', gene_name := gene_id]

setcolorder(trans, c('cancer', 'gene_id', 'gene_name', 'transition', 'switch_type', 'reading',
  'src', 'dst', 'dst_side_log2FC', 'dst_dIF', 'dst_regFDR', 'dst_meanT',
  'src_side_log2FC', 'src_dIF', 'src_regFDR', 'src_meanN', 'n_sig_buckets'))
setorder(trans, cancer, switch_type, -dst_side_log2FC)

# ---- write per-cancer + combined --------------------------------------------
for (s in slugs) {
  d <- trans[cancer == s]
  fwrite(d, file.path(ROOT, s, 'tier1', sprintf('%s_tier1_filtered_hits.csv', s)))
}
fwrite(trans, file.path(COH, 'tier1_filtered_hits_ALL.csv'))

# ---- markdown report --------------------------------------------------------
md_tab <- function(df) { df <- as.data.frame(df); df[] <- lapply(df, as.character); cn <- names(df)
  paste(c(paste0('| ', paste(cn, collapse = ' | '), ' |'),
          paste0('|', paste(rep('---', length(cn)), collapse = '|'), '|'),
          apply(df, 1, function(r) paste0('| ', paste(r, collapse = ' | '), ' |'))), collapse = '\n') }
pc <- function(a, b) ifelse(b > 0, sprintf('%.1f%%', 100 * a / b), '-')

Ftot <- F[, lapply(.SD, sum), .SDcols = -'cancer'][, cancer := 'TOTAL']
Fall <- rbind(F[order(-gene_switch)], Ftot, use.names = TRUE)
gene_funnel <- Fall[, .(cancer, tested = gene_tested,
  `FDR<0.05` = gene_sig, `valid FC` = gene_fc, `expr>=10` = gene_expr,
  `>=2 buckets` = gene_ge2, `switch (final)` = gene_switch,
  `% kept` = pc(gene_switch, gene_tested))]
feat_funnel <- Fall[, .(cancer, tested = feat_tested, `FDR<0.05` = feat_sig,
  `valid FC` = feat_fc, `expr>=10 (hit buckets)` = feat_expr)]
sw  <- trans[, .N, by = switch_type][order(-N)]
tr  <- trans[, .N, by = .(transition, reading)][order(-N)]

rep <- c(
  '# Localization-switch filtering report', '',
  sprintf('_Generated %s · %d cancers · FDR < %.2g · min mean count = %g (both conditions)_',
          Sys.Date(), length(slugs), FDR, MINC), '',
  '## Definition', '',
  'A gene is called a **localization switch** in a cancer when it reallocates isoform usage',
  'across subcellular compartments (**TM** = transmembrane/surface, **Sec** = secreted,',
  '**Int** = internal) between matched normal and tumor samples.', '',
  '## Filtering steps', '',
  '1. **Significant**, `regular_FDR < 0.05`: satuRn’s test for a tumor-versus-normal change in',
  '   usage of that bucket, where `regular_FDR` is the Benjamini-Hochberg adjustment of the raw',
  '   satuRn p-value.',
  '2. **Valid effect**, `side_log2FC` finite and non-zero. Direction is read from `side_log2FC`,',
  '   so a zero or undefined value carries none.',
  sprintf('3. **Expressed**, mean count ≥ %g in both normal and tumor. Applied per bucket, which is', MINC),
  '   where a near-empty compartment inflates the fold change through the +1 pseudocount.',
  '4. **≥ 2 surviving buckets per gene**, since a switch is a reallocation between compartments;',
  '   a single surviving bucket is up- or down-regulation of one compartment.',
  '5. **Directional**, `src ≠ dst`. Destination = the bucket with the largest `side_log2FC` (gained',
  '   the most share), source = the smallest (lost the most), since share_i = FC_i / FC_gene.', '',
  '**switch_type** describes the expression background and is not a filter: *cross-compartment*',
  '(buckets move in opposite directions at roughly constant gene total), *upregulated and switch*',
  '(all buckets up), or *downregulated and switch* (all buckets down).', '',
  '## Gene-level funnel (how many genes survive each step)', '',
  md_tab(gene_funnel), '',
  '## Feature-level funnel (gene × bucket rows)', '',
  md_tab(feat_funnel), '',
  '## Final switch calls', '',
  sprintf('**%d switch genes** across %d cancers (combined table: `_cohort/tier1_filtered_hits_ALL.csv`).',
          nrow(trans), length(slugs)), '',
  '### By switch type', '', md_tab(sw), '',
  '### By transition (src→dst)', '', md_tab(tr), '',
  '## Outputs', '',
  '- `<cancer>/tier1/<cancer>_tier1_filtered_hits.csv`: per-cancer gene-level switch calls',
  '- `_cohort/tier1_filtered_hits_ALL.csv`: all cancers combined',
  '- `_cohort/tier1_filtering_report.md`: this report')
writeLines(rep, file.path(COH, 'tier1_filtering_report.md'))

logf('\nFINAL: %d switch genes across %d cancers', nrow(trans), length(slugs))
logf('Per-cancer files + combined + report written under the output root')
print(gene_funnel)
