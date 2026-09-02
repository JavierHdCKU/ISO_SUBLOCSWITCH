# ISO_SUBLOCSWITCH

**Pan-cancer subcellular-localization isoform switching.**

Tumours can change *where* a gene's protein product ends up without changing how
much of that gene is expressed. If a receptor's transcript output is redirected
from a transmembrane isoform to a secreted or intracellular one, the gene still
looks expressed in bulk RNA-seq, but the antigen leaves the cell surface. This
repository implements a test for that effect across TCGA.

Every protein-coding isoform is assigned to one of three compartments with
DeepLoc 2.1, a gene's counts are pooled within each compartment, and the
resulting compartment composition is tested for a shift between each tumour and
its matched normal with satuRn. Genes whose expression is reallocated from one
compartment to another are reported as localization switches.

| | |
|---|---|
| Compartments | **TM** (transmembrane / cell surface), **Sec** (secreted), **Int** (intracellular) |
| Expression data | TCGA, UCSC Toil recompute (Kallisto, GENCODE v23 / Ensembl 81, GRCh38) |
| Cohort | 14 cancer types, 644 matched tumour/normal patient pairs (>= 10 pairs per cancer) |
| Design | paired randomized complete block, `~ 0 + condition + patient_id` |
| Result | **3,723** gene x cancer localization switches |

---

## What is in this repository

The code is here. The data is on Zenodo (see [Getting the data](#getting-the-data)).

```
ISO_SUBLOCSWITCH/
├── README.md                     this file
├── LICENSE                       MIT (code only; the Zenodo data is licensed separately)
├── CITATION.cff                  machine-readable citation metadata
├── environment.yml               conda environment, versions pinned to the run manifests
├── .gitignore                    excludes data/ (several GB, archived on Zenodo)
│
├── scripts/
│   ├── cohort_select.R           which cancers have >= N paired patients; writes the cohort manifest
│   ├── make_tier1_input.R        per cancer: recover counts, filter, assign compartments,
│   │                             build the gene x compartment SummarizedExperiment
│   ├── run_saturn_dtu_tier1.R    per cancer: fit and test with satuRn
│   ├── filter_switch_hits.R      across all cancers: apply the switch filter, call src -> dst
│   └── lib/
│       ├── module1_prep.R        sample pairing, count recovery, expression filter,
│       │                         DeepLoc compartment assignment, SE builder
│       ├── saturn_dtu.R          the satuRn fit/test engine and the side log2FC descriptor
│       ├── qc_report.R           input QC, assertions, manifest for the SE build
│       └── qc_results.R          result QC, assertions, manifest for the DTU step
│
└── data/                         NOT in Git. See data/README.md for the Zenodo layout.
    ├── inputs/                   the 5 files the pipeline consumes (~1 GB)
    └── outputs/                  everything the pipeline produces (~2.8 GB)
```

Four entrypoints, one shared library. Each entrypoint is a standalone
command-line tool: it logs to stderr, takes `--help` and `--version`, returns a
non-zero exit code on failure, and writes a JSON manifest next to its output
recording software versions, parameters, thresholds, input provenance (filename,
size, md5) and a QC verdict, plus a table of machine-checkable assertions. A
failed hard assertion stops the pipeline rather than producing a quiet wrong
answer.

### Scope

Two upstream stages sit outside this repository:

- **Reference preparation**, which builds `deeploc.csv` and `tx2gene.tsv.gz` from
  raw Ensembl files and cleans the raw Toil matrix into the two count matrices.
  [Getting the data](#getting-the-data) records the public files it consumes and
  what it does with them; `data/inputs/` is its output.
- **Figures and downstream annotation**, which produced the contents of
  `data/outputs/figures/` and `data/outputs/hits_topology/`.

---

## Getting the data

### Option A: download the prepared inputs

All inputs and outputs are archived on Zenodo. Download them into `data/`,
preserving the layout in [`data/README.md`](data/README.md).

> **Zenodo DOI:** `10.5281/zenodo.XXXXXXX` (added on deposit)

`data/inputs/` holds five files, which is everything the pipeline reads:

| file | contents |
|---|---|
| `tumor_log2counts.tsv.gz` | TCGA Primary Tumor, log2(est_counts + 1), transcripts x samples |
| `normal_log2counts.tsv.gz` | TCGA Solid Tissue Normal, same scale and transcript set |
| `metadata.tsv.gz` | per-sample phenotype: `sample`, `detailed_category`, `_sample_type`, `_primary_site`, `_gender`, `_study` |
| `tx2gene.tsv.gz` | Ensembl 81 map: `transcript_id`, `chr`, `gene_id`, `gene_name`, `transcript_name`, `transcript_biotype`, `gene_biotype` |
| `deeploc.csv` | DeepLoc 2.1 output per isoform: `Protein_ID`, ten compartment probabilities, `Localizations`, `Signals`, `Membrane types` |

### Option B: rebuild the inputs from public sources

Everything in `data/inputs/` derives from three public resources.

**1. Expression matrices and metadata, from the UCSC Xena Toil hub.**
The Toil recompute (Vivian et al., 2017) uniformly requantified TCGA, TARGET and
GTEx with Kallisto v0.42.4 against the GENCODE v23 transcriptome. Browse the hub
at <https://xenabrowser.net/datapages/?hub=https://toil.xenahubs.net:443> and
take these datasets (each dataset page carries its own download link; the direct
pattern is `https://toil.xenahubs.net/download/<filename>`):

| file | used for |
|---|---|
| `tcga_Kallisto_est_counts.gz` | the transcript x sample count matrix, log2(est_counts + 1) |
| `TcgaTargetGTEX_phenotype.txt.gz` | sample type, cancer type, primary site |
| `probeMap/gencode.v23.annotation.transcript.probemap` | transcript to gene mapping for the matrix rows |

**2. Annotation, from the Ensembl FTP server.**
Use **release 81 (GRCh38)**, not a current release: GENCODE v23 corresponds to
Ensembl 81, so this is what makes transcript IDs join between the expression
matrix and the annotation. A newer Ensembl release will silently fail to match a
large fraction of the matrix.

```bash
BASE=https://ftp.ensembl.org/pub/release-81
wget $BASE/gtf/homo_sapiens/Homo_sapiens.GRCh38.81.gtf.gz          # annotation
wget $BASE/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz  # protein translations
```

The GTF yields `tx2gene.tsv.gz` and the protein-coding keep-list. The peptide
FASTA supplies the sequences given to DeepLoc.

**3. Localization predictions, from DeepLoc 2.1.**
Obtain DeepLoc 2.1 from
<https://services.healthtech.dtu.dk/services/DeepLoc-2.1/> and install it in its
own environment (it is not on conda and is not part of `environment.yml`).

Prepare the input by taking every transcript annotated `protein_coding` on the
primary chromosomes (1-22, X, Y, MT), pulling its translation from the peptide
FASTA, and writing one FASTA per chromosome with **plain unversioned ENST
identifiers as headers** (`>ENST00000379228`), which is what DeepLoc 2.1
expects. Ensembl peptide headers carry versioned IDs while the GTF does not, so
strip the version suffix when joining. Transcripts present in the GTF with no
peptide translation are excluded.

Then run the **Accurate** model per chromosome and concatenate the results:

```bash
deeploc2 -f <chr>.faa -o <outdir> -m Accurate -d cpu
```

The Accurate model derives its embeddings from the ProtT5-XL-UniRef50 protein
language model rather than the ESM-1b embeddings of the Fast model. The weights
are frozen and used in inference mode, so predictions are deterministic and
re-running on the same input reproduces the output exactly.

**Preparing the two count matrices.** Restrict the Toil matrix to
`protein_coding` transcripts on the primary chromosomes using the GTF keep-list,
split it into Primary Tumor and Solid Tissue Normal, and intersect the two so
both carry the same transcript set. This yields **75,144 transcripts**.

---

## Method, in brief

Full detail is in the manuscript. This is enough to read the code.

### Stage 1: cohort and testable feature set

`cohort_select.R`, then `make_tier1_input.R` per cancer.

Samples are restricted to `Primary Tumor` and `Solid Tissue Normal`. A patient is
the first 12 characters of the TCGA barcode, and is **paired** if they contribute
at least one of each. A cancer is kept if it has **>= 10 paired patients**, a
power floor fixed in advance. Where a patient
contributed several vials of one type, the choice is deterministic (sort by
barcode, take the first), so the cohort is stable across reruns.

Toil ships log2(est_counts + 1), so raw counts are recovered as `2^x - 1`. These
are kept as **fractional expected counts and are not rounded**: satuRn's
quasi-binomial quasi-likelihood model specifies only a mean and a variance and
accepts continuous counts directly.

A conservative expression filter then runs on the raw counts, mirroring the
IsoformSwitchAnalyzeR pre-filter. For each feature and condition it uses the
lower bound of the 95% confidence interval of the mean, `CI = mu - 1.96 * sd /
sqrt(n)`, which is robust to a single high-count sample dragging an otherwise
silent feature over the threshold. A gene is kept if it has >= 2 transcripts,
its gene-summed CI exceeds 10 reads in **both** conditions, and at least one of
its transcripts exceeds 10 reads in **either**.

Each transcript is assigned the compartment DeepLoc ranks highest among its ten
classes (no probability cutoff, first match on ties), then collapsed to three
bins:

- **TM** requires *both* that the top compartment is `Cell membrane` *and* that
  DeepLoc's `Membrane types` field reports a transmembrane helix. This dual
  condition matters in both directions: cell-membrane proteins that are
  peripheral, lipid-anchored or soluble carry no transmembrane domain and are
  demoted to Int, while transmembrane proteins of *intracellular* membranes (ER,
  Golgi, mitochondrion, lysosome) already fall to Int through their top
  compartment and are not pulled back in.
- **Sec** is `Extracellular`.
- **Int** is the remaining eight compartments plus the demoted cell-membrane
  transcripts.

A gene's transcripts are then summed within each compartment into a
`gene__COMPARTMENT` pseudo-feature. Because each transcript belongs to exactly
one compartment, the sums are a disjoint partition of the gene and conserve its
total, which the QC asserts. Only genes occupying **>= 2 compartments** are kept:
a gene with one destination cannot redistribute.

### Stage 2: testing

`run_saturn_dtu_tier1.R`.

satuRn fits a quasi-binomial GLM of each compartment's *share* of its gene's
total under `~ 0 + condition + patient_id`, a randomized complete block design in
which `patient_id` absorbs per-patient baseline variation, and tests the contrast
`conditionTumor - conditionNormal`. satuRn is driven through its native interface
rather than through IsoformSwitchAnalyzeR's `importRdata`, which injects SVA
surrogate variables that on a paired design correlate with `condition` and cancel
the signal.

Because satuRn's coefficient describes a change in *share* and not the direction
of the absolute count mass, a model-free descriptor is computed alongside it on
the raw compartment counts:

```
side_log2FC = log2( (mean_tumour + 1) / (mean_normal + 1) )
```

### Stage 3: switch calling

`filter_switch_hits.R`, run once across all cancers.

Per gene and compartment, in order:

1. **significant**, `regular_FDR < 0.05`;
2. **valid effect**, `side_log2FC` finite and non-zero;
3. **expressed**, mean count >= 10 in **both** conditions. This is the step that
   removes the near-empty-compartment artefact, where the +1 pseudocount inflates
   a fold change for a compartment with no real product. It has to be applied per
   compartment, because a gene-level filter cannot see it;
4. **reallocation**, the gene must retain >= 2 compartments after steps 1 to 3,
   since a single surviving compartment is up- or down-regulation, not a switch.

Direction is then read off rather than filtered. A compartment's share fold
change equals its count fold change divided by the gene-total fold change, so
ranking compartments by `side_log2FC` ranks them by share gained or lost: the
**destination** is the compartment that gained most, the **source** the one that
lost most, giving a `src -> dst` call such as TM -> Sec.

Each call also carries a `switch_type` describing its expression background,
which is descriptive and not a filter: *cross-compartment* (compartments move in
opposite directions at roughly constant gene total), *upregulated and switch*
(all move up) or *downregulated and switch* (all move down).

---

## Reproducing the analysis

```bash
conda env create -f environment.yml
conda activate sublocswitch

# download data/ from Zenodo first (see data/README.md)

Rscript scripts/cohort_select.R                                   # -> data/outputs/_cohort/
Rscript scripts/make_tier1_input.R    --cancer "Esophageal Carcinoma"
Rscript scripts/run_saturn_dtu_tier1.R --cancer-slug Esophageal_Carcinoma
# repeat the two per-cancer steps for every cancer with pass == TRUE in
# data/outputs/_cohort/cohort_manifest.tsv, then:
Rscript scripts/filter_switch_hits.R                              # -> tier1_filtered_hits_ALL.csv
Rscript scripts/cohort_select.R --summarize                       # fold per-cancer counts into the summary
```

Every script takes `--input-dir` / `--output-dir` (or `--root`) and defaults to
`data/inputs` and `data/outputs`, so the commands above work from a fresh clone
with no configuration. `--help` lists the options for any of them.

**Resources.** Each step runs as a plain command on one core; no scheduler or
cluster is needed. `make_tier1_input.R` is the heaviest, peaking at about 3.2 GB
of RAM and taking roughly 4 minutes for a 13-pair cancer, most of it spent
decompressing the tumour matrix. Larger cohorts scale with the number of samples,
not with the size of the matrix on disk.

---

## Output files

Per cancer, under `data/outputs/<Cancer>/tier1/`:

| file | contents |
|---|---|
| `*_tier1_se.rds` | the gene x compartment `SummarizedExperiment` that was tested |
| `*_tier1_qc.md` | filter funnel, compartment counts, depth balance, sparsity, DeepLoc confidence, assertions |
| `*_tier1_manifest.json` | versions, thresholds, cohort membership, input provenance, QC verdict |
| `*_tier1_dtu_results.tsv` / `.rds` | one row per gene x compartment: `estimates`, `se`, `df`, `t`, `pval`, `regular_FDR`, `q_BH`, `side_log2FC` |
| `*_tier1_dtu_qc.md`, `*_qc_pfdr.png` | p-value diagnostics and yield |
| `*_tier1_dtu_manifest.json` | engine, formula, contrast, versions, summary counts |
| `*_tier1_filtered_hits.csv` | this cancer's switch calls |

Across the cohort, under `data/outputs/_cohort/`:

| file | contents |
|---|---|
| `cohort_manifest.tsv`, `cohort_summary.md` | per-cancer paired-patient counts and which passed |
| **`tier1_filtered_hits_ALL.csv`** | **the result: 3,723 switch calls across 14 cancers** |
| `tier1_filtering_report.md` | the filter funnel per cancer, and totals by transition and switch type |
| `tier1_recurrence_ranked.csv` | one row per gene x transition, with the number of cancers it recurs in |
| `isoform_cpm_*.tsv.gz`, `isoform_tpm_*.tsv.gz` | per-isoform expression caches used by the figures |

`data/outputs/figures/` holds the manuscript figures, and
`data/outputs/hits_topology/` the per-gene panels (MAFFT alignment, DeepTMHMM
topology, DeepLoc probabilities and per-sample expression) for 72 candidate
genes, including the three shown in the paper: IL6ST, IL1R1 and TRPM2.

Columns in the switch tables: `cancer`, `gene_id`, `gene_name`, `transition`,
`switch_type`, `reading`, `src`, `dst`, then the destination and source
statistics (`*_side_log2FC`, `*_dIF`, `*_regFDR`, `dst_meanT`, `src_meanN`) and
`n_sig_buckets`. `dIF` is the change in that compartment's mean share between
conditions. `reading` is a curated biological interpretation of the transition
and is annotation, not a computed quantity.

---

## Significance

Significance is assessed on satuRn's **`regular_FDR`**, the Benjamini-Hochberg
adjustment of the raw satuRn p-value, at 0.05. `q_BH` recomputes that adjustment
with `p.adjust(pval, "BH")` as a cross-check that the two agree.

Genomic inflation of the raw p-values (`lambda_raw`), the same statistic on the
per-gene minimum p-value (`lambda_dedup`) and the estimated null fraction (`pi0`)
are recorded per cancer in the DTU QC and manifests. They are descriptive
monitoring: they do not gate, filter or reweight anything.

---

## Software

R 4.5.3 with satuRn 1.18.0, edgeR 4.8.2, SummarizedExperiment 1.40.0,
BiocParallel and data.table 1.17.8. DeepLoc 2.1 for localization, DeepTMHMM 1.0
via BioLib for topology, MAFFT v7.525 for alignment.

Versions in `environment.yml` are pinned to those recorded in the shipped
manifests, so a rebuilt environment matches the one that produced the results.
DeepLoc and DeepTMHMM are external prerequisites and are not installable through
conda; neither is needed to reproduce the analysis from `data/inputs/`.

## Citation

_Manuscript in preparation._ Citation metadata is in
[`CITATION.cff`](CITATION.cff).

## License

Code is released under the MIT License ([`LICENSE`](LICENSE)). The archived data
on Zenodo is licensed separately on the deposit; see the Zenodo record.
