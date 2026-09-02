# Data

The data for this project is **not stored on GitHub** (it is several GB, and two
files exceed GitHub's 100 MB per-file limit). It is archived on Zenodo and should
be downloaded here, preserving the layout below.

> **Zenodo DOI:** `10.5281/zenodo.XXXXXXX` (added on deposit)

```
data/
├── inputs/      # pipeline inputs (~1 GB)
│   ├── tumor_log2counts.tsv.gz     TCGA Primary-Tumor log2(est_counts + 1), Toil/Kallisto, Ensembl 81
│   ├── normal_log2counts.tsv.gz    TCGA Solid-Tissue-Normal log2(est_counts + 1)
│   ├── metadata.tsv.gz             per-sample phenotype (TcgaTargetGTEX)
│   ├── tx2gene.tsv.gz              transcript -> gene / gene_name map (Ensembl 81)
│   └── deeploc.csv                 DeepLoc 2.1 (Accurate model) predictions per isoform
└── outputs/     # regenerable from inputs/ by the pipeline in ../scripts/
    ├── _cohort/                    cohort manifest, combined switch calls (tier1_filtered_hits_ALL.csv),
    │                               recurrence ranking, per-isoform CPM/TPM caches
    │                               (isoform_cpm_persample.tsv.gz is ~2.7 GB)
    ├── <Cancer>/tier1/             per-cancer SummarizedExperiment, satuRn DTU results, switch calls, QC
    ├── figures/                    manuscript figures
    └── hits_topology/              per-gene DeepLoc / DeepTMHMM / expression panels
```

These are exactly the five files `scripts/make_tier1_input.R` expects in
`--input-dir`; the JSON manifests written next to every output record which of
them was used, by filename, size and md5.

`inputs/` derives from the public UCSC Toil TCGA recompute. Everything in
`outputs/` is regenerable from `inputs/` (see the top-level README).
