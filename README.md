# Multi-Omics Integration Analysis Pipeline

A reproducible Snakemake pipeline for the integrative analysis of **RNA-Seq** (differential expression) and **WGBS** (differential methylation) data, followed by functional enrichment analysis of shared biological signals. The study compares *Mycobacterium tuberculosis*-infected (TB) vs. non-infected (NI) human dendritic cells.

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Requirements](#requirements)
- [Setup](#setup)
- [Running the Pipeline](#running-the-pipeline)
- [Pipeline Modules](#pipeline-modules)
  - [1 — RNA-Seq Analysis](#1--rna-seq-analysis)
  - [2 — WGBS Analysis](#2--wgbs-analysis)
  - [3 — Integration Analysis](#3--integration-analysis)
- [Outputs](#outputs)
- [Data Sources](#data-sources)

---

## Overview

```
Raw RNA-Seq reads (SRA)
        │
        ▼
  Download & Merge  ──►  QC (FastQC / MultiQC)
        │
        ▼
  Trimming (fastp)  ──►  QC (FastQC / MultiQC)
        │
        ▼
  Quantification (Salmon)
        │
        ▼
  DE Analysis (DESeq2)           Pre-processed WGBS data
        │                                  │
        │                                  ▼
        │                      DM Analysis (methylKit)
        │                                  │
        └──────────────┬───────────────────┘
                       ▼
             Integration Analysis (R)
                       │
                       ▼
         Functional Enrichment (clusterProfiler)
```

The pipeline is fully automated via Snakemake, with each section managed by dedicated Conda environments to ensure reproducibility.

---

## Project Structure

```
Multi-omics-analysis/
│
├── Snakefile
├── setup.sh
│
├── envs/
│   ├── rna_seq.yml
│   ├── r_DE_analysis.yml
│   ├── r_DM_analysis.yml
│   └── r_integration_analysis.yml
│
├── 01_rnaseq/
│   ├── SRR_Acc_List.txt             # List of SRA accession codes (48 runs, 12 samples)
│   ├── scripts/
│   │   └── DE_analysis.R
│   ├── references/                  # Populated by setup.sh
│   ├── fastq/                       # Populated at runtime (temp)
│   ├── fastq_trimmed/               # Populated at runtime (temp)
│   ├── fastqc_results/              # QC reports
│   ├── salmon_output/               # Salmon quantification directories
│   └── DE_results/                  # Final DE outputs
│
├── 02_wgbs/
│   ├── dataset/                     # Pre-processed bisulfite sequencing files (GEO)
│   ├── scripts/
│   │   └── DMR_analysis.R
│   ├── references/                  # Populated by setup.sh (liftOver chain)
│   └── DM_results/                  # Final DM outputs
│
└── 03_integration_analysis/
    ├── scripts/
    │   ├── Integration_analysis.R
    │   └── FE_analysis.R
    └── integration_results/         # Final integration outputs
```

---

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) ≥ 7.0
- [Conda](https://docs.conda.io/) / [Mamba](https://mamba.readthedocs.io/) (recommended)
- Internet access for reference downloads and SRA data retrieval

All bioinformatics tools (FastQC, MultiQC, fastp, Salmon, SRA-toolkit, R packages) are managed automatically through the Conda environments defined in `envs/`.

---

## Setup

Before running the pipeline, execute the setup script **once** to create all required directories, download reference files, and build the Salmon index:

```bash
bash setup.sh
```

This script performs the following steps:

1. Creates the directory structure under `~/Multi-omics-analysis/`
2. Downloads the **GENCODE v49** human transcriptome FASTA (for Salmon indexing) and GTF annotation (for DESeq2)
3. Builds the **Salmon index** for hg38 (GENCODE v49) using the `rna_seq` Conda environment
4. Downloads the **hg19→hg38 liftOver chain** file required for coordinate conversion in the WGBS analysis

> !! The Salmon indexing step requires the `rna_seq` Conda environment to be available. The script creates it automatically from `envs/rna_seq.yml` before indexing.
> Conda activation may not work within the scripts: in that case it is better to activate the environment directly with "> conda activate rna_seq", and then run the download parts of "setup.sh".

---

## Running the Pipeline

After setup, launch the full pipeline from the project root:

```bash
# Dry run (recommended first)
snakemake -n --use-conda

# Full run (adjust --cores as needed)
snakemake --use-conda --cores 4
```

Snakemake will automatically resolve dependencies and execute rules in the correct order, activating the appropriate Conda environment for each step.

---

## Pipeline Modules

### 1 — RNA-Seq Analysis

**Scripts:** `01_rnaseq/scripts/DE_analysis.R`  
**Environment:** `envs/rna_seq.yml`, `envs/r_DE_analysis.yml`

| Step | Tool | Description |
|------|------|-------------|
| Download | `prefetch` + `fasterq-dump` | Downloads 48 SRR runs from NCBI SRA |
| Merge | `cat` | Merges 4 runs per sample → 12 merged FASTQ files |
| QC (raw) | FastQC + MultiQC | Quality report on raw reads |
| Trimming | fastp | Quality filter (Phred ≥ 20, min length 40 bp) |
| QC (trimmed) | FastQC + MultiQC | Quality report on trimmed reads |
| Quantification | Salmon | Quasi-mapping against GENCODE v49 hg38 index |
| DE Analysis | DESeq2 | Paired design (6 replicates, TB vs NI); filtering ≥ 10 counts in ≥ 6 samples; threshold: padj < 0.05 & \|log2FC\| > 1 |

**Exploratory visualizations produced:** PCA, heatmap, per-sample expression boxplots, volcano plot (with apeglm shrinkage).

---

### 2 — WGBS Analysis

**Script:** `02_wgbs/scripts/DMR_analysis.R`  
**Environment:** `envs/r_DM_analysis.yml`

The bisulfite sequencing data (6 paired donors, TB vs NI) are already pre-processed and methylation-called (GEO accession GSE63409). Input files contain four columns: chromosome, position, methylated counts (M), and total coverage (M+U).

| Step | Tool | Description |
|------|------|-------------|
| Object construction | methylKit | Custom `methylRaw` objects from simplified bismark-like format |
| Filtering | methylKit | Coverage filter: lo ≥ 5×, hi ≤ 99.99th percentile |
| Merging | methylKit | `unite()` across all 12 samples |
| Explorative analysis | methylKit | Hierarchical clustering, PCA |
| DM analysis | methylKit | Chi-squared test with multinomial overdispersion correction (MN); thresholds: qvalue < 0.01 & \|meth.diff\| > 25% |
| Coordinate conversion | rtracklayer | LiftOver from hg19 → hg38 |
| Annotation | ChIPseeker | Genomic annotation relative to TSS (promoter: −2000 to +200 bp) using TxDb hg38 |

**Visualizations produced:** methylation/coverage stats, hierarchical clustering dendrogram, PCA, volcano plot, annotation barplot, pie chart, distance-to-TSS plot.

---

### 3 — Integration Analysis

**Scripts:** `03_integration_analysis/scripts/Integration_analysis.R`, `FE_analysis.R`  
**Environment:** `envs/r_integration_analysis.yml`

#### 3a — Integration Analysis

Intersects DE genes and DM genes (restricted to a shared universe of genes present in both analyses). The analysis is stratified by genomic region (promoter, exon, intron, 3′ UTR; intergenic sites are excluded from the primary intersection).

Statistical questions addressed (hypergeometric tests and Spearman correlations):

| Question | Test |
|----------|------|
| Is the DE ∩ DM intersection significant? | Hypergeometric |
| Do DM genes trend towards up/downregulation? | Hypergeometric |
| Is the DE ∩ DM subset more directionally biased than all DM? | Hypergeometric |
| Are hypomethylated/hypermethylated sites associated with up/downregulation? | Hypergeometric |
| Are intersecting genes enriched at specific genomic regions? | Hypergeometric |
| Does methylation magnitude correlate with expression magnitude? | Spearman |
| Do genes with ≥2 DM sites have a higher DE probability? | Hypergeometric |
| Does the number of DM sites correlate with log2FC? | Spearman |

**Visualizations produced:** Euler diagram (DE/DM overlap), up/down regulation bar charts, hypo/hypermethylation vs. expression direction chart, meth.diff vs. log2FC scatter, region distribution bar chart, methylation magnitude vs. expression scatter plots, DM site multiplicity bar chart.

#### 3b — Functional Enrichment Analysis

Over-Representation Analysis (ORA) on Biological Process (BP) GO terms, performed independently on the DE gene set and the DM gene set, followed by integration of shared significant terms.

| Step | Tool |
|------|------|
| ORA (DE genes) | clusterProfiler `enrichGO` |
| ORA (DM genes) | clusterProfiler `enrichGO` |
| Integration | Intersection of significant BP terms (padj < 0.05) |

**Visualizations produced:** barplots and dotplots of top enriched terms per omic, Euler diagram of GO term overlap, dotplot of shared top-ranked terms (scored by combined −log10(p.adjust)).

---

## Outputs

| Path | Description |
|------|-------------|
| `01_rnaseq/fastqc_results/full_qc_report_not_trimmed.html` | MultiQC report — raw reads |
| `01_rnaseq/fastqc_results/full_qc_report_trimmed.html` | MultiQC report — trimmed reads |
| `01_rnaseq/DE_results/DE_results.csv` | Significant DE genes (padj < 0.05, \|log2FC\| > 1) |
| `01_rnaseq/DE_results/RNAseq_universe.csv` | Full gene universe with log2FC and padj |
| `01_rnaseq/DE_results/RNA-Seq Graphs.pdf` | All RNA-Seq exploratory and DE visualizations |
| `01_rnaseq/DE_results/session_info.txt` | R session info |
| `02_wgbs/DM_results/DM_sites.csv` | Significant DM sites (qvalue < 0.01, \|meth.diff\| > 25%) with genomic annotation |
| `02_wgbs/DM_results/WGBS_universe.csv` | All covered gene universe |
| `02_wgbs/DM_results/WGBS Graphs.pdf` | All WGBS exploratory and DM visualizations |
| `02_wgbs/DM_results/session_info.txt` | R session info |
| `03_integration_analysis/integration_results/All_intersecting_genes.csv` | Full DE ∩ DM intersection table |
| `03_integration_analysis/integration_results/Statistics_table.csv` | Integration statistical results (by region) |
| `03_integration_analysis/integration_results/Statistics_expr_table.csv` | Expression directionality statistics |
| `03_integration_analysis/integration_results/Integration analysis Graphs.pdf` | All integration visualizations |
| `03_integration_analysis/integration_results/Common_GO_terms.csv` | Shared significant GO terms (DE + DM) |
| `03_integration_analysis/integration_results/Functional Enrichment Graphs.pdf` | All FE visualizations |
| `03_integration_analysis/integration_results/session_info_int.txt` | R session info (integration) |
| `03_integration_analysis/integration_results/session_info_FE.txt` | R session info (FE) |

---

## Data Sources

| Data | Accession | Description |
|------|-----------|-------------|
| RNA-Seq | GSE64179 | NCBI SRA (see `SRR_Acc_List.txt`) | 12 samples (6 TB, 6 NI), 4 runs/sample |
| WGBS | GSE64177 | 12 pre-processed methylation files (6 donors × 2 conditions) |
| Reference genome | GENCODE v49 (hg38) | Transcriptome FASTA + GTF annotation |
| LiftOver chain | UCSC | hg19 → hg38 coordinate conversion |
