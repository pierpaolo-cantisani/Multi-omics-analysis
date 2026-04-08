#!/bin/bash

### Creating conda environments ###

#Environment for the preprocessing on bash:
conda env create -f ~/Multi-omics-analysis/01_rnaseq/rna_seq.yml
#Environment for the DE analysis on R:
conda env create -f ~/Multi-omics-analysis/01_rnaseq/r_DE_analysis.yml


### Creating directories ###
mkdir -p ~/Multi-omics-analysis/01_rnaseq/{fastq,fastqc_results/not_trimmed,fastq_trimmed,fastqc_results/trimmed,salmon_output,references,DE_results}


### Downloading and indexing reference transcriptome for salmon quantification ###

#Downloading
wget -P ~/Multi-omics-analysis/01_rnaseq/references/ \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz

# Indexing
salmon index \
  -t ~/Multi-omics-analysis/01_rnaseq/references/gencode.v49.transcripts.fa.gz \
  -i ~/Multi-omics-analysis/01_rnaseq/references/salmon_index_hg38_49 \
  --threads 4


### Downloading the reference for the R analysis:
wget -P ~/Multi-omics-analysis/01_rnaseq/references/ \
  ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz