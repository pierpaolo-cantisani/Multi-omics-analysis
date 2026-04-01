#!/bin/bash
set -e

### Creating conda environment ###

conda env create -f ~/Multi-omics-analysis/01_rnaseq/rna_seq.yml


### Creating directories ###
mkdir -p ~/Multi-omics-analysis/01_rnaseq/{fastq,fastqc_results/not_trimmed,fastq_trimmed,fastqc_results/trimmed,salmon_output}
mkdir -p ~/Multi-omics-analysis/references


### Downloading and indexing reference transcriptome for salmon quantification ###

#Downloading
wget -P ~/Multi-omics-analysis/references/ \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz

# Indexing
salmon index \
  -t ~/Multi-omics-analysis/references/gencode.v49.transcripts.fa.gz \
  -i ~/Multi-omics-analysis/references/salmon_index_hg38_49 \
  --threads 4