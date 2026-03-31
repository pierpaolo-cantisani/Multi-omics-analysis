#!/bin/bash

### Creating conda environment ###

conda create -n project_env python=3.10 

conda activate project_env

conda install -c bioconda -c conda-forge \
  sra-tools \
  fastqc \
  multiqc \
  fastp \
  salmon

conda deactivate



### Creating directories ###

mkdir -p ~/multi_omics/rna_seq
mkdir -p ~/multi_omics/rna_seq/fastq
mkdir -p ~/multi_omics/rna_seq/fastqc_results/not_trimmed
mkdir -p ~/multi_omics/rna_seq/fastq_trimmed
mkdir -p ~/multi_omics/rna_seq/fastqc_results/trimmed
mkdir -p ~/multi_omics/rna_seq/salmon_output
mkdir -p ~/references



### Downloading and indexing reference transcriptome for salmon quantification ###

#Downloading
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz

# Indexing
salmon index -t ~/references/gencode.v49.transcripts.fa.gz -i ~/references/salmon_index_hg38_49
