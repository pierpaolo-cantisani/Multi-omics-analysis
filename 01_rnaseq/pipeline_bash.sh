#!/bin/bash

### 1. Downloading and merging files ###
# SRR codes are in "SRR_Acc_List.txt", obtained from the SRA website.

# Downloading
while read id; do
  echo "Downloading and compressing sample: $id"

  #Downloading
  prefetch "$id"
  fasterq-dump "$id" --outdir ~/multi_omics/rna_seq/fastq
  #Compressing
  gzip "$HOME/multi_omics/rna_seq/fastq/${id}.fastq"

  # Removing the prefetch directory
  rm -r "$id"
done < ~/multi_omics/rna_seq/SRR_Acc_List.txt


# Merging
#The experiment is composed of 48 runs, 4 for each of 12 samples. Merging each quadruple of files:
mapfile -t srr_codes < ~/multi_omics/rna_seq/SRR_Acc_List.txt

for ((i=0; i<${#srr_codes[@]}; i+=4)); do
   j=$((i/4 + 1))
   echo "Merging sample $j"
   cat \
      "$HOME/multi_omics/rna_seq/fastq/${srr_codes[i]}.fastq.gz" \
      "$HOME/multi_omics/rna_seq/fastq/${srr_codes[i+1]}.fastq.gz" \
      "$HOME/multi_omics/rna_seq/fastq/${srr_codes[i+2]}.fastq.gz" \
      "$HOME/multi_omics/rna_seq/fastq/${srr_codes[i+3]}.fastq.gz" \
      > "$HOME/multi_omics/rna_seq/fastq/sample_${j}_merged.fastq.gz"
  
   #Removing non-merged files
   rm "$HOME/multi_omics/rna_seq/fastq/${srr_codes[i]}.fastq.gz" "$HOME/multi_omics/rna_seq/fastq/${srr_codes[i+1]}.fastq.gz" "$HOME/multi_omics/rna_seq/fastq/${srr_codes[i+2]}.fastq.gz" "$HOME/multi_omics/rna_seq/fastq/${srr_codes[i+3]}.fastq.gz"
done


### 2. Quality Control

#Performing QC:
fastqc -t 2 -o ~/multi_omics/rna_seq/fastqc_results/not_trimmed ~/multi_omics/rna_seq/fastq/*.fastq.gz
#Merging QC results:
multiqc ~/multi_omics/rna_seq/fastqc_results/not_trimmed -n full_qc_report -o ~/multi_omics/rna_seq/fastqc_results/not_trimmed


### 3. Trimming and second Quality Control
for ((i=1; i<13; i+=1)); do
  echo "Trimming sample: $i"

  fastp \
    -t 2 \
    --in1 "$HOME/multi_omics/rna_seq/fastq/sample_${i}_merged.fastq.gz" \
    --out1 "$HOME/multi_omics/rna_seq/fastq_trimmed/sample_${i}_trimmed.fastq.gz" \
    --qualified_quality_phred 20 \
    --length_required 40

  rm "$HOME/multi_omics/rna_seq/fastq/sample_${i}_merged.fastq.gz"

done

#Performing QC:
fastqc -t 2 -o ~/multi_omics/rna_seq/fastqc_results/trimmed ~/multi_omics/rna_seq/fastq_trimmed/*.fastq.gz
#Merging QC results:
multiqc ~/multi_omics/rna_seq/fastqc_results/trimmed -n full_qc_report_trimmed -o ~/multi_omics/rna_seq/fastqc_results/trimmed



### 4. Mapping and quantification ###
##Using Salmon

mkdir -p ~/multi_omics/rna_seq/salmon_output

for ((i=1; i<13; i+=1)); do
  echo "Quantifying sample: ${i}"

  salmon quant -l A \
        -i ~/references/salmon_index_hg38_49 \
        -r "$HOME/multi_omics/rna_seq/fastq_trimmed/sample_${i}_trimmed.fastq.gz" \
        -p 2 \
        --validateMappings \
        -o "$HOME/multi_omics/rna_seq/salmon_output/sample_${i}_quant"

done

#Check: Statistics and file format
less ~/multi_omics/rna_seq/salmon_output/SRR1725922_quant/lib_format_counts.json
head -n 20 ~/multi_omics/rna_seq/salmon_output/SRR1725922_quant/quant.sf | column -t | less

#Renaming directories for differential analysis:
n=1
for ((i=1; i<13; i+=2)); do

  mv "$HOME/multi_omics/rna_seq/salmon_output/sample_${i}_quant" "$HOME/multi_omics/rna_seq/salmon_output/rep_${n}_TB_quant"
  mv "$HOME/multi_omics/rna_seq/salmon_output/sample_$((i+1))_quant" "$HOME/multi_omics/rna_seq/salmon_output/rep_${n}_NI_quant"
  ((n++))

done
