conda activate project_env

conda install -c bioconda \
  sra-tools \
  fastqc \
  multiqc \
  trim-galore \
  bowtie2 \
  samtools \
  bismark


mkdir -p ~/references/genome

wget ftp://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O ~/references/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#genome preparation
bismark_genome_preparation ~/references/genome/


conda activate project_env

conda install -c bioconda \
  sra-tools \
  fastqc \
  multiqc \
  trim-galore \
  bowtie2 \
  samtools \
  bismark

mkdir -p ~/references/genome

wget ftp://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O ~/references/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#genome preparation
bismark_genome_preparation ~/references/genome/

mkdir -p ~/PROJECTS/THESIS/BS_seq/fastq                       #For section 1 and 2
mkdir -p ~/PROJECTS/THESIS/BS_seq/fastqc_results/not_trimmed  #For section 3
mkdir -p ~/PROJECTS/THESIS/BS_seq/fastq_trimmed               #For section 4
mkdir -p ~/PROJECTS/THESIS/BS_seq/fastqc_results/trimmed      #For section 4
mkdir -p ~/PROJECTS/THESIS/BS_seq/mapped                      #For section 5


mapfile -t srr_codes_bs < ~/PROJECTS/THESIS/BS_seq/SRR_Acc_List_BS_seq.txt

#Every sample is made up of 4 files, all successive in the list, for a total of 48 files for 12 samples
#This cycle is for the samples
for((s=0; s<48; s+=4)); do   

        ### 1. Downloading and compressing fastq files ###

        #This cycle is for the single files
        for((i=s; i<s+4; i++)); do
                echo "Downloading file: $i"
                #Prefetch and download
                prefetch ${srr_codes_bs[$i]} --output-directory ~/PROJECTS/THESIS/BS_seq/fastq
                fasterq-dump ${srr_codes_bs[$i]} --outdir ~/PROJECTS/THESIS/BS_seq/fastq --progress
                #Compressing
                gzip "$HOME/PROJECTS/THESIS/BS_seq/fastq/${srr_codes_bs[$i]}.fastq"
                
                #removing the prefetch dir
                rm -r ~/PROJECTS/THESIS/BS_seq/fastq/${srr_codes_bs[$i]}
        done


        ### 2. Merging ###
        j=$((s/4 + 1))
        echo "Merging sample $j"

        zcat \
        "$HOME/PROJECTS/THESIS/BS_seq/fastq/${srr_codes_bs[s]}.fastq.gz" \
        "$HOME/PROJECTS/THESIS/BS_seq/fastq/${srr_codes_bs[s+1]}.fastq.gz" \
        "$HOME/PROJECTS/THESIS/BS_seq/fastq/${srr_codes_bs[s+2]}.fastq.gz" \
        "$HOME/PROJECTS/THESIS/BS_seq/fastq/${srr_codes_bs[s+3]}.fastq.gz" \
        | gzip > "$HOME/PROJECTS/THESIS/BS_seq/fastq/sample_${j}.fastq.gz"

        rm "$HOME/PROJECTS/THESIS/BS_seq/fastq/${srr_codes_bs[s]}.fastq.gz" "$HOME/PROJECTS/THESIS/BS_seq/fastq/${srr_codes_bs[s+1]}.fastq.gz" \
        "$HOME/PROJECTS/THESIS/BS_seq/fastq/${srr_codes_bs[s+2]}.fastq.gz" "$HOME/PROJECTS/THESIS/BS_seq/fastq/${srr_codes_bs[s+3]}.fastq.gz"


        ### 3. Quality Control (QC) ###

        fastqc -t 2 -o ~/PROJECTS/THESIS/BS_seq/fastqc_results/not_trimmed ~/PROJECTS/THESIS/BS_seq/fastq/sample_${j}.fastq.gz


        ### 4. Trimming ###
        trim_galore --output_dir ~/PROJECTS/THESIS/BS_seq/fastq_trimmed ~/PROJECTS/THESIS/BS_seq/fastq/sample_${j}.fastq.gz
        #Removing the original file after trimming
        rm ~/PROJECTS/THESIS/BS_seq/fastq/sample_${j}.fastq.gz

        #QC post-trimming
        fastqc -t 2 -o ~/PROJECTS/THESIS/BS_seq/fastqc_results/trimmed ~/PROJECTS/THESIS/BS_seq/fastq_trimmed/sample_${j}_trimmed.fastq.gz


        ### 5. Mapping and methylation calling ###
        #mapping
        bismark --genome ~/references/genome/ \
                -o ~/PROJECTS/THESIS/BS_seq/mapped \
                --gzip \
                ~/PROJECTS/THESIS/BS_seq/fastq_trimmed/sample_${j}_trimmed.fastq.gz

        #Deduplication
        deduplicate_bismark --bam \
                --single \
                -o ~/PROJECTS/THESIS/BS_seq/mapped \
                ~/PROJECTS/THESIS/BS_seq/mapped/sample_${j}_trimmed_bismark_bt2.bam

        # Methylation calling
        bismark_methylation_extractor \
                --single-end \
                --comprehensive \
                --CX_context \
                --cytosine_report \
                --genome_folder ~/references/genome \
                --gzip \
                -o ~/PROJECTS/THESIS/BS_seq/mapped \
                ~/PROJECTS/THESIS/BS_seq/mapped/sample_${j}_trimmed_bismark_bt2.deduplicated.bam
done

# Merging qc files:
#raw
multiqc ~/PROJECTS/THESIS/BS_seq/fastqc_results/not_trimmed -n multiqc_report_not_trimmed -o ~/PROJECTS/THESIS/BS_seq/fastqc_results/not_trimmed
#trimmed
multiqc ~/PROJECTS/THESIS/BS_seq/fastqc_results/trimmed -n multiqc_report_trimmed -o ~/PROJECTS/THESIS/BS_seq/fastqc_results/trimmed

# 5. Report
cd ~/PROJECTS/THESIS/BS_seq/mapped/
bismark2report
bismark2summary