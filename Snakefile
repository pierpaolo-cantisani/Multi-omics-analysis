#Snakefile

#### !!Important:!!
# Before running the Snakefile: run "setup.sh"


### Global variables ###

# Reading SRR codes 
SRR_codes = open("01_rnaseq/SRR_Acc_List.txt").read().splitlines()
# Variable for merging:
runs_per_sample = [SRR_codes[i:i+4] for i in range(0, len(SRR_codes), 4)]
# Variables for file/dir names
SAMPLES = range(1, 13)
REPS = range(1,7)
CONDS = ["TB", "NI"]


### Flow start ###

rule all:
    input: 
        "01_rnaseq/DE_results/DE_results.csv",
        "01_rnaseq/DE_results/RNA-Seq Graphs.pdf",
        "01_rnaseq/DE_results/session_info.txt",
        "01_rnaseq/DE_results/RNAseq_universe.csv",
        "01_rnaseq/fastqc_results/full_qc_report_not_trimmed.html",
        "01_rnaseq/fastqc_results/full_qc_report_trimmed.html",
        "02_wgbs/DM_results/DM_sites.csv",
        "02_wgbs/DM_results/WGBS Graphs.pdf",
        "02_wgbs/DM_results/session_info.txt",
        "02_wgbs/DM_results/WGBS_universe.csv",
        "03_integration_analysis/integration_results/All_intersecting_genes.csv",
        "03_integration_analysis/integration_results/Integration analysis Graphs.pdf",
        "03_integration_analysis/integration_results/Statistics_table.csv",
        "03_integration_analysis/integration_results/Statistics_expr_table.csv",
        "03_integration_analysis/integration_results/session_info_int.txt",
        "03_integration_analysis/integration_results/Common_GO_terms.csv",
        "03_integration_analysis/integration_results/Functional Enrichment Graphs.pdf",
        "03_integration_analysis/integration_results/session_info_FE.txt"


### Rules ###


##### ----    Section 1: RNA-Seq    ---- #####

### 1. Downloading and merging fastq files

rule get_files:
    output: temp("01_rnaseq/fastq/{srr}.fastq.gz")
    wildcard_constraints:
        srr = "|".join(SRR_codes)
    conda: "envs/rna_seq.yml"
    shell:
        '''
        echo "Downloading sample: {wildcards.srr}"
        prefetch {wildcards.srr} -O 01_rnaseq/fastq/
        fasterq-dump 01_rnaseq/fastq/{wildcards.srr} --outdir 01_rnaseq/fastq/
        gzip 01_rnaseq/fastq/{wildcards.srr}.fastq
        rm -r 01_rnaseq/fastq/{wildcards.srr}
        '''


rule merge_files:
    input:
        lambda wildcards: expand("01_rnaseq/fastq/{srr}.fastq.gz", srr = runs_per_sample[int(wildcards.i)-1])
    output:
        temp("01_rnaseq/fastq/sample_{i}_merged.fastq.gz")
    conda: "envs/rna_seq.yml"
    shell: 
        '''
          echo "Merging sample {wildcards.i}"
          cat {input} > {output}
        '''

### 2. Quality control

rule fastqc_raw:
    input:
        "01_rnaseq/fastq/sample_{i}_merged.fastq.gz"
    output:
        html = "01_rnaseq/fastqc_results/not_trimmed/sample_{i}_merged_fastqc.html",
        zip = "01_rnaseq/fastqc_results/not_trimmed/sample_{i}_merged_fastqc.zip"
    threads: 2
    conda: "envs/rna_seq.yml"
    shell:
        "fastqc -t {threads} -o 01_rnaseq/fastqc_results/not_trimmed {input}"



rule multiqc_raw:
    input:
        expand("01_rnaseq/fastqc_results/not_trimmed/sample_{i}_merged_fastqc.zip", i=SAMPLES)
    output:
        "01_rnaseq/fastqc_results/full_qc_report_not_trimmed.html"
    conda: "envs/rna_seq.yml"
    shell:
        "multiqc 01_rnaseq/fastqc_results/not_trimmed -n full_qc_report_not_trimmed -o 01_rnaseq/fastqc_results/"

### 3. Trimming and second Quality Control

rule trimming:
    input: 
        "01_rnaseq/fastq/sample_{i}_merged.fastq.gz"
    output:
        fastq = "01_rnaseq/fastq_trimmed/sample_{i}_trimmed.fastq.gz",
        html = "01_rnaseq/fastq_trimmed/sample_{i}_fastp.html",
        json = "01_rnaseq/fastq_trimmed/sample_{i}_fastp.json"
    threads: 2
    conda: "envs/rna_seq.yml"
    shell:
        '''
        echo "Trimming sample: {wildcards.i}"
        fastp \
            -w {threads} \
            --in1 "{input}" \
            --out1 "{output.fastq}" \
            --html "{output.html}" \
            --json "{output.json}" \
            --qualified_quality_phred 20 \
            --length_required 40
        '''


rule fastqc_trimmed:
    input:
        "01_rnaseq/fastq_trimmed/sample_{i}_trimmed.fastq.gz"
    output:
        html = "01_rnaseq/fastqc_results/trimmed/sample_{i}_trimmed_fastqc.html",
        zip = "01_rnaseq/fastqc_results/trimmed/sample_{i}_trimmed_fastqc.zip"
    threads: 2
    conda: "envs/rna_seq.yml"
    shell:
         "fastqc -t {threads} -o 01_rnaseq/fastqc_results/trimmed {input}"


rule multiqc_trimmed:
    input:
        expand("01_rnaseq/fastqc_results/trimmed/sample_{i}_trimmed_fastqc.zip", i=SAMPLES)
    output:
        "01_rnaseq/fastqc_results/full_qc_report_trimmed.html"
    conda: "envs/rna_seq.yml"
    shell:
        "multiqc 01_rnaseq/fastqc_results/trimmed -n full_qc_report_trimmed -o 01_rnaseq/fastqc_results/"

### 4. Mapping and quantification ###

rule quantification:
    input: 
        "01_rnaseq/fastq_trimmed/sample_{i}_trimmed.fastq.gz"
    output: 
        directory("01_rnaseq/salmon_output/sample_{i}_quant")
    threads: 2
    conda: "envs/rna_seq.yml"
    shell:
        '''
        echo "Quantifying sample: {wildcards.i}"

        salmon quant -l A \
            -i 01_rnaseq/references/salmon_index_hg38_49 \
            -r {input} \
            -p {threads} \
            --validateMappings \
            -o {output}
        '''

rule Rename_dir:
    input:
        expand("01_rnaseq/salmon_output/sample_{i}_quant", i=SAMPLES)
    output:
        expand("01_rnaseq/salmon_output/rep_{n}_{cond}_quant/quant.sf", n=REPS, cond=CONDS)
    run:
        import os
        n = 1
        for j in range(1, 13, 2):
            src_tb = f"01_rnaseq/salmon_output/sample_{j}_quant"
            src_ni = f"01_rnaseq/salmon_output/sample_{j+1}_quant"
            if os.path.exists(src_tb):
                os.rename(src_tb, f"01_rnaseq/salmon_output/rep_{n}_TB_quant")
            if os.path.exists(src_ni):
                os.rename(src_ni, f"01_rnaseq/salmon_output/rep_{n}_NI_quant")
            n += 1


### 5. DE analysis on R

rule rnaseq_analysis:
    input:
        expand("01_rnaseq/salmon_output/rep_{n}_{cond}_quant/quant.sf", n=REPS, cond=CONDS)
    output:
        "01_rnaseq/DE_results/DE_results.csv",
        "01_rnaseq/DE_results/RNA-Seq Graphs.pdf",
        "01_rnaseq/DE_results/session_info.txt",
        "01_rnaseq/DE_results/RNAseq_universe.csv"
    conda:
        "envs/r_DE_analysis.yml"
    script:
        "01_rnaseq/scripts/DE_analysis.R"


##### ----    Section 2: WGBS    ---- #####
#WGBS data are already preprocessed and methylation-called: proceeding directly with the DM analysis

### 6. DM analysis on R

rule wgbs_analysis:
    input:
        "02_wgbs/dataset/GSM1565939_DC81_MTB_5mC.txt.gz", "02_wgbs/dataset/GSM1565940_DC81_NI_5mC.txt.gz", "02_wgbs/dataset/GSM1565941_DC82_MTB_5mC.txt.gz",
        "02_wgbs/dataset/GSM1565942_DC82_NI_5mC.txt.gz", "02_wgbs/dataset/GSM1565943_DC83_MTB_5mC.txt.gz", "02_wgbs/dataset/GSM1565944_DC83_NI_5mC.txt.gz",
        "02_wgbs/dataset/GSM1565945_DC87_MTB_5mC.txt.gz", "02_wgbs/dataset/GSM1565946_DC87_NI_5mC.txt.gz", "02_wgbs/dataset/GSM1565947_DC89_MTB_5mC.txt.gz",
        "02_wgbs/dataset/GSM1565948_DC89_NI_5mC.txt.gz", "02_wgbs/dataset/GSM1565949_DC91_MTB_5mC.txt.gz", "02_wgbs/dataset/GSM1565950_DC91_NI_5mC.txt.gz"
    output:
        "02_wgbs/DM_results/DM_sites.csv",
        "02_wgbs/DM_results/WGBS Graphs.pdf",
        "02_wgbs/DM_results/session_info.txt",
        "02_wgbs/DM_results/WGBS_universe.csv"
    conda:
        "envs/r_DM_analysis.yml"
    script:
        "02_wgbs/scripts/DMR_analysis.R"


##### ----    Section 3: Integration    ---- #####

### 7. Integration analysis of RNA-Seq DE analysis and WGBS DMR analysis on R

rule integration_analysis:
    input:
        "01_rnaseq/DE_results/DE_results.csv",
        "02_wgbs/DM_results/DM_sites.csv",
        "01_rnaseq/DE_results/RNAseq_universe.csv",
        "02_wgbs/DM_results/WGBS_universe.csv"
    output:
        "03_integration_analysis/integration_results/All_intersecting_genes.csv",
        "03_integration_analysis/integration_results/Integration analysis Graphs.pdf",
        "03_integration_analysis/integration_results/Statistics_table.csv",
        "03_integration_analysis/integration_results/Statistics_expr_table.csv",
        "03_integration_analysis/integration_results/session_info_int.txt"
    conda:
        "envs/r_integration_analysis.yml"
    script:
        "03_integration_analysis/scripts/Integration_analysis.R"

### 8. Functional enrichment analysis (BP) and integration of resulting GO terms

rule FE_analysis:
    input:
        "01_rnaseq/DE_results/DE_results.csv",
        "02_wgbs/DM_results/DM_sites.csv",
        "01_rnaseq/DE_results/RNAseq_universe.csv",
        "02_wgbs/DM_results/WGBS_universe.csv"
    output:
        "03_integration_analysis/integration_results/Common_GO_terms.csv",
        "03_integration_analysis/integration_results/Functional Enrichment Graphs.pdf",
        "03_integration_analysis/integration_results/session_info_FE.txt"
    conda:
        "envs/r_integration_analysis.yml"
    script:
        "03_integration_analysis/scripts/FE_analysis.R"