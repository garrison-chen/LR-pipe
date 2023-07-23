configfile: "config.yaml"

rule all:
   input:
        config["path"]["data"]

rule qc_reporting_minionqc:
    input:
        config["path"]["data"]
    output:
        directory("{input}/Results/minion_qc_results")
    shell:
        "Rscript packages/MinIONQC.R -i {input}/sequencing_summary*.txt -o {input}/Results/minion_qc_results"
        
rule qc_reporting_fastqc:
    input:
        #os.path.join(config["path"]["data"], 'fastq_pass/barcode01/{sample}.fastq.gz'),
        config["path"]["data"]
    output:
        directory("{input}/Results/fastqc_results")
    conda: 'fastqc'
    shell:
        "mkdir {output} && fastqc -t 8 -o {input}/Results/fastqc_results {input}/fastq_pass/barcode01/*.fastq.gz"

'''
rule qc_reporting:
    input:
        "{input}/Results/minion_qc_results",
        "{input}/Results/fastqc_results"
    output:
        "{input}/Results/qc_results",
        directory("{input}/Results/qc_results")
    run:
        shell("mkdir {input}/Results/qc_results")
        shell("mv {input}/Results/minion_qc_results {output}")
        shell("mv {input}/Results/fastqc_results {output}")
'''

rule chopper:
    input:
        config["path"]["data"]
    output:
        directory("{input}/Results/results_chopper")
    run:
        shell("for i in 01 02 03 04 05 06 07 08 09 10 11 12; do gunzip -c {input}/fastq_pass/barcode$i/*.fastq.gz | ~/chopper --minlength 2000 --maxlength 50000 | gzip > {input}/fastq_pass/filtered_reads_barcode$i.fastq.gz; done")
        shell("mkdir {output}")
        shell("mv {input}/fastq_pass/filtered_reads_barcode*.fastq.gz {output}")

'''
rule taxonomic_profiling_sourmash
'''

rule taxonomic_profiling_mmseqs2:
    input:
        config["path"]["data"]
    output:
        directory("{input}/Results/results_mmseqs")
    conda: 'mmseqs'
    shell:
        "mkdir {input}/Results/results_mmseqs && for i in 01 02 03 04 05 06 07 08 09 10 11 12; do mmseqs easy-taxonomy {input}/Results/results_chopper/filtered_reads_barcode$i.fastq.gz /scratch/SCRATCH_SAS/chen/Meta-data/Taxonomy-databases/mmseqs-databases/NR_db {input}/Results/results_mmseqs/result_barcode$i tmp --threads 48; done"
