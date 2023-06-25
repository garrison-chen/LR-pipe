configfile: "config.yaml"

rule all:
    input:
        "sans_greedytree.new"

rule merge_input:
    input:
        config["path"]["data"]
    output:
        "input.fasta"
    shell:
        "cat {input}/* > {output}"

rule clustering:
    input:
        "input.fasta"
    output:
        "clusterRes_cluster.tsv"
    run:
        shell("mmseqs easy-cluster input.fasta clusterRes tmp --min-seq-id {config[params][sequence_identity]} -c 0.8 --cov-mode 1 --threads {config[params][threads]}")  
        
rule u_shape_plotting:
    input:
        clusterRes="clusterRes_cluster.tsv"
    output:
        "U-shape_plot.png"
    script:
        "Scripts/u-shape-plot.py"
        
rule select_proteins:
    input:
        "clusterRes_cluster.tsv",
        config["path"]["data"]
    output:
        "Selected_Prot.fasta"
    params:
        config["threshold"]["select_threshold"]
    script:
        "Scripts/select-proteins.py"

rule separate_proteins:
    input:
        "Selected_Prot.fasta"
    output:
        directory("toSeparate/")
    script:
        "Scripts/separate-proteins.py"

rule run_sans:
    input:
        directory("toSeparate/")
    output:
        "sans_greedytree.new"
    run:
        shell("ls -1 toSeparate/* > list.txt")
        shell("SANS -i list.txt -k 31 -f strict -N sans_greedytree.new")
        
rule run_protspam:
    input:
        directory("toSeparate/")
    output:
        "DMat"
    run:
        shell("ls -1 toSeparate/* > filelist")
        shell("~/ProtSpaM/bin/Debug/protspam -l filelist > DMat")
        
rule run_multispam:
    input:
        directory("toSeparate/")
    output:
        "multi-spam.new"
    run:
        shell("python ~/multi-SpaM/multispam.py -t <number of threads> -i <input file> -o multi-spam.new")
