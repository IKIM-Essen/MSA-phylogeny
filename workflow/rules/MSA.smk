rule gff2bed:
    input:
        gff="resources/annotation.gff.gz",
    output:
        bed="resources/annotation.bed",
    log:
        "logs/gff2bed/main.log",
    conda:
        "../envs/bedops.yaml"
    shell:
        "gff2bed < <(zcat {input.gff} | "
        "awk '{{if ($3 ~ /gene/) print $0}}') > {output.bed} 2> {log}"


rule extract_regions:
    input:
        sequence="resources/genomes/MN908947.fasta",
        bed="resources/annotation.bed",
    output:
        regions="resources/MN908947.regions.fasta",
    log:
        "logs/extract_regions/main.log",
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk subseq {input.sequence} {input.bed} > {output.regions} 2> {log}"


rule rename_regions:
    input:
        regions="resources/MN908947.regions.fasta",
        bed="resources/annotation.bed",
    output:
        renamed_regions="resources/MN908947.renamed-regions.fasta",
    log:
        "logs/rename-regions/main.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename-regions.py"


rule align_genomes:
    input:
        query=lambda wildcards: get_genomes(wildcards.sample),
        target="resources/MN908947.renamed-regions.fasta",
    output:
        "results/{tag}/aligned~main/{sample}.sam",
    log:
        "logs/{tag}/aligned~main/{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -a {input.target} {input.query} -o {output} 2> {log}"


rule extract_sequence:
    input:
        bam="results/{tag}/aligned~main/{sample}.sorted.bam",
        index="results/{tag}/aligned~main/{sample}.sorted.bam.bai",
    output:
        seq="results/{tag}/region-of-interest/{sample}~{region}.fasta",
    log:
        "logs/{tag}/extract-sequence/{sample}~{region}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract-sequences.py"


rule cat_genomes_for_MSA:
    input:
        lambda wildcards: get_fastas_by_region(wildcards),
    output:
        all_genomes="results/{tag}/msa-unaligned-seqs/{region}.fasta",
    log:
        "logs/{tag}/cat-sequences/{region}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cat {input} > {output}"


rule mafft_alignment:
    input:
        "results/{tag}/msa-unaligned-seqs/{region}.fasta",
    output:
        "results/{tag}/msa-aligned/{region}~aligned.fasta",
    log:
        "logs/{tag}/aligned~mafft/{region}.log",
    conda:
        "../envs/mafft.yaml"
    shell:
        "mafft {input} > {output} 2> {log}"
