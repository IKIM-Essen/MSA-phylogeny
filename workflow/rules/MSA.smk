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
        query=lambda wildcards: get_genomes(wildcards),
        target="resources/MN908947.renamed-regions.fasta",
    output:
        "results/aligned~main/{sample}.sam",
    log:
        "logs/aligned~main/{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -a {input.target} {input.query} -o {output} 2> {log}"


checkpoint extract_sequence:
    input:
        bam="results/aligned~main/{sample}.sorted.bam",
        index="results/aligned~main/{sample}.sorted.bam.bai",
    output:
        seq="results/region-of-interest/{sample}~{region}.fasta",
    log:
        "logs/extract-sequence/{sample}~{region}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract-sequences.py"


rule cat_genomes_for_MSA:
    input:
        sample_genomes=lambda wildcards: expand(
            "results/region-of-interest/{sample}~{{region}}.fasta",
            sample=get_samples(),
        )
    output:
        all_genomes="results/region-of-interest/{region}.fasta",
    log:
        "logs/cat-sequences/{region}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cat {input} > {output}"


rule mafft_alignment:
    input:
        "results/region-of-interest/{region}.fasta",
    output:
        "results/region-of-interest/{region}~aligned.fasta",
    log:
        "logs/aligned~mafft/{region}.log",
    conda:
        "../envs/mafft.yaml"
    shell:
        "mafft {input} > {output}"


rule get_output:
    input:
        "results/region-of-interest/S~aligned.fasta",