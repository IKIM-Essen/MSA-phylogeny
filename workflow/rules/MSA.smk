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
        sections="resources/MN908947.sections.fasta",
    log:
        "logs/extract_regions/main.log",
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk subseq {input.sequence} {input.bed} > {output.sections} 2> {log}"
        
    
rule rename_regions:
    input:
        sections="resources/MN908947.sections.fasta",
        bed="resources/annotation.bed",
    output:
        renamed_sections="resources/MN908947.renamed-sections.fasta",
    log:
        "logs/rename-regions/main.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename-regions.py"


rule align_genomes:
    input:
        query=lambda wildcards: get_genomes(wildcards),
        target="resources/MN908947.renamed-sections.fasta",
    output:
        "results/aligned~main/{sample}.sam",
    log:
        "logs/aligned~main/{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -a {input.target} {input.query} -o {output} 2> {log}"


rule freebayes:
    input:
        ref="resources/genomes/MN908947.fasta",
        ref_idx="resources/genomes/MN908947.fasta.fai",
        samples="results/aligned~main/{sample}.sorted.bam",
        index="results/aligned~main/{sample}.sorted.bam.bai",
    output:
        "results/var-calls/ref~main/{sample}.bcf",
    params:
        extra=("--min-alternate-count 1"),
    log:
        "logs/freebayes/ref~main/{sample}.log",
    wrapper:
        "0.68.0/bio/freebayes"


# rule bed_intersect:
#     input:
#         gff="resources/annotation.gff.gz",
#         bcf="results/var-calls/ref~main/{sample}.bcf",
#     output:
#         region="results/region-of-interest/{sample}~{region}.vcf",
#     log:
#         "logs/bedtools-intersect/{sample}~{region}.log",
#     conda:
#         "../envs/bedtools.yaml"
#     shell:
#         "bedtools intersect -a <(bcftools view {input.bcf}) -b <(zcat {input.gff} | "
#         "grep 'ame={wildcards.region}' |  "
#         "awk -F '\\t' '{{print $4}}' | "
#         "uniq | "
#         "grep -f - <(zcat {input.gff}) | "
#         "bedtools merge -i /dev/stdin) -bed > {output} 2> {log}"


rule bed_intersect:
    input:
        gff="resources/annotation.gff.gz",
        bam="results/aligned~main/{sample}.sorted.bam",
    output:
        region="results/region-of-interest/{sample}~{region}.bed",
    log:
        "logs/bedtools-intersect/{sample}~{region}.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -abam {input.bam} -b <(zcat {input.gff} | "
        "grep 'ame={wildcards.region}' |  "
        "awk -F '\\t' '{{print $4}}' | "
        "uniq | "
        "grep -f - <(zcat {input.gff}) | "
        "bedtools merge -i /dev/stdin) -u > {output} 2> {log}"


# rule bed_intersect:
#     input:
#         gff="resources/annotation~S.gff",
#         bcf="results/var-calls/ref~main/{sample}.bcf",
#     output:
#         region="results/region-of-interest/{sample}.vcf",
#     log:
#         "logs/bedtools-intersect/{sample}.log",
#     conda:
#         "../envs/bedtools.yaml"
#     shell:
#         "bedtools intersect -a <(bcftools view {input.bcf}) -b {input.gff}  > {output} 2> {log} "
#         # "bedtools merge -i /dev/stdin > {output} 2> {log}"


