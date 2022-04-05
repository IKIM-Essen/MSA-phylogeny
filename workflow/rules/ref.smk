rule get_genome:
    output:
        "resources/genomes/MN908947.fasta",
    params:
        accession=config["virus-reference-genome"],
    log:
        "logs/genomes/get-genome/main.log",
    conda:
        "../envs/entrez.yaml"
    resources:
        ncbi_api_requests=1,
    shell:
        "((esearch -db nucleotide -query '{params.accession}' | "
        "efetch -format fasta > {output}) && [ -s {output} ]) 2> {log}"


rule get_genome_annotation:
    output:
        "resources/annotation.gff.gz",
    log:
        "logs/get-annotation.log",
    conda:
        "../envs/tabix.yaml"
    shell:
        # download, sort and bgzip gff (see https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html)
        "(curl -sSL ftp://ftp.ensemblgenomes.org/pub/viruses/gff3/sars_cov_2/Sars_cov_2.ASM985889v3.101.gff3.gz | "
        "zcat | grep -v '#' | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > {output}) 2> {log}"


