use rule align_genomes as align_genomes2 with:
    input:
        query=lambda wildcards: get_genomes(wildcards.sample),
        target="resources/genomes/MN908947.fasta",
    output:
        "results/{tag}/aligned~main~genome/{sample}.sam",
    log:
        "logs/{tag}/whole-genome/{sample}.log"


rule bcftools_mpileup:
    input:
        alignments="results/{tag}/aligned~main~genome/{sample}.sorted.bam",
        ref="resources/genomes/MN908947.fasta",
        index="resources/genomes/MN908947.fasta.fai",
    output:
        pileup="results/{tag}/varcalls~main/{sample}.pileup.bcf",
    params:
        uncompressed_bcf=False,
        extra="--max-depth 100 --min-BQ 15",
    log:
        "logs/{tag}/variant-calling/{sample}.pileup.log"
    wrapper:
        "v1.12.2/bio/bcftools/mpileup"
    

rule bcftools_call:
    input:
        pileup="results/{tag}/varcalls~main/{sample}.pileup.bcf",
    output:
        calls="results/{tag}/varcalls~main/{sample}.calls.bcf",
    params:
        uncompressed_bcf=False,
        caller="-m",  # valid options include -c/--consensus-caller or -m/--multiallelic-caller
        extra="--ploidy 1",
    log:
        "logs/{tag}/variant-calling/{sample}.calls.log"
    wrapper:
        "v1.12.2/bio/bcftools/call"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    params:
        release=100,
    log:
        "logs/vep-plugins.log",
    wrapper:
        "v1.12.2/bio/vep/plugins"


rule annotate_variants:
    input:
        calls="results/{tag}/varcalls~main/{sample}.calls.bcf",
        plugins="resources/vep/plugins",
        fasta="resources/genomes/MN908947.fasta",
        fai="resources/genomes/MN908947.fasta.fai",
        gff="resources/annotation.gff.gz",
        csi="resources/annotation.gff.gz.tbi",
    output:
        calls="results/{tag}/varcalls~main/{sample}.calls.annotated.bcf",  # .vcf, .vcf.gz or .bcf
        stats="results/{tag}/varcalls~main/{sample}.calls.annotated.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["LoFtool"],
        extra="--everything",  # optional: extra arguments
    log:
        "logs/{tag}/vep/{sample}.log",
    threads: 4
    wrapper:
        "v1.12.2/bio/vep/annotate"


rule plot_variants:
    input:
        calls="results/{tag}/varcalls~main/{sample}.calls.annotated.bcf",
    output:
        plot="results/{tag}/plots/{sample}.varplot.png",
    log:
        "logs/{tag}/varplot/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-variants.py"
