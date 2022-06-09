from snakemake.utils import validate
import pandas as pd


configfile: "config/config.yaml"


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_genomes(sample):
    return pep.sample_table.loc[sample]["file"]


def get_all_tags():
    return list(pep.sample_table["tag"].unique())


def get_samples_for_tag(tag):
    df = pep.sample_table
    df = df[df["tag"] == tag]

    return list(df["sample_name"].values)


def get_fastas_by_region(wildcards):
    if wildcards.region == "genome":
        return get_genomes(sample=get_samples_for_tag(wildcards.tag))
    else:
        return expand(
            "results/{{tag}}/region-of-interest/{sample}~{{region}}.fasta",
            sample=get_samples_for_tag(wildcards.tag),
        )

def get_output_stem(wildcards, output):
    return output[0].split("_")[0]


def get_msa_input(wildcards):
    if wildcards.state == "aligned":
        return "results/{tag}/msa-aligned/{region}~aligned.fasta"
    elif wildcards.state == "cleaned":
        return "results/{tag}/msa-cleaned/{region}~aligned_cleaned.fasta",