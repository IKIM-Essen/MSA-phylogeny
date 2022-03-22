from snakemake.utils import validate
import pandas as pd

configfile: "config/config.yaml"

def get_samples(wildcards):
    return list(pep.sample_table["sample_name"].values)

def get_genomes(wildcards):
    return pep.sample_table.loc[wildcards.sample]["file"]
