from Bio import SeqIO
import pandas as pd

bed = pd.read_csv(snakemake.input.bed, delimiter="\t", header=None)

bed[1] += 1
bed["name"] = bed[0] + ":" + bed[1].astype(str) + "-" + bed[2].astype(str)
bed["feature"] = bed[9].str.extract(r';Name=([A-Z]{1,3}\d*[a-z]{0,2});')[0]
bed.set_index("name", inplace=True)

with open(snakemake.input.regions) as input_handle, open(snakemake.output.renamed_regions, "w") as output_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        record.id = record.id + "_" + bed.loc[str(record.id)]["feature"]
        record.name = ""
        record.description = ""
        SeqIO.write(record, output_handle, "fasta")