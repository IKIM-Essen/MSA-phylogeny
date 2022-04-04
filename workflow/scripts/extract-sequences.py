import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

with pysam.AlignmentFile(snakemake.input.bam, "rb") as samfile, open(snakemake.output.seq, "w") as out_fasta:
    for record in samfile.fetch():
        if record.reference_name.endswith(snakemake.wildcards.region):
            out_seq = SeqRecord(
                Seq(record.query_sequence),
                id="%s_%s" %(record.query_name, record.reference_name.split("_")[-1]),
                description=""
            )
            # print(record.reference_name, record.reference_start, record.query_alignment_start, record.query_alignment_end)
            SeqIO.write(out_seq, out_fasta, "fasta")
