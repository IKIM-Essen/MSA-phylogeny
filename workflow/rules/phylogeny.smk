rule make_phylogeny:
    input:
        fasta="results/{tag}/msa-aligned/{region}~aligned.fasta",
    output:
        inter_fasta=temp("results/{tag}/phylogeny/{region}/{region}~phylo"),
        treefile="results/{tag}/phylogeny/{region}/{region}~phylo.iqtree",
        nwk_tree="results/{tag}/phylogeny/{region}/{region}~phylo.treefile",
        ml_dist="results/{tag}/phylogeny/{region}/{region}~phylo.mldist",
        iqtree_log="results/{tag}/phylogeny/{region}/{region}~phylo.log",
    log:
        "logs/{tag}/iqtree/{region}.log",
    conda:
        "../envs/iqtree.yaml"
    shell:
        "cp {input.fasta} {output.inter_fasta} && "
        "iqtree -s {output.inter_fasta} > {log} 2>&1"