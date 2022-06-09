rule make_phylogeny:
    input:
        fasta=get_msa_input,
    output:
        inter_fasta=temp("results/{tag}/phylogeny~{state}/{region}/{region}~phylo"),
        treefile="results/{tag}/phylogeny~{state}/{region}/{region}~phylo.iqtree",
        nwk_tree="results/{tag}/phylogeny~{state}/{region}/{region}~phylo.treefile",
        ml_dist="results/{tag}/phylogeny~{state}/{region}/{region}~phylo.mldist",
        iqtree_log="results/{tag}/phylogeny~{state}/{region}/{region}~phylo.log",
    log:
        "logs/{tag}/iqtree/{region}~{state}.log",
    conda:
        "../envs/iqtree.yaml"
    shell:
        "cp {input.fasta} {output.inter_fasta} && "
        "iqtree -s {output.inter_fasta} > {log} 2>&1"
