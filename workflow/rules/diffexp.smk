def get_salmon_quant(wildcards):
    res = []
    for unit in units.itertuples():
        if is_paired_end(unit.sample_name):
            res.append("results/salmon/{sample}/quant.sf".format(sample=unit.sample_name)) 
        else:
            res.append("results/salmon_SE/{sample}/quant.sf".format(sample=unit.sample_name))
    return res

rule count_matrix:
    input:
        get_salmon_quant
    output:
        "results/counts/all.tsv"
    params:
        units=units,
        transcript_map=config["ref"]["transcript_map"]
    script:
        "../scripts/count-matrix.py"

rule get_TPM:
    input:
        inpt=get_salmon_quant
    output:
        "results/TPM/all.tsv"
    params:
        units=units,
        transcript_map=config["ref"]["transcript_map"]
    script:
        "../scripts/tpm-matrix.py"

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

rule deseq2_init:
    input:
        inpt=get_salmon_quant
    output:
        "results/deseq2/init.rds",
        "results/deseq2/normcounts.csv",
        "results/deseq2/rawcounts.csv"
    params:
        species=config["ref"]["biomart_species"],
        condition=config["diffexp"]["model"],
        samples=config["samples"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: 1
    script:
        "../scripts/deseq2_init.R"

rule pca:
    input:
        "results/deseq2/init.rds"
    output:
        report("results/plots/pca/pca.svg", "../report/pca.rst")
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.log"
    script:
        "../scripts/plot-pca.R"

rule deseq2:
    input:
        "results/deseq2/init.rds"
    output:
        table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
        table_LRT=report("results/diffexp/{contrast}.LRT.diffexp.tsv", "../report/diffexp_LRT.rst"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        contrast=lambda wildcards: config["diffexp"]["contrasts"][wildcards.contrast],
        reduced_model=config["diffexp"]["reduced_model"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: 1
    script:
        "../scripts/deseq2.R"
