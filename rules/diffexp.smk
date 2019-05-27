rule count_matrix:
    input:
        expand("salmon/{unit.sample}-{unit.unit}/quant.sf", unit=units.itertuples())
    output:
        "counts/all.tsv"
    params:
        units=units
    script:
        "../scripts/count-matrix.py"


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        inpt=expand(
    "salmon/{unit.sample}-{unit.unit}/quant_mouse.sf", unit=units.itertuples())
    output:
        "deseq2/init.rds",
        "deseq2/normcounts.tsv",
        "deseq2/rawcounts.tsv"
    params:
        samples = config["samples"],
        species = config["species"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: 1
    script:
        "../scripts/deseq2_init.R"

rule pca:
    input:
        "deseq2/init.rds"
    output:
        report("plots/pca/pca.svg", "../report/pca.rst")
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.log"
    script:
        "../scripts/plot-pca.R"


#def get_contrast(wildcards):
#    return config["diffexp"]["contrasts"][wildcards.contrast]


rule deseq2:
    input:
        "deseq2/init.rds"
    output:
        table=report("diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
        ma_plot=report("diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        contrast=lambda wildcards: config["diffexp"]["contrasts"][wildcards.contrast]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: 1
    script:
        "../scripts/deseq2.R"
