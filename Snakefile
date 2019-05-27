import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample","unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
#validate(units, schema="schemas/units.schema.yaml")


##### target rules #####

rule all:
    input:
        expand(["salmon/{sample}-{unit}/quant_human.sf",
        "salmon/{sample}-{unit}/quant_mouse.sf"], sample=units["sample"], unit=units["unit"]),
        "deseq2/init_mouse.rds"
#        expand(["diffexp/{contrast}.diffexp.tsv",
#                "diffexp/{contrast}.ma-plot.svg"],
#               contrast=config["diffexp"]["contrasts"]),
#        "plots/pca/pca.svg"


##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
#singularity: "docker://continuumio/miniconda3"


##### setup report #####

report: "report/workflow.rst"


##### load rules #####

include: "rules/common.smk"
#include: "rules/trim.smk"
include: "rules/salmoncount.smk"
include: "rules/diffexp.smk"
