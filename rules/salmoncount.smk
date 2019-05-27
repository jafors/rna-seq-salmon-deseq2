#def get_fastq1(wildcards):
#    if not is_single_end(**wildcards):
#        # paired-end sample
#        return expand("fastq/{sample}.{group}.fastq.gz",
#                      group=[1], **wildcards)
#    # single end sample
#    return "fastq/{sample}-{unit}.fastq.gz".format(**wildcards)

#def get_fastq2(wildcards):
#    if not is_single_end(**wildcards):
#        # paired-end sample
#        return expand("fastq/{sample}.{group}.fastq.gz",
#                      group=[2], **wildcards)
#    # single end sample
#    return "fastq/{sample}.fastq.gz".format(**wildcards)


rule salmon_quant_reads:
    input:
        r1 = get_fastqs1,
        r2 = get_fastqs2,
        index = config["reference"]["salmon_index"]
    output:
        quant = "salmon/{sample}-{unit}/quant.sf",
        lib = "salmon/{sample}-{unit}/lib_format_counts.json"
    log:
        "logs/salmon/{sample}-{unit}.log"
    conda:
        "../envs/salmon.yaml"
    params:
        # path to STAR reference genome index
        libtype = "A",
        # optional parameters
        extra = ""
    threads: 12
    wrapper:
        "0.31.1/bio/salmon/quant"

rule salmon_split_counts:
    input:
        "salmon/{sample}-{unit}/quant.sf"
    output:
        "salmon/{sample}-{unit}/quant_human.sf",
        "salmon/{sample}-{unit}/quant_mouse.sf"
    script:
        "../scripts/split_organisms.py"
