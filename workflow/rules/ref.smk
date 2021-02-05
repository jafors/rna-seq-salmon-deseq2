rule get_transcriptome:
    output:
        "resources/transcriptome.fasta"
    log:
        "logs/get-transcriptome.log"
    params:
        species=config["ref"]["species"],
        datatype="cdna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    cache: True
    wrapper:
        "0.59.2/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf"
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="" # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    cache: True
    log:
        "logs/get_annotation.log"
    wrapper:
        "0.59.2/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        "resources/genome.fasta"
    output:
        "resources/genome.fasta.fai"
    log:
        "logs/genome-faidx.log"
    cache: True
    wrapper:
        "0.59.2/bio/samtools/faidx"


rule salmon_index:
    input:
        "resources/transcriptome.fasta"
    output:
        directory("resources/salmon_index")
    log:
        "logs/salmon_index.log"
    threads: 2
    params:
        # optional parameters
        extra=""
    wrapper:
        "0.68.0/bio/salmon/index"
