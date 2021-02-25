rule salmon_quant_reads_PE:
    input:
        r1=get_salmon_input_R1,
        r2=get_salmon_input_R2,
        index="resources/salmon_index"
    output:
        quant="results/salmon/{sample}/quant.sf",
        lib="results/salmon/{sample}/lib_format_counts.json"
    log:
        "logs/salmon/{sample}.log"
    params:
        # path to STAR reference genome index
        libtype="A",
        # optional parameters
        extra=""
    threads: 2
    wrapper:
        "0.68.0/bio/salmon/quant"
        
rule salmon_quant_reads_SE:
    input:
        r=get_salmon_input_R1,
        index="resources/salmon_index"
    output:
        quant="results/salmon_SE/{sample}/quant.sf",
        lib="results/salmon_SE/{sample}/lib_format_counts.json"
    log:
        "logs/salmon/{sample}.log"
    params:
        # path to STAR reference genome index
        libtype="A",
        # optional parameters
        extra=""
    threads: 2
    wrapper:
        "0.68.0/bio/salmon/quant"

