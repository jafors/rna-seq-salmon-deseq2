def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_fastqs1(wildcards):
    """Get raw FASTQ files from unit sheet."""
    return units.loc[
        (wildcards.sample, wildcards.unit), ["fq1"]].dropna()

def get_fastqs2(wildcards):
    """Get raw FASTQ files from unit sheet."""
    return units.loc[
        (wildcards.sample, wildcards.unit), ["fq2"]].dropna()
