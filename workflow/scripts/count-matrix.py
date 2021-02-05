import pandas as pd

counts = [pd.read_csv(f, index_col=0, usecols=[0, 4], header=None, skiprows=[0], sep='\t')
          for f in snakemake.input]

for t, (sample, unit) in zip(counts, snakemake.params.units.index):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "TranscriptID"
matrix=matrix.reset_index()
matrix["TranscriptID"] = matrix["TranscriptID"].str.split(".",expand=True)[0]
#matrix[["GeneID","Symbol","TranscriptID"]] = matrix["TranscriptID"].str.split("|", expand=True)
#matrix.to_csv(snakemake.output[0], sep="\t")
t2g = pd.read_csv(snakemake.params["transcript_map"], sep='\t')

matrix = matrix.merge(t2g, on="TranscriptID", how="left").drop("TranscriptID", axis=1)
## collapse technical replicates
matrix = matrix.groupby(["GeneID","Symbol"]).sum()
#print(matrix)
matrix.to_csv(snakemake.output[0], sep="\t")
