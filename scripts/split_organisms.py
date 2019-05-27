import pandas as pd

table = pd.read_csv(snakemake.input[0],sep='\t')

human = table[table.Name.str.contains("ENST")]
mouse = table[table.Name.str.contains("ENSMUST")]

human.to_csv(snakemake.output[0], index=False, sep = '\t')
mouse.to_csv(snakemake.output[1], index=False, sep = '\t')
