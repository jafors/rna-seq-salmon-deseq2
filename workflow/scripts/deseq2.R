log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

dds<-readRDS(snakemake@input[[1]])

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
print(dds)
contrast <- c("condition", snakemake@params[["contrast"]])
res<-results(dds, contrast=contrast)
# shrink fc for lowly expressed genes
res<-lfcShrink(dds, contrast=contrast)
# sort by p_value
res<-res[order(res$padj),]

dds2 <- DESeq(dds, test="LRT", reduced=snakemake@params[["reduced_model"]])
res2 <- results(dds2)

# store results
svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(data.frame("ID" = rownames(res), res) ,file=snakemake@output[["table"]], row.names=F, sep='\t')
write.table(data.frame("ID" = rownames(res2), res2) ,file=snakemake@output[["table_LRT"]], row.names=F, sep='\t')
