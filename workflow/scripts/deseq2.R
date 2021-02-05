log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

dds<-readRDS(snakemake@input[[1]])

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

contrast <- c("condition", snakemake@params[["contrast"]])
res<-results(dds, contrast=contrast)
# shrink fc for lowly expressed genes
res<-lfcShrink(dds, contrast=contrast)
# sort by p_value
res<-res[order(res$padj),]

# store results
svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

write.csv(as.data.frame(res),file=snakemake@output[["table"]])

