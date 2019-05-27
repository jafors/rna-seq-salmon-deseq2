#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

print("Start")
library(tximport)
library(DESeq2)
library(rjson)
#library(AnnotationHub)
library(ensembldb)
library("biomaRt")


#print("Loaded")
#ah = AnnotationHub()
#edb<-ah[["AH57770"]] #mm10
##edb <- ah[["AH57757"]] #grch37

#Tx<-transcripts(edb, return.type = "DataFrame")
#tx2gene<-as.data.frame(Tx[,c("tx_id","gene_id")])
#gs2sym<-as.data.frame(genes(edb, return.type = "DataFrame")[,c("gene_id","symbol")])

#tx2sym<-merge(tx2gene,gs2sym,by="gene_id")[,-1]


mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = paste(snakemake@params[["species"]], "_gene_ensembl", sep = ""),
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, tx_id = ensembl_transcript_id, symbol = external_gene_name)


samplefile<-read.table(snakemake@params[["samples"]],sep='\t',header=TRUE)
sal_files<-snakemake@input[["inpt"]]
samples<- as.character(samplefile$sample)
condition<- as.character(samplefile$condition)
#line<- as.character(samplefile$line)
#batch<- as.character(samplefile$batch)

### tximport ###
print("import")
names(sal_files)<-paste0(samples)
txi.salmon_PR<-tximport(sal_files,type="salmon",txOut=F,tx2gene=t2g)

### prepare sample table ###
sampleTable_PR<-data.frame(condition=condition)#,line=line,batch=batch)
rownames(sampleTable_PR)<-colnames(txi.salmon_PR$counts)

### DESeq ###

### salmon
sal_ddsTxi_PR<-DESeqDataSetFromTximport(txi.salmon_PR,colData=sampleTable_PR,design= ~1)
sal_ddsTxi_PR<-sal_ddsTxi_PR[rowSums(counts(sal_ddsTxi_PR)) > 1,]
sal_deseq_PR<-DESeq(sal_ddsTxi_PR)
sal_counts<-counts(sal_deseq_PR, normalized=T)
raw_counts<-counts(sal_deseq_PR, normalized=F)
write.table(sal_counts, snakemake@output[[2]],sep='\t')
write.table(raw_counts, snakemake@output[[3]],sep='\t')
saveRDS(sal_deseq_PR, snakemake@output[[1]])


