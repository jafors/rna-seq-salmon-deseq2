#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

print("Start")
library(tximport)
library(DESeq2)
library(ensembldb)
library("biomaRt")
library("tidyverse")

# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because 
# of very frequent "Internal Server Error"s)
mart <- "useast"
rounds <- 0
while ( class(mart)[[1]] != "Mart" ) {
  mart <- tryCatch(
    {
      # done here, because error function does not
      # modify outer scope variables, I tried
      if (mart == "www") rounds <- rounds + 1
      # equivalent to useMart, but you can choose
      # the mirror instead of specifying a host
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = str_c(snakemake@params[["species"]], "_gene_ensembl"),
        mirror = mart
      )
    },
    error = function(e) {
      # change or make configurable if you want more or
      # less rounds of tries of all the mirrors
      if (rounds >= 3) {
        stop(
          str_c(
            "Have tried all 4 available Ensembl biomaRt mirrors ",
            rounds,
            " times. You might have a connection problem, or no mirror is responsive."
          )
        )
      }
      # hop to next mirror
      mart <- switch(mart,
                     useast = "uswest",
                     uswest = "asia",
                     asia = "www",
                     www = {
                       # wait before starting another round through the mirrors,
                       # hoping that intermittent problems disappear
                       Sys.sleep(30)
                       "useast"
                     }
              )
    }
  )
}


t2g <- biomaRt::getBM(
            attributes = c( "ensembl_transcript_id",
                            "ensembl_gene_id",
                            "external_gene_name"),
            mart = mart,
            ) %>%
        rename( target_id = ensembl_transcript_id,
                ens_gene = ensembl_gene_id,
                ext_gene = external_gene_name
                )


#t2g <- dplyr::rename(t2g, tx_id = ensembl_transcript_id, symbol = external_gene_name)


samplefile<-read.table(snakemake@params[["samples"]],sep='\t',header=TRUE)
print(samplefile)
sal_files<-snakemake@input[["inpt"]]
print(sal_files)
samples<- as.character(samplefile$sample_name)
condition<- as.character(samplefile$condition)
#line<- as.character(samplefile$line)
batch<- as.character(samplefile$batch)
print(batch)

### tximport ###
print("import")
names(sal_files)<-paste0(samples)
txi.salmon_PR<-tximport(sal_files,type="salmon",txOut=F,tx2gene=t2g, ignoreTxVersion=TRUE)

### prepare sample table ###
sampleTable_PR<-data.frame(condition=condition, batch=batch)
rownames(sampleTable_PR)<-colnames(txi.salmon_PR$counts)


### DESeq ###

formula = as.formula("~condition")

### salmon
sal_ddsTxi_PR<-DESeqDataSetFromTximport(txi.salmon_PR,colData=sampleTable_PR,design= formula)
sal_ddsTxi_PR<-sal_ddsTxi_PR[rowSums(counts(sal_ddsTxi_PR)) > 1,]
sal_deseq_PR<-DESeq(sal_ddsTxi_PR)
sal_counts<-counts(sal_deseq_PR, normalized=T)
raw_counts<-counts(sal_deseq_PR, normalized=F)
write.table(sal_counts, snakemake@output[[2]],sep='\t')
write.table(raw_counts, snakemake@output[[3]],sep='\t')
saveRDS(sal_deseq_PR, snakemake@output[[1]])


