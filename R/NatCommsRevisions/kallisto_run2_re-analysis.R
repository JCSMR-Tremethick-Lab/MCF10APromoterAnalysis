library(sleuth)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

setwd("/home/sebastian/Development/JCSMR-Tremethick-Lab/MCF10APromoterAnalysis/R")

# prepare annotation
TxDbUCSC <- TxDb.Hsapiens.UCSC.hg19.knownGene
knownGenes <- genes(TxDbUCSC)
knownExons <- exonsBy(TxDbUCSC, by = "gene")
knownTranscripts <- transcriptsBy(TxDbUCSC, by = 'gene')

promoterIDs <- readLines('Promoters_UCSC_IDs.txt')
head(promoterIDs)
res <- select(TxDbUCSC, promoterIDs, columns = c('GENEID', 'TXNAME'), keytype = 'TXNAME')
res2 <- select(org.Hs.eg.db, res[!is.na(res$GENEID), 'GENEID'], columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'ENTREZID')
head(res2)
dim(res2)
dim(res)
res3 <- merge(res, res2, by.x = 'GENEID', by.y = 'ENTREZID', all.x = T, all.y = F)
head(res3)
dim(res3)
res3 <- res3[isUnique(res3$TXNAME),]
head(res3)
t2g <- res3
t2g <- data.table::as.data.table(t2g)
setnames(t2g, 'TXNAME', 'target_id')

kal_dir <- '~/Data/Tremethick/Breast/RNA-Seq_run2/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/processed_data/GRCh37_hg19_UCSC_ERCC/kallisto'
list.dirs(kal_dir, full.names = F)
sample_id <- list.dirs(kal_dir, full.names = F)[-1]
kal_dirs <- paste(kal_dir, sample_id, sep = '/')

s2c <- data.table::data.table(sample = sample_id, 
                              condition = c(rep('WT', 6), rep('shZ', 3), rep('TGFb', 3)), 
                              day = c(rep('D6', 3), rep('D8', 6), rep('D6', 3)))
s2c$path <- kal_dirs
s2c$condition <- as.factor(s2c$condition)
s2c$condition <- relevel(s2c$condition, ref = "WT")


design <- model.matrix(~ condition, data = s2c[s2c$day == 'D8'])

so <- sleuth_prep(s2c[s2c$day == 'D8'], extra_bootstrap_summary = TRUE, target_mapping = t2g, aggregation_column = 'SYMBOL')
so <- sleuth::sleuth_fit(so, formula = design[,c(1:2)])
so <- sleuth::sleuth_fit(so, ~1, "reduced")
so <- sleuth::sleuth_lrt(so, "reduced", "full")
for (i in colnames(design)[grep("Intercept", colnames(design[,c(1:2)]), invert = T)]){
  so <- sleuth::sleuth_wt(so, i)  
}

