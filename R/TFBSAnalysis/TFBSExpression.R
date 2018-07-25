#
library("sleuth")
library("Homo.sapiens")
library("TFBSTools")
library("JASPAR2014")
library("TFutils")
library("data.table")

Hsap <- Homo.sapiens
analysisDir <- "/home/sebastian/Data/Tremethick/Breast/RNA-Seq/combined/R_Analysis"
setwd(analysisDir)

# load gene expression results
load("/home/sebastian/Data/Tremethick/Breast/RNA-Seq/combined/R_Analysis/sleuthResults_GRCh37_hg19_UCSC_V1.rda")

kTGenes <- results$kallisto_table_genes
setkey(kTGenes, "target_id")

# GO:0003700 is top parent term for "DNA binding transcription factor activity"
putativeTFs <- data.table(select(Hsap, "GO:0003700", columns = c("SYMBOL", "GO"), keytype = "GO"))
setkey(putativeTFs, "SYMBOL")
putativeTFs <- putativeTFs[,c("SYMBOL", "GO")]
putativeTFs <- putativeTFs[!duplicated(SYMBOL)]
expressedTFs <- putativeTFs[which(putativeTFs$SYMBOL %in% results$sleuth_object_genes$obs_norm$target_id)]$SYMBOL
length(expressedTFs)
DTexpressedTFs <- kTGenes[expressedTFs]
# search in JASPAR for TFBS data
opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- expressedTFs[1:915]
opts[["all_versions"]] <- T
PFMatrixList <- getMatrixSet(JASPAR2014, opts)
