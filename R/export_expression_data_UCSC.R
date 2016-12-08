require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(data.table)

ucscIDsFile <- "UCSC_IDs_mapped_to_Ensembl_transcript_IDs.rda"
if (!file.exists(ucscIDsFile)){
  ucscIDs <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                           "ucsc"),
                            mart = mart,
                            filter = "ensembl_transcript_id",
                            values = ensTranscripts$ensembl_transcript_id)
  save(ucscIDs, file = ucscIDsFile)
} else {
  load(ucscIDsFile, .GlobalEnv)
}

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#########################################
# rename colnames to resolve sample mixup
# and prepare data for export
tab1 <- resultsCompressed[[1]]$kallisto_table_wide
colnames(tab1)[2:6] <- as.character(s2c.mcf10a_vs_mcf10a_shZ$condition)[match(colnames(tab1)[2:6], as.character(s2c.mcf10a_vs_mcf10a_shZ$sample))]
colnames(tab1)[seq(2,5,2)] <- paste(colnames(tab1)[seq(2,5,2)], "rep1", sep = "_") 
colnames(tab1)[seq(3,5,2)] <- paste(colnames(tab1)[seq(3,5,2)], "rep2", sep = "_")
colnames(tab1)[6] <- paste(colnames(tab1)[6], "rep3", sep = "_")

tab2 <- resultsCompressed[[2]]$kallisto_table_wide
colnames(tab2)[2:6] <- as.character(s2c.mcf10a_vs_mcf10a_TGFb$condition)[match(colnames(tab2)[2:6], as.character(s2c.mcf10a_vs_mcf10a_TGFb$sample))]
colnames(tab2)[seq(2,5,2)] <- paste(colnames(tab2)[seq(2,5,2)], "rep1", sep = "_") 
colnames(tab2)[seq(3,5,2)] <- paste(colnames(tab2)[seq(3,5,2)], "rep2", sep = "_")
colnames(tab2)[6] <- paste(colnames(tab2)[6], "rep3", sep = "_")
colnames(tab2)[2:3] <-c("MCF10A_wt_rep3", "MCF10A_wt_rep4")

tab3 <- merge(tab1, tab2, by.x = "target_id", by.y = "target_id", all.x = T)
tab_exportFile <- paste("MCF10A_RNA-Seq_results_", runConfig$references[[refVersion]]$version[2], "_transcripts_UCSC.csv", sep = "")
write.csv(tab3, tab_exportFile)
tab3 <- data.table(tab3)
pairs(tab3[,2:11])

ucscIDs <- data.table(ucscIDs)
tab3 <- merge(tab3, ucscIDs, all.x = T, by.x = "target_id", by.y = "ensembl_transcript_id")
tab3 <- tab3[!which(is.na(tab3$ucsc)),]
tab3 <- tab3[which(tab3$ucsc != ""),]
tab3 <- tab3[,list(ucsc, MCF10A_wt_rep1, MCF10A_wt_rep2, MCF10A_wt_rep3, MCF10A_wt_rep4, 
                   MCF10A_shZ_rep1, MCF10A_shZ_rep2, MCF10A_shZ_rep3, 
                   MCF10A_TGFb_rep1, MCF10A_TGFb_rep2, MCF10A_TGFb_rep3)]
setkey(tab3, "ucsc")
tab3 <- tab3[,lapply(.SD,sum),by=ucsc]
# load("tab3.rda")
tab3 <- as.data.frame(tab3)
tab4 <- data.frame(cbind(tab3[,"ucsc"],
                         apply(tab3[,c(2:5)], 1, mean),
                         apply(tab3[,c(6:8)], 1, mean),
                         apply(tab3[,c(9:11)], 1, mean)))

colnames(tab4) <- c("ucsc", "MCF10A_WT", "MCF10A_shZ", "MCF10A_TGFb")

tab4$ucsc <- as.character(tab4$ucsc)

res <- select(txdb, keys = tab4$ucsc, columns = c("TXNAME","TXCHROM", "TXSTART", "TXEND", "TXSTRAND"), keytype = "TXNAME")
tab4 <- merge(tab4, res, by.x = "ucsc", by.y = "TXNAME")
sapply(colnames(tab4)[2:4], function(x) {
  df <- data.frame(seqnames = tab4$TXCHROM, 
                   starts = as(tab4$TXSTART, "integer"), 
                   ends = as(tab4$TXEND, "integer"), 
                   names = tab4$ucsc, 
                   scores = tab4[,x], 
                   strand = tab4$TXSTRAND)
  write.table(df, 
              file = paste("~/Data/WorkingData/", x, "_UCSC_IDs.bed", sep = ""), 
              quote = F, 
              sep = "\t", 
              row.names = F, 
              col.names = F)})


# data filtered for mean express > 1 --------------------------------------
meanExp <- apply(tab3[,2:11], 1, mean)
tab4 <- tab3[meanExp > 1,]
tab5 <- data.frame(cbind(tab4[,"ucsc"],
                         apply(tab4[,c(2:5)], 1, mean),
                         apply(tab4[,c(6:8)], 1, mean),
                         apply(tab4[,c(9:11)], 1, mean)))

colnames(tab5) <- c("ucsc", "MCF10A_WT", "MCF10A_shZ", "MCF10A_TGFb")

tab5$ucsc <- as.character(tab5$ucsc)
tab_exportFile <- paste("MCF10A_RNA-Seq_results_", runConfig$references[[refVersion]]$version[2], "_transcripts_UCSC_meanExpFiltered.csv", sep = "")
write.csv(tab5, file = tab_exportFile)

res <- select(txdb, keys = tab5$ucsc, columns = c("TXNAME","TXCHROM", "TXSTART", "TXEND", "TXSTRAND"), keytype = "TXNAME")
tab5 <- merge(tab5, res, by.x = "ucsc", by.y = "TXNAME")
sapply(colnames(tab5)[2:4], function(x) {
  df <- data.frame(seqnames = tab5$TXCHROM, 
                   starts = as(tab5$TXSTART, "integer"), 
                   ends = as(tab5$TXEND, "integer"), 
                   names = tab5$ucsc, 
                   scores = tab5[,x], 
                   strand = tab5$TXSTRAND)
  write.table(df, 
              file = paste("~/Data/WorkingData/", x, "_UCSC_IDs_meanExpFiltered.bed", sep = ""), 
              quote = F, 
              sep = "\t", 
              row.names = F, 
              col.names = F)})


