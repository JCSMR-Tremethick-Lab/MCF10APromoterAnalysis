require(TxDb.Hsapiens.UCSC.hg19.knownGene)
ucscIDs <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                         "ucsc"),
                          mart = mart,
                          filter = "ensembl_transcript_id",
                          values = ensTranscripts$ensembl_transcript_id)

#########################################
# rename colnames to resolve sample mixup
# and prepare data for export
tab1 <- resultsCompressed[[1]]$kallisto_table_wide
colnames(tab1)[2:9] <- as.character(s2c.mcf10a$condition)[match(colnames(tab1)[2:9], as.character(s2c.mcf10a$sample))]
tab1 <- tab1[, c(1,3,6,4,5,2,7:9)]
colnames(tab1)[seq(2,9,2)] <- paste(colnames(tab1)[seq(2,9,2)], "rep1", sep = "_") 
colnames(tab1)[seq(3,9,2)] <- paste(colnames(tab1)[seq(3,9,2)], "rep2", sep = "_")

tab2 <- resultsCompressed[[2]]$kallisto_table_wide
# colnames(tab2)[2:5] <- as.character(s2c.mcf10Ca1a$condition)[match(colnames(tab2)[2:5], as.character(s2c.mcf10Ca1a$sample))]
# colnames(tab2)[seq(2,5,2)] <- paste(colnames(tab2)[seq(2,5,2)], "rep1", sep = "_") 
# colnames(tab2)[seq(3,5,2)] <- paste(colnames(tab2)[seq(3,5,2)], "rep2", sep = "_")

tab3 <- merge(tab1, tab2[, c(1, 4:5)], by.x = "target_id", by.y = "target_id", all.x = T)
tab_exportFile <- paste("MCF10A_RNA-Seq_results_", runConfig$references[[refVersion]]$version, "_transcripts_UCSC.csv", sep = "")
write.csv(tab3, tab_exportFile)
tab3 <- data.table(tab3)
ucscIDs <- data.table(ucscIDs)
tab3 <- merge(tab3, ucscIDs, all.x = T, by.x = "target_id", by.y = "ensembl_transcript_id")
tab3 <- tab3[!which(is.na(tab3$ucsc)),]
tab3 <- tab3[which(tab3$ucsc != ""),]
tab3 <- tab3[,list(ucsc, MCF10A_wt_rep1, MCF10A_wt_rep2, MCF10A_shZ_rep1, MCF10A_shZ_rep2, MCF10A_TGFb_rep1, MCF10A_TGFb_rep2, MCF10Ca1a_wt_rep1, MCF10Ca1a_wt_rep2, MCF10Ca1a_shZ_rep1, MCF10Ca1a_shZ_rep2)]
setkey(tab3, "ucsc")
tab3 <- tab3[,lapply(.SD,sum),by=ucsc]

# load("tab3.rda")
tab3 <- as.data.frame(tab3)
tab4 <- data.frame(cbind(tab3[,"ucsc"],
                         apply(tab3[,c(2,3)], 1, mean),
                         apply(tab3[,c(4,5)], 1, mean),
                         apply(tab3[,c(6,7)], 1, mean),
                         apply(tab3[,c(8,9)], 1, mean),
                         apply(tab3[,c(10,11)], 1, mean)))

tab4$ucsc <- as.character(tab4$ucsc)
colnames(tab4)[1] <- "ucsc"
colnames(tab4)[2:6] <- unlist(lapply(strsplit(colnames(tab3)[seq(2,11,2)], "_"), function(x) paste(x[1:2], collapse = "_")))
res <- select(txdb, keys = tab4$ucsc, columns = c("TXNAME","TXCHROM", "TXSTART", "TXEND", "TXSTRAND"), keytype = "TXNAME")
tab4 <- merge(tab4, res, by.x = "ucsc", by.y = "TXNAME")
sapply(colnames(tab4)[2:6], function(x) {
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
