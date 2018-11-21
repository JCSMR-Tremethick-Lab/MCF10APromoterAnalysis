# volcano plots
library(data.table)
library(ggplot2)
library(sleuth)
library(ggrepel)

# import data
sigEMTCellLines <- data.table::fread("https://raw.githubusercontent.com/skurscheid/GeneSets/master/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_cellLine.txt")
sigEMTTumors <- data.table::fread("https://raw.githubusercontent.com/skurscheid/GeneSets/master/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_tumor.txt")

# load pre-processed expression data
load("~/Data/Tremethick/Breast/RNA-Seq_run2/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/R_Analysis/sleuthResults_GRCh37_hg19_ensembl75_ERCC_V2.rda")

rTshH2AZ <- as.data.table(results$MCF10A_vs_shZ$sleuth_results_genes$conditionMCF10A_shZ)
setkey(rTshH2AZ, target_id)
rTshH2AZ <- rTshH2AZ[!is.na(pval)]
rTTGFb <- as.data.table(results$MCF10A_vs_TGFb$sleuth_results_genes$conditionMCF10A_TGFb)
setkey(rTTGFb, target_id)
rTTGFb <- rTTGFb[!is.na(pval)]

setkey(sigEMTCellLines, cellLine_sig)
setkey(sigEMTTumors, tumor_sig)
sigEMTCellLines$source <- "cells"
sigEMTTumors$source <- "tumors"

colnames(sigEMTTumors)[1] <- "gene_symbol"
colnames(sigEMTCellLines)[1] <- "gene_symbol"

sigEMTTumors[intersect(sigEMTCellLines$gene_symbol, sigEMTTumors$gene_symbol)]
sigEMTCellLines[intersect(sigEMTCellLines$gene_symbol, sigEMTTumors$gene_symbol)]

sigEMT <- rbind(sigEMTTumors,
                sigEMTCellLines[!intersect(sigEMTCellLines$gene_symbol, sigEMTTumors$gene_symbol)])
sigEMT <- rbind(sigEMT, list("TGFB1", "mes", "manual"))
sigEMT <- rbind(sigEMT, list("H2AFZ", "mes", "manual"))

rTshH2AZ[sigEMT$gene_symbol][!is.na(pval)]
rTshH2AZ$epi_mes <- NA

merge(rTshH2AZ, sigEMT, by.x = "target_id", by.y = "gene_symbol")[target_id == "H2AFZ"]


# volcano plot for WT vs shH2AZ -------------------------------------------
labelDatshH2AZ <- rbind(merge(rTshH2AZ, 
                            sigEMT, 
                            by.x = "target_id", 
                            by.y = "gene_symbol"))
setkey(labelDatshH2AZ, target_id)

volcPlotEMTshH2AZ <- ggplot(data = rTshH2AZ, 
                            aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", size = 0.8) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)
volcPlotEMTshH2AZ <- volcPlotEMTshH2AZ + geom_point(data = labelDatshH2AZ, aes(colour = epi_mes), size = 0.8)
volcPlotEMTshH2AZ <- volcPlotEMTshH2AZ + geom_text_repel(data = rbind(labelDatshH2AZ[qval < 0.1][abs(b) > 1.5], labelDatshH2AZ[target_id %in% c("H2AFZ", "TGFB1")]), 
                                  aes(label = target_id), 
                                  min.segment.length = 0,
                                  size = 3)
volcPlotEMTshH2AZ
# highlight cfos, cmy, stat3
targetGenes <- c("FOS", "MYC", "STAT3")

volcPlotEMTshH2AZ <- volcPlotEMTshH2AZ + geom_text_repel(data = rTshH2AZ[target_id %in% targetGenes], 
                                                         aes(label = target_id), 
                                                         min.segment.length = 0,
                                                         size = 3)


volcPlotEMTshH2AZ
ggsave("~/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/DifferentialExpression_shH2AZ_volcanoPlot.pdf", volcPlotEMTshH2AZ, scale = 2)
write.csv(rTTGFb[order(qval, decreasing = F)], file = "~/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/DifferentialExpression_shH2AZ_Table.csv")

# volcano plot for WT vs TGFB ---------------------------------------------
rm(labelDatTGFb)
labelDatTGFb <- rbind(merge(rTTGFb, 
                            sigEMT, 
                            by.x = "target_id", 
                            by.y = "gene_symbol"))
setkey(labelDatTGFb, target_id)



volcPlotEMTTGFb <- ggplot(data = rTTGFb,
                          aes(x = b, y = -(log10(qval)))) +
                    geom_point(color = "grey", size = 0.8) +
                    geom_vline(xintercept = 0) +
                    geom_hline(yintercept = 0)
volcPlotEMTTGFb <- volcPlotEMTTGFb + geom_point(data = labelDatTGFb, aes(colour = epi_mes), size = 1.3)
volcPlotEMTTGFb <- volcPlotEMTTGFb + geom_text_repel(data = rbind(labelDatTGFb[qval < 0.1][abs(b) > 2.75], labelDatTGFb[target_id %in% c("H2AFZ", "TGFB1")]), 
                                                     aes(label = target_id), 
                                                     min.segment.length = 0,
                                                     size = 3)


volcPlotEMTTGFb <- volcPlotEMTTGFb + geom_text_repel(data = rTTGFb[target_id %in% targetGenes], 
                                                         aes(label = target_id), 
                                                         min.segment.length = 0,
                                                         size = 3)
volcPlotEMTTGFb
ggsave("~/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/DifferentialExpression_TGFb_volcanoPlot.pdf", volcPlotEMTTGFb, scale = 2)
write.csv(rTTGFb[order(qval, decreasing = F)], file = "~/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/DifferentialExpression_TGFb_Table.csv")


