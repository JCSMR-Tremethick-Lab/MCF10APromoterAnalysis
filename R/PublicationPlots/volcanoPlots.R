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
rm(results)

setwd("~/Data/Tremethick/Breast/RNA-Seq/combined/R_Analysis/")


# common plotting parameters ----------------------------------------------
cexRedDot <- 0.4
segment.alpha <- 0.5

# plot shH2A.Z data -------------------------------------------------------
targetGenes <- c("FOS", "MYC", "STAT3", "TGFB1")

N1 = nrow(rTshH2AZ[qval < 0.1 & abs(b) > 0.25])
write.csv(file = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/rTshH2AZ.csv", rTshH2AZ)
vPshH2AZ <- ggplot(data = rTshH2AZ, 
              aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", cex = 0.3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("B") + # shH2AZ
  xlab("Fold change (log2)") +
  ylab("- log10(qvalue)") +
  geom_point(data = rTshH2AZ[qval < 0.1 & abs(b) > 0.25], aes(color = "red"), cex = cexRedDot) +
  theme(legend.position = c(0.85,0.9), legend.title = element_blank(), axis.title.x = element_blank()) +
  labs(color = "q-value < 0.1 & |log2 FC| > 0.25", cex = 1) +
  scale_color_manual(labels = paste("N = ", N1, sep = ""), values = c("red")) + 
  geom_text_repel(data = rTshH2AZ[target_id %in% targetGenes], aes(label = target_id), min.segment.length = 0,size = 3, nudge_y = 75, segment.alpha = segment.alpha) +
  lims(x = c(-10,10))
vPshH2AZ
ggsave(plot = vPshH2AZ, filename = "volcanoPlot_WT_vs_shH2AZ.pdf", device = "pdf")


# plot TGFb data ----------------------------------------------------------
N2 = nrow(rTTGFb[qval < 0.1 & abs(b) > 0.25])
write.csv(file = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/rTTGFb.csv", rTTGFb)

vPTGFb <- ggplot(data = rTTGFb, 
                   aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", cex = 0.3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("A") + #TFGb 
  xlab("Fold change (log2)") +
  ylab("- log10(qvalue)") +
  geom_point(data = rTTGFb[qval < 0.1 & abs(b) > 0.25], aes(color = "red"), cex = cexRedDot) +
  theme(legend.position = c(0.85,0.9), legend.title = element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(labels = paste("N = ", N2, sep = ""), values = c("red")) +
  geom_text_repel(data = rTTGFb[target_id %in% targetGenes], aes(label = target_id), min.segment.length = 0,size = 3, nudge_y = 75, segment.alpha = segment.alpha) +
  lims(x = c(-10,10))
vPTGFb
ggsave(plot = vPTGFb, filename = "volcanoPlot_WT_vs_TGFb.pdf", device = "pdf")


# plot MCF10Ca1a data -----------------------------------------------------
load("~/Data/Tremethick/Breast/RNA-Seq/combined/R_Analysis/sleuthResults_GRCh37_hg19_UCSC_V1.rda")
rTCa1a <- as.data.table(results$sleuth_results_genes$conditionMCF10Ca1a_wt)
rTCa1a <- rTCa1a[!is.na(pval)]
write.csv(file = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/rTCa1a.csv", rTCa1a)
rm(results)

N3 = nrow(rTCa1a[qval < 0.1 & abs(b) > 0.25])
vPCa1a <- ggplot(data = rTCa1a, 
                 aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", cex = 0.3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("C") + # MCF10Ca1a
  xlab("Fold change (log2)") +
  ylab("- log10(qvalue)") +
  geom_point(data = rTCa1a[qval < 0.05 & abs(b) > 1], aes(color = "red"), cex = cexRedDot) +
  theme(legend.position = c(0.85,0.9), legend.title = element_blank()) +
  scale_color_manual(labels = paste("N = ", N3, sep = ""), values = c("red")) +
  geom_text_repel(data = rTCa1a[target_id %in% targetGenes], aes(label = target_id), min.segment.length = 0,size = 3, nudge_y = 75, segment.alpha = segment.alpha) +
  lims(x = c(-10,10))
vPCa1a
ggsave(plot = vPCa1a, filename = "volcanoPlot_WT_vs_Ca1a.pdf", device = "pdf", width = 210, height = 297, units = "mm")


combinedVPs <- gridExtra::grid.arrange(vPTGFb, vPshH2AZ, vPCa1a, ncol = 1)
ggsave(plot = combinedVPs, filename = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/volcanoPlots_combined.pdf",  width = 210, height = 297, units = "mm")
# looking for expression of some of the TFs with discovered motifs --------
tomtom <- data.table::fread("http://alternate.meme-suite.org/opal-jobs/appTOMTOM_5.0.1_1532328493389-1789556112/tomtom.tsv")


# differential expression -------------------------------------------------
volcanoPlotGenes <- c("^SP[1234]$|NFY[ABC]$|CEBPZ$|^EGR|FOS$|MYC$|STAT3")
vp1 <- ggplot(data = rTshH2AZ[qval < 0.1], 
              aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey") +
  #  geom_point(aes(color = target_id)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)
vp1 + geom_point(data = rTshH2AZ[grep(volcanoPlotGenes, target_id)], aes(color = target_id))


# expression data ---------------------------------------------------------
NFY <- ggplot(data = kTGenes[grep("D|MCF10Ca1a", sample)][grep("NFY[ABC]$|CEBPZ$", target_id)], 
             aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
NFY

SP <- ggplot(data = kTGenes[grep("D|MCF10Ca1a", sample)][grep("^SP[1234]$", target_id)], 
              aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
SP

EGR <- ggplot(data = kTGenes[grep("D|MCF10Ca1a", sample)][grep("^EGR", target_id)], 
              aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
EGR

p1 <- ggplot(data = kTGenes[grep("D|MCF10Ca1a", sample)][grep("FOS$|MYC$|STAT3|^SP1$", target_id)], 
       aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
p1

p2 <- ggplot(data = kTGenes[grep("D|MCF10Ca1a", sample)][grep("H2AFZ", target_id)], 
             aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
p2

ggplot(data = kTGenes[grep("D|MCF10Ca1a", sample)][grep("STAT3", target_id)], 
       aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  ggtitle("STAT3")

ggplot(data = kTGenes[grep("D", sample)][grep("MYC$", target_id)], 
       aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  ggtitle("MYC")

ggplot(data = kTGenes[grep("D6|MCF10Ca1a|shZD8", sample)][grep("BRCA1", target_id)], 
       aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
       
