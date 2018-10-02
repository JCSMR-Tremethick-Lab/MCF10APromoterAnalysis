library(data.table)
library(ggplot2)

load("~/Data/Tremethick/Breast/RNA-Seq/combined/R_Analysis/sleuthResults_GRCh37_hg19_UCSC_V1.rda")
setwd("~/Data/Tremethick/Breast/RNA-Seq/combined/R_Analysis/")

lapply(results$sleuth_results_genes, function(x) {
  setkey(x, "target_id")
  fc <- round(x["H2AFZ"]$b, 2)
  perc <- round(2^fc, 2)
  y <- list("fold-change" = fc, "percent reduction" = perc)
  return(y)
})

kTGenes <- results$kallisto_table_genes
setkey(kTGenes, "target_id")

rTshH2AZ <- results$sleuth_results_genes$conditionMCF10A_shZ
N = nrow(rTshH2AZ[qval < 0.05 & abs(b) > 1])
vPshH2AZ <- ggplot(data = rTshH2AZ, 
              aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", cex = 0.3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("Differential expression - MCF10A vs MCF10A + shH2AZ") +
  xlab("Fold change (log2)") +
  ylab("- log10(qvalue)") +
  geom_point(data = rTshH2AZ[qval < 0.05 & abs(b) > 1], aes(color = "red"), cex = 0.6) +
  theme(legend.position = c(0.85,0.9)) +
  labs(color = "q-value < 0.05 and absolute log2 FC > 1") +
  scale_color_manual(labels = paste("N = ", N, sep = ""), values = c("red"))
vPshH2AZ
ggsave(plot = vPshH2AZ, filename = "volcanoPlot_WT_vs_shH2AZ.pdf", device = "pdf")

rTTGFb <- results$sleuth_results_genes$conditionMCF10A_TGFb
N = nrow(rTTGFb[qval < 0.05 & abs(b) > 1])
vPTGFb <- ggplot(data = rTTGFb, 
                   aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", cex = 0.3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("Differential expression - MCF10A vs MCF10A + TGFb")+ 
  xlab("Fold change (log2)") +
  ylab("- log10(qvalue)") +
  geom_point(data = rTTGFb[qval < 0.05 & abs(b) > 1], aes(color = "red"), cex = 0.6) +
  theme(legend.position = c(0.85,0.9)) +
  labs(color = "q-value < 0.05 and absolute log2 FC > 1") +
  scale_color_manual(labels = paste("N = ", N, sep = ""), values = c("red"))
vPTGFb
ggsave(plot = vPTGFb, filename = "volcanoPlot_WT_vs_TGFb.pdf", device = "pdf")

rTCa1a <- results$sleuth_results_genes$conditionMCF10Ca1a_wt
N = nrow(rTCa1a[qval < 0.05 & abs(b) > 1])
vPCa1a <- ggplot(data = rTCa1a, 
                 aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", cex = 0.3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("Differential expression - MCF10A vs MCF10Ca1a") + 
  xlab("Fold change (log2)") +
  ylab("- log10(qvalue)") +
  geom_point(data = rTCa1a[qval < 0.05 & abs(b) > 1], aes(color = "red"), cex = 0.6) +
  theme(legend.position = c(0.85,0.9)) +
  labs(color = "q-value < 0.05 and absolute log2 FC > 1") +
  scale_color_manual(labels = paste("N = ", N, sep = ""), values = c("red"))
vPCa1a
ggsave(plot = vPCa1a, filename = "volcanoPlot_WT_vs_Ca1a.pdf", device = "pdf")

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
       