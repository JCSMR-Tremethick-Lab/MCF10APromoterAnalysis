load("~/Data/Tremethick/Breast/RNA-Seq/combined/R_Analysis/sleuthResults_GRCh37_hg19_UCSC_V1.rda")

lapply(results$sleuth_results_genes, function(x) {
  setkey(x, "target_id")
  fc <- round(x["H2AFZ"]$b, 2)
  perc <- round(2^fc, 2)
  y <- list("fold-change" = fc, "percent reduction" = perc)
  return(y)
})

kTGenes <- results$kallisto_table_genes
setkey(kTGenes, "target_id")

kTGenes[grep("D|MCF10Ca1a", sample)][grep("^SP1$", target_id)]
rTshH2AZ <- results$sleuth_results_genes$conditionMCF10A_shZ

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
       