## load libraries
library(ggparallel)
library(alluvial)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(msigdbr)
library(bitr)
library(sleuth)

# remove all variables from environment
rm(list = ls())

dataDir <- "./alluvial_plots_sensitivity_h2az/"

# load expression data ----------------------------------------------------
load(file = file.path(dataDir, 'so_wt_ca1a.rda'))
rT <- sleuth::sleuth_results(so, test = 'conditionca1a', test_type = 'wt', show_all = F)
setDT(rT)
setkey(rT, target_id)  

kT <- sleuth::kallisto_table(so, normalized = T, include_covariates = T)                  
setDT(kT)
kTMean <- kT[, lapply(.SD, mean), by = list(condition, target_id), .SDcols='tpm']
kTMeanWide <- dcast(kTMean, target_id ~ condition, value.var = 'tpm')

# load parallel plot data -------------------------------------------------
l1 <- lapply(list.files(path = dataDir, pattern="Log2"), function(x){
  dt1 <- data.table::fread(paste(dataDir, x, sep = "/"))
  return(dt1)
})
n1 <- gsub(x = list.files(path = dataDir, pattern="Log2"), pattern = ".tsv", replacement =  "") 
names(l1) <- unlist(lapply(strsplit(n1, "_"), function(x) paste(x[2:3], collapse = "_")))

l1 <- lapply(l1, function(x) {
  ucscID <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[4]))
  extGene <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[6]))
  x$ucscID <- ucscID
  x$extGene <- extGene
  x$chr <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[1]))
  x$start <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[2]))
  x$end <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[3]))
  return(x)
})
l1

# For h2az (Figure X - sensitivity)
geneList <- unique(l1$WT_H2AZ$extGene)

lapply(names(l1), function(x){
  setkey(l1[[x]], "extGene")
})
mcf10awtCategories <- unique(l1$WT_H2AZ[,c('group1', 'color')])
setorder(mcf10awtCategories, group1)

# merge all tables into single data.table ---------------------------------
dt1 <- l1$WT_H2AZ[!duplicated(extGene), c("extGene" ,"group1")]
colnames(dt1)[2] <- "wt.group"
dt1 <- merge(dt1, l1$TFGb_H2AZ[!duplicated(extGene), c("extGene", "group1")], by.x = "extGene", by.y = "extGene")
colnames(dt1)[length(colnames(dt1))] <- "tgfb.group"
dt1 <- merge(dt1, l1$CA1a_H2AZ[!duplicated(extGene), c("extGene","group1")], by.x = "extGene", by.y = "extGene")
colnames(dt1)[length(colnames(dt1))] <- "ca1a.group"
dt1
save(dt1, file = file.path(dataDir, "alluvialPlotsSensitivityH2AZData.rda"))

tab.ca1a <- table("WT" = dt1$wt.group, "Ca1a" = dt1$ca1a.group)
rownames(tab.ca1a) <- paste('WT_H2AZ_cluster_', c(1:7), sep = '')
colnames(tab.ca1a) <- paste('Ca1a_H2AZ_cluster_', c(1:7), sep = '')
write.csv(tab.ca1a, file = file.path(dataDir, "cross_table_WT_Ca1a.csv"))

# alluvial plots WT -> CA1A 5,7 to 1,3,7---------------------
# 5,7 to 1,3,7
fig_wt_ca1a <- data.table::as.data.table(table("WT" = dt1$wt.group, "Ca1a" = dt1$ca1a.group))
fig_wt_ca1a %>% group_by(WT, Ca1a) %>% summarise(n = sum(N)) -> fig_wt_ca1a

mcf10awtCategories$color4 <- mcf10awtCategories$color
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(5,7)]$color4 <- "grey"
mcf10awtCategories$alpha4 <- 1
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(5,7)]$alpha4 <- 0.1


png(file = file.path(dataDir, "alluvial_plot_sensitivity_h2az_wt_5_7_ca1a_1_3_7.png"))
alluvial(fig_wt_ca1a[,c(1:2)], freq = fig_wt_ca1a$n,
         col = mcf10awtCategories$color4[match(as.integer(fig_wt_ca1a$WT), mcf10awtCategories$group)],
         hide = !fig_wt_ca1a$Ca1a %in% c(1,3,7),
         alpha = mcf10awtCategories$alpha4[match(as.integer(fig_wt_ca1a$WT), mcf10awtCategories$group)])
dev.off()

pdf(file = file.path(dataDir, "alluvial_plot_sensitivity_h2az_wt_5_7_ca1a_1_3_7.pdf"))
alluvial(fig_wt_ca1a[,c(1:2)], freq = fig_wt_ca1a$n,
         col = mcf10awtCategories$color4[match(as.integer(fig_wt_ca1a$WT), mcf10awtCategories$group)],
         hide = !fig_wt_ca1a$Ca1a %in% c(1,3,7),
         alpha = mcf10awtCategories$alpha4[match(as.integer(fig_wt_ca1a$WT), mcf10awtCategories$group)])
dev.off()

selectedGenes1 <- dt1[dt1$ca1a.group %in% c(1,3,7) & (dt1$wt.group %in% c(5,7))]$extGene # list of genes
geneTable1 <- AnnotationDbi::select(org.Hs.eg.db, keys = selectedGenes1, keytype = 'SYMBOL', columns = c('SYMBOL', 'GENENAME'))
setDT(geneTable1)
geneTable1 <- merge(dt1[,.(extGene, wt.group, ca1a.group)], geneTable1, by.x = 'extGene', by.y = 'SYMBOL')
geneTable1 <- merge(geneTable1, kTMeanWide, by.x = 'extGene', by.y = 'target_id', all.x = T, all.y = F)
geneTable1 <- merge(geneTable1, rT[,.(target_id, b, qval)], by.x = 'extGene', by.y = 'target_id', all.x = T, all.y = F)

table(geneTable1$ca1a.group)
table(geneTable1$wt.group)
write.csv(geneTable1, file = file.path(dataDir, 'gene_table_sensitivity_h2az_wt_5_7_ca1a_1_3_7.csv'))
write.csv(geneTable1[(qval <= 0.1 & b > 0)], file = file.path(dataDir, 'gene_table_sensitivity_h2az_wt_5_7_ca1a_1_3_7_upregulated.csv'))

# differential expression of these
rT1 <- rT[selectedGenes1]
rT1 <- rT1[!is.na(pval)]
write.csv(rT1, file = file.path(dataDir, 'diffGenes_wt_5_7_ca1a_1_3_7.csv'))
hist(rT1$qval)
plot(rT1$b, -log10(rT1$qval))
setkey(kT, target_id)
kT1 <- kT[selectedGenes1]
kT1 <- kT1[!is.na(kT1$tpm)]
kT1 <- merge(kT1[kT1$condition %in% c('WT', 'ca1a')], dt1[,c('extGene','wt.group','ca1a.group')], by.x = 'target_id', by.y = 'extGene', all.x = T, all.y = F)
table(kT1$ca1a.group)
table(kT1$wt.group)

bpWT1 <- ggplot2::ggplot(kT1, aes(x = as.character(wt.group), y = log2(tpm + 1), group = wt.group)) + 
  geom_boxplot() 
vpWT1 <- ggplot2::ggplot(kT1, aes(x = as.character(wt.group), y = log2(tpm + 1), group = wt.group)) + 
  geom_violin() 
ggsave(bpWT1, filename = file.path(dataDir, 'boxplot_sensitivity_h2az_wt_5_7.png'))

bpca1a1 <- ggplot2::ggplot(kT1, aes(x = as.character(ca1a.group), y = log2(tpm + 1), group = ca1a.group)) + 
  geom_boxplot()
vptca1a1 <- ggplot2::ggplot(kT1, aes(x = as.character(ca1a.group), y = log2(tpm + 1), group = ca1a.group)) + 
  geom_violin()
ggsave(bpca1a1, filename = file.path(dataDir, 'boxplot_sensitivity_h2az_ca1a_1_3_7.png'))

# enrichment analysis
m_t2g <- msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name, human_gene_symbol)

geneList <- rT[order(b,decreasing=TRUE),]$b
names(geneList) <- rT[order(b,decreasing=TRUE),]$target_id
gene <- geneList[abs(geneList) > 0.5]

em1 <- enricher(rT1[rT1$qval < 0.1 & rT1$b > 0]$target_id, pvalueCutoff = 0.01, TERM2GENE = m_t2g)
png(file = file.path(dataDir, "msigdb_dotplot_sensitivity_h2az_wt_5_7_ca1a_1_3_7_upregulated.png"), height = 1024, width = 768)
dotplot(em1)
dev.off()

pdf(file = file.path(dataDir, "msigdb_dotplot_sensitivity_h2az_wt_5_7_ca1a_1_3_7_upregulated.pdf"), height = 15, width = 20)
dotplot(em1)
dev.off()

write.csv(rT1[rT1$qval < 0.1 & rT1$b > 0, .(target_id, b, qval, mean_obs),][order(b, decreasing = T)], 
          file = file.path(dataDir, 'upregulated_genes_sensitivity_h2az_wt_5_7_ca1a_1_3_7.csv'))

# alluvial plots WT -> CA1A 1,3,4,6 to 4,6---------------------------------
# 1,3,4,6 to 4,6
mcf10awtCategories$color5 <- mcf10awtCategories$color
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(1,3,4,6)]$color5 <- "grey"
mcf10awtCategories$alpha5 <- 1
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(1,3,4,6)]$alpha5 <- 0.1

png(file = file.path(dataDir, "alluvial_plot_sensitivity_h2az_wt_1_3_4_6_ca1a_4_6.png"))
alluvial(fig_wt_ca1a[,c(1:2)], freq = fig_wt_ca1a$n,
         col = mcf10awtCategories$color5[match(as.integer(fig_wt_ca1a$WT), mcf10awtCategories$group)],
         hide = !fig_wt_ca1a$Ca1a %in% c(4,6),
         alpha = mcf10awtCategories$alpha5[match(as.integer(fig_wt_ca1a$WT), mcf10awtCategories$group)])
dev.off()

pdf(file = file.path(dataDir, "alluvial_plot_sensitivity_h2az_wt_1_3_4_6_ca1a_4_6.pdf"))
alluvial(fig_wt_ca1a[,c(1:2)], freq = fig_wt_ca1a$n,
         col = mcf10awtCategories$color5[match(as.integer(fig_wt_ca1a$WT), mcf10awtCategories$group)],
         hide = !fig_wt_ca1a$Ca1a %in% c(4,6),
         alpha = mcf10awtCategories$alpha5[match(as.integer(fig_wt_ca1a$WT), mcf10awtCategories$group)])
dev.off()

selectedGenes2 <- dt1[dt1$ca1a.group %in% c(4,6) & (dt1$wt.group %in% c(1,3,4,6))]$extGene # list of genes
geneTable2 <- AnnotationDbi::select(org.Hs.eg.db, keys = selectedGenes2, keytype = 'SYMBOL', columns = c('SYMBOL', 'GENENAME'))
setDT(geneTable2)
geneTable2 <- merge(dt1[,.(extGene, wt.group, ca1a.group)], geneTable2, by.x = 'extGene', by.y = 'SYMBOL')
geneTable2 <- merge(geneTable2, kTMeanWide, by.x = 'extGene', by.y = 'target_id', all.x = T, all.y = F)
geneTable2 <- merge(geneTable2, rT[,.(target_id, b, qval)], by.x = 'extGene', by.y = 'target_id', all.x = T, all.y = F)
table(geneTable2$wt.group)
table(geneTable2$ca1a.group)
write.csv(geneTable2, file = file.path(dataDir, 'gene_table_sensitivity_h2az_wt_1_3_4_6_ca1a_4_6.csv'))
write.csv(geneTable2[(qval <= 0.1 & b < 0)], file = file.path(dataDir, 'gene_table_sensitivity_h2az_wt_1_3_4_6_ca1a_4_6_downregulated.csv'))

# differential expression of these
rT2 <- rT[selectedGenes2]
rT2 <- rT2[!is.na(pval)]
write.csv(rT2, file = file.path(dataDir, 'diffGenes_wt_1_3_4_6_ca1a_4_6.csv'))
hist(rT2$qval)
plot(rT2$b, -log10(rT2$qval))
setkey(kT, target_id)
kT2 <- kT[selectedGenes2]
kT2 <- kT2[!is.na(kT2$tpm)]
kT2 <- merge(kT2[kT2$condition %in% c('WT', 'ca1a')], dt1[,c('extGene','wt.group','ca1a.group')], by.x = 'target_id', by.y = 'extGene', all.x = T, all.y = F)
table(kT2$ca1a.group)
table(kT2$wt.group)

bpWT2 <- ggplot2::ggplot(kT2, aes(x = as.character(wt.group), y = log2(tpm + 1), group = wt.group)) + 
  geom_boxplot() 
vpWT2 <- ggplot2::ggplot(kT2, aes(x = as.character(wt.group), y = log2(tpm + 1), group = wt.group)) + 
  geom_violin() 
ggsave(bpWT2, filename = file.path(dataDir, 'boxplot_sensitivity_h2az_wt_1_3_4_6.png'))

bptgfb2 <- ggplot2::ggplot(kT2, aes(x = as.character(ca1a.group), y = log2(tpm + 1), group = ca1a.group)) + 
  geom_boxplot()
vptgfb2 <- ggplot2::ggplot(kT2, aes(x = as.character(ca1a.group), y = log2(tpm + 1), group = ca1a.group)) + 
  geom_violin()
ggsave(vptgfb2, filename = file.path(dataDir, 'boxplot_sensitivity_h2az_ca1a_4_6.png'))

# enrichment analysis
m_t2g <- msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name, human_gene_symbol)

geneList <- rT[order(b,decreasing=TRUE),]$b
names(geneList) <- rT[order(b,decreasing=TRUE),]$target_id
gene <- geneList[abs(geneList) > 0.5]

em2 <- enricher(rT2[rT2$qval < 0.1 & rT2$b < 0]$target_id, pvalueCutoff = 0.01, TERM2GENE = m_t2g)
png(file = file.path(dataDir, "msigdb_dotplot_sensitivity_h2az_wt_1_3_4_6_ca1a_4_6_downregulated.png"), height = 1024, width = 768)
dotplot(em2)
dev.off()

pdf(file = file.path(dataDir, "msigdb_dotplot_sensitivity_h2az_wt_1_3_4_6_ca1a_4_6_downregulated.pdf"), height = 14, width = 20)
dotplot(em2)
dev.off()

merge(rT2[rT2$qval < 0.1 & rT2$b < 0, .(target_id, b, qval, mean_obs),][order(b, decreasing = T)], 
write.csv(rT2[rT2$qval < 0.1 & rT2$b < 0, .(target_id, b, qval, mean_obs),][order(b, decreasing = T)], 
          file = file.path(dataDir, 'downregulated_genes_sensitivity_h2az_wt_1_3_4_6_ca1a_4_6.csv'))


