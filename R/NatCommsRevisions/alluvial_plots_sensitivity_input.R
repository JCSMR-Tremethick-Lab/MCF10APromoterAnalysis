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


# alluvial plots
# remove all variables from environment
rm(list = ls())

dataDir <- "./alluvial_plots_sensitivity_input/"

# load expression data ----------------------------------------------------
load(file = file.path(dataDir, 'so_wt_shz_d8.rda'))
rT <- sleuth::sleuth_results(so, test = 'conditionshZ', test_type = 'wt', show_all = F)
setDT(rT)
setkey(rT, target_id)  

kT <- sleuth::kallisto_table(so, normalized = T, include_covariates = T)                  
setDT(kT)
kTMean <- kT[, lapply(.SD, mean), by = list(condition, target_id), .SDcols='tpm']
setkey(kTMean, target_id)
kTMeanWide <- dcast(kTMean, target_id ~ condition, value.var = 'tpm')

# load parallel plot data -------------------------------------------------
l1 <- lapply(list.files(path = dataDir, pattern="Log2"), function(x){
  dt <- data.table::fread(paste(dataDir, x, sep = "/"))
  return(dt)
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

# For input (Figure X - sensitivity)
# ok we have duplicated gene IDs in this list
# a total of 224 which in successive merging lead to increase in total categories
geneList <- unique(l1$WT_Inp$extGene)

lapply(names(l1), function(x){
  setkey(l1[[x]], "extGene")
})
mcf10awtCategories <- unique(l1$WT_Inp[,c('group1', 'color')])
mcf10awtCategories <- mcf10awtCategories[order(group1)]

# merge all tables into single data.table ---------------------------------
dt1 <- l1$WT_Inp[!duplicated(extGene), c("extGene" ,"group1")]
colnames(dt1)[2] <- "wt.group"
dt1 <- merge(dt1, l1$TGFb_Inp[!duplicated(extGene), c("extGene", "group1")], by.x = "extGene", by.y = "extGene")
colnames(dt1)[length(colnames(dt1))] <- "tgfb.group"
dt1 <- merge(dt1, l1$shZ_k7[!duplicated(extGene), c("extGene", "group1")], by.x = "extGene", by.y = "extGene")
colnames(dt1)[length(colnames(dt1))] <- "shZ.group"
dt1 <- merge(dt1, l1$CA1a_Inp[!duplicated(extGene), c("extGene","group1")], by.x = "extGene", by.y = "extGene")
colnames(dt1)[length(colnames(dt1))] <- "CA1a.group"
dt1
save(dt1, file = file.path(dataDir, "alluvialPlotSensitivityInputData.rda"))

# alluvial plots WT -> shH2AZ -------------------------------
tab <- table("WT" = dt1$wt.group, "shH2AZ" = dt1$shZ.group)
rownames(tab) <- paste('WT_cluster_', c(1:7), sep = '')
colnames(tab) <- paste('shH2AZ_cluster_', c(1:7), sep = '')
write.csv(tab, file = file.path(dataDir, "cross_table_WT_shH2AZ.csv"))

fig_wt_shz <- data.table::as.data.table(table("WT" = dt1$wt.group, "shH2AZ" = dt1$shZ.group))
fig_wt_shz %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> fig_wt_shz
png(file = file.path(dataDir, "alluvial_plot_sensitivity_input_wt_shz.png"))
alluvial::alluvial(fig_wt_shz[,c(1:2)], freq = fig_wt_shz$n,
         col = mcf10awtCategories$color[match(as.integer(fig_wt_shz$WT), mcf10awtCategories$group)])
dev.off()

# alluvial plots WT -> shH2AZ --------------------
# 2,5 to 1,2,3,4,7 - generally: inactive to active
#fig_wt_shz <- data.table::as.data.table(table("WT" = dt1$wt.group, "shH2AZ" = dt1$shZ.group))
#fig_wt_shz %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> fig_wt_shz

mcf10awtCategories$color1 <- mcf10awtCategories$color
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(2,5)]$color1 <- 'grey'
mcf10awtCategories$alpha1 <- 1
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(2,5)]$alpha1 <- 0.1

png(file = file.path(dataDir, "alluvial_plot_sensitivity_input_wt_2_5_shz_1_2_3_4_7.png"))
alluvial(fig_wt_shz[,c(1:2)], freq = fig_wt_shz$n,
         col = mcf10awtCategories$color1[match(as.integer(fig_wt_shz$WT), mcf10awtCategories$group1)],
         alpha = mcf10awtCategories$alpha1[match(as.integer(fig_wt_shz$WT), mcf10awtCategories$group1)],
         hide = fig_wt_shz$shH2AZ %in% c(5,6))
dev.off()

# add table of genes going from 2 & 5 WT to 1,2,3,4,7 shZ
selectedGenes1 <- dt1[dt1$shZ.group %in% c(1,2,3,4,7) & (dt1$wt.group %in% c(2,5))]$extGene # list of genes
geneTable1 <- AnnotationDbi::select(org.Hs.eg.db, keys = selectedGenes1, keytype = 'SYMBOL', columns = c('SYMBOL', 'GENENAME'))
setDT(geneTable1)
geneTable1 <- merge(dt1[,.(extGene, wt.group, shZ.group)], geneTable1, by.x = 'extGene', by.y = 'SYMBOL')
geneTable1 <- merge(geneTable1, kTMeanWide, by.x = 'extGene', by.y = 'target_id', all.x = T, all.y = F)
geneTable1 <- merge(geneTable1, rT[,.(target_id, b, qval)], by.x = 'extGene', by.y = 'target_id', all.x = T, all.y = F)
write.csv(geneTable1, file = file.path(dataDir, 'gene_table_sensitivity_input_wt_2_5_shz_1_2_3_4_7.csv'))

# differential expression of these
rT1 <- rT[selectedGenes1]
rT1 <- rT1[!is.na(pval)]
write.csv(rT1, file = file.path(dataDir, 'diffGenes_WT_2_5_shH2AZ_1_2_3_4_7.csv'))
hist(rT1$qval)
plot(rT1$b, -log10(rT1$qval))
setkey(kT, target_id)
kT1 <- kT[selectedGenes1]
kT1 <- kT1[!is.na(kT1$tpm)]
kT1 <- merge(kT1[kT1$condition %in% c('WT', 'shZ')], dt1[,c('extGene','wt.group','shZ.group')], by.x = 'target_id', by.y = 'extGene', all.x = T, all.y = F)
table(kT1$shZ.group)
table(kT1$wt.group)

bpWT1 <- ggplot2::ggplot(kT1, aes(x = as.character(wt.group), y = log2(tpm + 1), group = wt.group)) + 
  geom_boxplot() 
vpWT1 <- ggplot2::ggplot(kT1, aes(x = as.character(wt.group), y = log2(tpm + 1), group = wt.group)) + 
  geom_violin() 
ggsave(bpWT1, filename = file.path(dataDir, 'boxplot_sensitivity_input_wt_2_5.png'))

bpshZ1 <- ggplot2::ggplot(kT1, aes(x = as.character(shZ.group), y = log2(tpm + 1), group = shZ.group)) + 
  geom_boxplot()
vpshZ1 <- ggplot2::ggplot(kT1, aes(x = as.character(shZ.group), y = log2(tpm + 1), group = shZ.group)) + 
  geom_violin()
ggsave(bpshZ1, filename = file.path(dataDir, 'boxplot_sensitivity_input_shz_1_2_3_4_7.png'))

#####TODO include enrichment analysis
m_t2g <- msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name, human_gene_symbol)

geneList <- rT[order(b,decreasing=TRUE),]$b
names(geneList) <- rT[order(b,decreasing=TRUE),]$target_id
gene <- geneList[abs(geneList) > 0.5]

em1 <- enricher(rT1[rT1$qval < 0.1 & rT1$b > 0]$target_id, pvalueCutoff = 0.01, TERM2GENE = m_t2g)
png(file = file.path(dataDir, "msigdb_dotplot_sensitivity_input_wt_2_5_shz_1_2_3_4_7_upregulated.png"), height = 1024, width = 768)
dotplot(em1)
dev.off()

write.csv(rT1[rT1$qval < 0.1 & rT1$b > 0, .(target_id, b, qval, mean_obs),][order(b, decreasing = T)], file = file.path(dataDir, 'upregulated_genes_sensitivit_input_wt_2_5_shz_1_2_3_4_7.csv'))

# alluvial plots WT -> shH2AZ --------
# 1,4,6,7  to 6
# active to inactive
mcf10awtCategories$color2 <- mcf10awtCategories$color
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(1,4,6,7)]$color2 <- 'grey'
mcf10awtCategories$alpha2 <- 1
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(1,4,6,7)]$alpha2 <- 0.1

png(file = file.path(dataDir, "alluvial_plot_sensitivity_input_wt_1_4_6_7_shz_6.png"))
alluvial(fig_wt_shz[,c(1:2)], freq = fig_wt_shz$n,
         col = mcf10awtCategories$color2[match(as.integer(fig_wt_shz$WT), mcf10awtCategories$group)],
         alpha = mcf10awtCategories$alpha2[match(as.integer(fig_wt_shz$WT), mcf10awtCategories$group)],
         hide = !fig_wt_shz$shH2AZ %in% c(6))
dev.off()

# add table of genes going from 1,4,6,7 WT to 6 shZ
selectedGenes2 <- dt1[dt1$shZ.group %in% c(6) & (dt1$wt.group %in% c(1,4,6,7))]$extGene # list of genes
geneTable2 <- AnnotationDbi::select(org.Hs.eg.db, keys = selectedGenes2, keytype = 'SYMBOL', columns = c('SYMBOL', 'GENENAME'))
setDT(geneTable2)
geneTable2 <- merge(dt1[,.(extGene, wt.group, shZ.group)], geneTable2, by.x = 'extGene', by.y = 'SYMBOL')
geneTable2 <- merge(geneTable2, kTMeanWide, by.x = 'extGene', by.y = 'target_id', all.x = T, all.y = F)
geneTable2 <- merge(geneTable2, rT[,.(target_id, b, qval)], by.x = 'extGene', by.y = 'target_id', all.x = T, all.y = F)
write.csv(geneTable2, file = file.path(dataDir, 'gene_table_sensitivity_input_wt_1_4_6_7_shz_6.csv'))

# differential expression of these
rT2 <- rT[selectedGenes2]
rT2 <- rT2[!is.na(pval)]
write.csv(rT2, file = file.path(dataDir, 'diffGenes_WT_1_4_6_7_shH2AZ_6.csv'))

plot(rT2$b, -log10(rT2$qval))

kT2 <- kT[selectedGenes2]
kT2 <- kT2[!is.na(kT2$tpm)]
kT2 <- merge(kT2[kT2$condition %in% c('WT', 'shZ')], dt1[,c('extGene','wt.group','shZ.group')], by.x = 'target_id', by.y = 'extGene', all.x = T, all.y = F)
table(kT2$shZ.group)
table(kT2$wt.group)

bpWT2 <- ggplot2::ggplot(kT2, aes(x = as.character(wt.group), y = log2(tpm + 1), group = wt.group)) + 
  geom_boxplot() 
ggsave(bpWT2, filename = file.path(dataDir, 'boxplot_sensitivity_input_wt_1_4_6_7.png'))

bpshZ2 <- ggplot2::ggplot(kT2, aes(x = as.character(shZ.group), y = log2(tpm + 1), group = shZ.group)) + 
  geom_boxplot()
ggsave(bpshZ2, filename = file.path(dataDir, 'boxplot_sensitivity_input_shz_6.png'))




#####TODO include enrichment analysis
em2 <- enricher(rT2[rT2$qval < 0.1 & rT2$b < 0]$target_id, pvalueCutoff = 0.01, TERM2GENE = m_t2g)

png(file = file.path(dataDir, "msigdb_dotplot_sensitivity_input_wt_1_4_6_7_shz_6_downregulated.png"), height = 1024, width = 768)
dotplot(em2)
dev.off()

write.csv(rT2[rT2$qval < 0.1 & rT2$b < 0, .(target_id, b, qval, mean_obs),][order(b, decreasing = T)], file = file.path(dataDir, 'downregulated_genes_sensitivity_input_wt_1_4_6_7_shz_6_downregulated.csv'))

