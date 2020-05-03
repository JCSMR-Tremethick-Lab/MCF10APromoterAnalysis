## load libraries
library(ggparallel)
library(alluvial)
library(tidyverse)
library(data.table)
# alluvial plots
# remove all variables from environment
rm(list = ls())


# load expression data ----------------------------------------------------
RDAs <- list.files("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures", 
                   pattern = "*.rda", 
                   full.names = T)
for (i in RDAs){
  load(i)
}
setkey(rT.shH2AZ, "target_id")
setkey(rT.MCF10Ca1a, "target_id")
setkey(rT.TGFbD6, "target_id")
setkey(kT1, "target_id")

# load parallel plot data -------------------------------------------------
dataDir <- "."
l2 <- lapply(list.files(path = dataDir, pattern="Log2"), function(x){
  dt1 <- data.table::fread(paste(dataDir, x, sep = "/"))
  return(dt1)
})
n2 <- gsub(x = list.files(path = dataDir, pattern="Log2"), pattern = ".tsv", replacement =  "") 
names(l2) <- unlist(lapply(strsplit(n2, "_"), function(x) paste(x[2:3], collapse = "_")))

l2 <- lapply(l2, function(x) {
  ucscID <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[4]))
  extGene <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[6]))
  x$ucscID <- ucscID
  x$extGene <- extGene
  x$chr <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[1]))
  x$start <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[2]))
  x$end <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[3]))
  return(x)
})
l2

# For h2az (Figure X - sensitivity)
# ok we have duplicated gene IDs in this list
# a total of 224 which in successive merging lead to increase in total categories
geneList <- unique(l2$WT_H2AZ$extGene)

lapply(names(l2), function(x){
  setkey(l2[[x]], "extGene")
})
mcf10awtCategories <- unique(l2$WT_H2AZ[,c('group1', 'color')])

# merge all tables into single data.table ---------------------------------
dt1 <- l2$WT_H2AZ[!duplicated(extGene), c("extGene" ,"group1")]
colnames(dt1)[2] <- "wt.group"
dt1 <- merge(dt1, l2$TFGb_H2AZ[!duplicated(extGene), c("extGene", "group1")], by.x = "extGene", by.y = "extGene")
colnames(dt1)[length(colnames(dt1))] <- "tgfb.group"
dt1 <- merge(dt1, l2$CA1a_H2AZ[!duplicated(extGene), c("extGene","group1")], by.x = "extGene", by.y = "extGene")
colnames(dt1)[length(colnames(dt1))] <- "ca1a.group"
dt1
save(dt1, file = file.path(dataDir, "alluvialPlotsSensitivityH2AZData.rda"))


# alluvial plots WT -> TGFb ---------------------------------
fig_wt_tgfb <- data.table::as.data.table(table("WT" = dt1$wt.group, "TGFb" = dt1$tgfb.group))
fig_wt_tgfb %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> fig_wt_tgfb
png(file = file.path(dataDir, "alluvial_plot_sensitivity_h2az_wt_tgfb.png"))
alluvial(fig_wt_tgfb[,c(1:2)], freq = fig_wt_tgfb$n,
         col = mcf10awtCategories$color[match(as.integer(fig_wt_tgfb$WT), mcf10awtCategories$group)])
dev.off()


# alluvial plots WT -> CA1A ---------------------------------
fig_wt_ca1a <- data.table::as.data.table(table("WT" = dt1$wt.group, "Ca1a" = dt1$ca1a.group))
fig_wt_ca1a %>% group_by(WT, Ca1a) %>% summarise(n = sum(N)) -> fig_wt_ca1a
png(file = file.path(dataDir, "alluvial_plot_sensitivity_h2az_wt_ca1a.png"), )
alluvial(fig_wt_ca1a[,c(1:2)], freq = fig_wt_ca1a$n,
         col = mcf10awtCategories$color[match(as.integer(fig_wt_ca1a$WT), mcf10awtCategories$group)])
dev.off()

