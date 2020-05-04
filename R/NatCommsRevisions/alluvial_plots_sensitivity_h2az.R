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
dataDir <- "./alluvial_plots_sensitivity_h2az/"
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
setorder(mcf10awtCategories, group1)

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
# 5,7 to 1,2,6
mcf10awtCategories$color2 <- mcf10awtCategories$color
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(5,7)]$color2 <- "grey"
mcf10awtCategories$alpha2 <- 1
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(5,7)]$alpha2 <- 0.1
fig_wt_tgfb <- data.table::as.data.table(table("WT" = dt1$wt.group, "TGFb" = dt1$tgfb.group))
fig_wt_tgfb %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> fig_wt_tgfb
png(file = file.path(dataDir, "alluvial_plot_sensitivity_h2az_wt_5_7_tgfb_1_2_6.png"))
alluvial(fig_wt_tgfb[,c(1:2)], freq = fig_wt_tgfb$n,
         col = mcf10awtCategories$color2[match(as.integer(fig_wt_tgfb$WT), mcf10awtCategories$group)],
         hide = !fig_wt_tgfb$TGFb %in% c(1,2,6),
         alpha = mcf10awtCategories$alpha2[match(as.integer(fig_wt_tgfb$WT), mcf10awtCategories$group)])
dev.off()

# 1,3,4,6    to     5,7 (just chose 2)
mcf10awtCategories$color3 <- mcf10awtCategories$color
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(1,3,4,6)]$color3 <- "grey"
mcf10awtCategories$alpha3 <- 1
mcf10awtCategories[!mcf10awtCategories$group1 %in% c(1,3,4,6)]$alpha3 <- 0.1
png(file = file.path(dataDir, "alluvial_plot_sensitivity_h2az_wt_1_3_4_6_tgfb_5_7.png"))
alluvial(fig_wt_tgfb[,c(1:2)], freq = fig_wt_tgfb$n,
         col = mcf10awtCategories$color3[match(as.integer(fig_wt_tgfb$WT), mcf10awtCategories$group)],
         hide = !fig_wt_tgfb$TGFb %in% c(5,7),
         alpha = mcf10awtCategories$alpha3[match(as.integer(fig_wt_tgfb$WT), mcf10awtCategories$group)])
dev.off()


# alluvial plots WT -> CA1A ---------------------------------
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
