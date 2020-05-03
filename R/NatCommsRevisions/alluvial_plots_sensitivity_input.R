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
dataDir <- "./alluvial_plots_sensitivity_input/"
l2 <- lapply(list.files(path = dataDir, pattern="Log2"), function(x){
  dt <- data.table::fread(paste(dataDir, x, sep = "/"))
  return(dt)
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

# For input (Figure X - sensitivity)
# ok we have duplicated gene IDs in this list
# a total of 224 which in successive merging lead to increase in total categories
geneList <- unique(l2$WT_Inp$extGene)

lapply(names(l2), function(x){
  setkey(l2[[x]], "extGene")
})
mcf10awtCategories <- unique(l2$WT_Inp[,c('group1', 'color')])
mcf10awtCategories <- mcf10awtCategories[order(group1)]
mcf10awtCategories$status <- c('active', 'inactive', 'active', 'active', 'active', 'active', 'active')
mcf10awtCategories$color2 <- mcf10awtCategories$color
mcf10awtCategories[mcf10awtCategories$status == 'inactive']$color2 <- 'grey'
mcf10awtCategories$color3 <- mcf10awtCategories$color
mcf10awtCategories[mcf10awtCategories$status == 'active']$color3 <- 'grey'
mcf10awtCategories$color4 <- c('grey', 'green', 'grey', 'grey', 'green', 'grey', 'grey')

mcf10ashzCategories <- unique(l2$shZ_k7[,c('group1', 'color')])
mcf10ashzCategories <- mcf10ashzCategories[order(group1)]
mcf10ashzCategories$status <- c('active', 'active', 'active', 'active', 'active', 'inactive', 'active')
mcf10ashzCategories$color2 <- mcf10ashzCategories$color
mcf10ashzCategories[mcf10ashzCategories$status == 'active']$color2 <- 'grey'

# merge all tables into single data.table ---------------------------------
dt1 <- l2$WT_Inp[!duplicated(extGene), c("extGene" ,"group1")]
colnames(dt1)[2] <- "wt.group"
dt1 <- merge(dt1, l2$TGFb_Inp[!duplicated(extGene), c("extGene", "group1")], by.x = "extGene", by.y = "extGene")
colnames(dt1)[length(colnames(dt1))] <- "tgfb.group"
dt1 <- merge(dt1, l2$shZ_k7[!duplicated(extGene), c("extGene", "group1")], by.x = "extGene", by.y = "extGene")
colnames(dt1)[length(colnames(dt1))] <- "shZ.group"
dt1 <- merge(dt1, l2$CA1a_Inp[!duplicated(extGene), c("extGene","group1")], by.x = "extGene", by.y = "extGene")
colnames(dt1)[length(colnames(dt1))] <- "CA1a.group"
dt1
save(dt1, file = file.path(dataDir, "alluvialPlotSensitivityInputData.rda"))


# alluvial plots WT -> TGFb ---------------------------------
fig_wt_tgfb <- data.table::as.data.table(table("WT" = dt1$wt.group, "TGFb" = dt1$tgfb.group))
fig_wt_tgfb %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> fig_wt_tgfb
pdf(file = file.path(dataDir, "alluvial_plot_sensitivity_input_wt_tgfb.pdf"))
alluvial(fig_wt_tgfb[,c(1:2)], freq = fig_wt_tgfb$n,
         col = mcf10awtCategories$color[match(as.integer(fig_wt_tgfb$WT), mcf10awtCategories$group)])
dev.off()


# alluvial plots WT -> shH2AZ -------------------------------
fig_wt_shz <- data.table::as.data.table(table("WT" = dt1$wt.group, "shH2AZ" = dt1$shZ.group))
fig_wt_shz %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> fig_wt_shz
png(file = file.path(dataDir, "alluvial_plot_sensitivity_input_wt_shz.png"))
alluvial(fig_wt_shz[,c(1:2)], freq = fig_wt_shz$n,
         col = mcf10awtCategories$color[match(as.integer(fig_wt_shz$WT), mcf10awtCategories$group)])
dev.off()


# alluvial plots WT -> CA1A ---------------------------------
fig_wt_ca1a <- data.table::as.data.table(table("WT" = dt1$wt.group, "Ca1a" = dt1$CA1a.group))
fig_wt_ca1a %>% group_by(WT, Ca1a) %>% summarise(n = sum(N)) -> fig_wt_ca1a
png(file = file.path(dataDir, "alluvial_plot_sensitivity_input_wt_ca1a.png"))
alluvial(fig_wt_ca1a[,c(1:2)], freq = fig_wt_ca1a$n,
         col = mcf10awtCategories$color[match(as.integer(fig_wt_ca1a$WT), mcf10awtCategories$group)])
dev.off()


# alluvial plots WT -> shH2AZ, active to inactive only
fig_wt_shz <- data.table::as.data.table(table("WT" = dt1$wt.group, "shH2AZ" = dt1$shZ.group))
fig_wt_shz %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> fig_wt_shz
png(file = file.path(dataDir, "alluvial_plot_sensitivity_input_wt_active_shz_inactive.png"))
alluvial(fig_wt_shz[,c(1:2)], freq = fig_wt_shz$n,
         col = mcf10awtCategories$color[match(as.integer(fig_wt_shz$WT), mcf10awtCategories$group)],
         hide = !fig_wt_shz$shH2AZ == 6)
dev.off()

# alluvial plots WT -> shH2AZ, inactive to active only
fig_wt_shz <- data.table::as.data.table(table("WT" = dt1$wt.group, "shH2AZ" = dt1$shZ.group))
fig_wt_shz %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> fig_wt_shz
png(file = file.path(dataDir, "alluvial_plot_sensitivity_input_wt_inactive_shz_active.png"))
alluvial(fig_wt_shz[,c(1:2)], freq = fig_wt_shz$n,
         col = mcf10awtCategories$color[match(as.integer(fig_wt_shz$shH2AZ), mcf10awtCategories$group)],
         hide = !fig_wt_shz$WT == 2)
dev.off()

