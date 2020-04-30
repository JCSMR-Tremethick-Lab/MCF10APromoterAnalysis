## load libraries
library(ggparallel)
library(alluvial)
library(tidyverse)
library(data.table)
# Figure S2 (supplement for Figure 2)
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
dataDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/Figure_2"
l2 <- lapply(list.files(path = dataDir, pattern=".tsv"), function(x){
  dt <- data.table::fread(paste(dataDir, x, sep = "/"))
  return(dt)
})
n2 <- gsub(x = list.files(path = dataDir, pattern=".tsv"), pattern = ".tsv", replacement =  "") 
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

# For input (Figure 2 - sensitivity)
# 1. Strong -1 sensitive (1 for 10A, 3 for TGFb, 3 for shH2AZ, 3 for CAC1).
# 2. neither sensitive or resistant across promoter (2 for 10A, 4 for TGFb, 6 for shH2AZ, 5 for CAC1).
# 3. Strong +1 (3 for 10A, 7 for TGFb, 2 & 1 for shH2AZ, 7 for CAC1)
# 4. resistant upstream (4 for 10A only seen here)
# 5. strong -2(5 for 10A, 2 for TGFb, 7 for shH2AZ, 6 for CAC1)
# 6. TSS resistant, sensitive upstream (6 for 10A, 1 for TGFb, 4 for CAC1)-not seen in shH2A.Z
# 7. resistant entire promoter (7 for 10A, 6 for TGFb, 2 for CAC1)-not seen in shH2AZ
# 8. sensitive -3 (5, only TGFb,)
# 9. sensitive TSS (5, only seen is shH2A.Z).

# WT
mcf10awtCategories <- data.table::data.table(cat = c("Strong -1 sensitive", 
                                                     "Neither sensitive or resistant across promoter", 
                                                     "Strong +1", 
                                                     "Resistant upstream", 
                                                     "Strong -2", 
                                                     "TSS resistant, sensitive upstream",
                                                     "Resistant entire promoter"),
                                             group = c(1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6,
                                                       7))
mcf10awtCategories$cat <- factor(mcf10awtCategories$cat, levels = mcf10awtCategories$cat, ordered = T)

#TGFb
mcf10atgfbCategories <- data.table::data.table(cat = c(
                                                       "TSS resistant, sensitive upstream", #1
                                                       "Strong -2", #2
                                                       "Strong -1 sensitive", #3
                                                       "Neither sensitive or resistant across promoter", #4
                                                       "Sensitive -3 ", #5
                                                       "Resistant entire promoter", #6
                                                       "Strong +1" #7
                                                       ), 
                                               group = c(1,
                                                         2,
                                                         3,
                                                         4,
                                                         5,
                                                         6,
                                                         7))
mcf10atgfbCategories$cat <- factor(mcf10atgfbCategories$cat, levels = mcf10atgfbCategories$cat, ordered = T)
mcf10awtCategories <- merge(mcf10awtCategories, unique(l2$log2ratio_10A[, c("group1", "color")]), by.x = "group", by.y = "group1", all.x = F, all.y = F)

# shH2AZ
# can't have duplicate categories so separated Strong +1 in to a & b
mcf10ashh2azCategories <- data.table::data.table(cat = c(
                                                         "Strong +1a", #1
                                                         "Strong +1b", #2
                                                         "Strong -1 sensitive", #3
                                                         "Sensitive +2", #4
                                                         "Sensitive TSS",#5
                                                         "Neither sensitive or resistant across promoter", #6
                                                         "Strong -2" #7
                                                         ),
                                                 group = c(1,
                                                           2,
                                                           3,
                                                           4,
                                                           5,
                                                           6,
                                                           7))
mcf10ashh2azCategories$cat <- factor(mcf10ashh2azCategories$cat, levels = mcf10ashh2azCategories$cat, ordered = T)
mcf10ashh2azCategories <- merge(mcf10ashh2azCategories, unique(l2$log2ratio_shZ[, c("group1", "color")]), by.x = "group", by.y = "group1", all.x = F, all.y = F)

# Ca1a
mcf10ca1aCategories <- data.table::data.table(cat = c(
                                                      "Sensitive + 1", #1
                                                      "Resistant entire promoter", #2
                                                      "Strong -1 sensitive", #3
                                                      "TSS resistant, sensitive upstream", #4
                                                      "Neither sensitive or resistant across promoter", #5
                                                      "Strong -2", #6
                                                      "Strong +1" #7
                                                      ),
                                              group = c(1,
                                                        2,
                                                        3,
                                                        4,
                                                        5,
                                                        6,
                                                        7))

mcf10ca1aCategories$cat <- factor(mcf10ca1aCategories$cat, levels = mcf10ca1aCategories$cat, ordered = T)
mcf10ca1aCategories <- merge(mcf10ca1aCategories, unique(l2$log2ratio_CA1a[, c("group1", "color")]), by.x = "group", by.y = "group1", all.x = F, all.y = F)


# merge clusters and categories -------------------------------------------
l2$log2ratio_10A$left.MCF10A_WT.category <- as.factor(mcf10awtCategories$cat[match(l2$log2ratio_10A$group1, mcf10awtCategories$group)])
l2$log2ratio_TGFb$right.MCF10A_TGFb.category <- as.factor(mcf10atgfbCategories$cat[match(l2$log2ratio_TGFb$group1, mcf10atgfbCategories$group)])
l2$log2ratio_shZ$right.MCF10A_shZ.category <- as.factor(mcf10ashh2azCategories$cat[match(l2$log2ratio_shZ$group1, mcf10ashh2azCategories$group)])
l2$log2ratio_CA1a$right.MCF10CA1A.category <- as.factor(mcf10ca1aCategories$cat[match(l2$log2ratio_CA1a$group1, mcf10ca1aCategories$group)])

# ok we have duplicated gene IDs in this list
# a total of 224 which in successive merging lead to increase in total categories
geneList <- unique(l2$log2ratio_10A$extGene)

lapply(names(l2), function(x){
  setkey(l2[[x]], "extGene")
})


# merge all tables into single data.table ---------------------------------
dt2 <- l2$log2ratio_10A[!duplicated(extGene), c("extGene", "left.MCF10A_WT.category" ,"group1")]
colnames(dt2)[3] <- "wt.group"
dt2 <- merge(dt2, l2$log2ratio_TGFb[!duplicated(extGene), c("extGene", "right.MCF10A_TGFb.category", "group1")], by.x = "extGene", by.y = "extGene")
colnames(dt2)[length(colnames(dt2))] <- "tgfb.group"
dt2 <- merge(dt2, l2$log2ratio_shZ[!duplicated(extGene), c("extGene", "right.MCF10A_shZ.category", "group1")], by.x = "extGene", by.y = "extGene")
colnames(dt2)[length(colnames(dt2))] <- "shZ.group"
dt2 <- merge(dt2, l2$log2ratio_CA1a[!duplicated(extGene), c("extGene", "right.MCF10CA1A.category", "group1")], by.x = "extGene", by.y = "extGene")
colnames(dt2)[length(colnames(dt2))] <- "CA1a.group"
dt2
FigS2Data <- dt2
save(FigS2Data, file = file.path(dataDir, "FigS2Data.rda"))
write.csv(FigS2Data, file = file.path(dataDir, "FigS2Data.csv"))


# alluvial plots - Figures S2a WT -> TGFb ---------------------------------
tab1 <- data.table::as.data.table(table(FigS2Data$wt.group, FigS2Data$tgfb.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2aTab <- data.table::as.data.table(table("WT" = FigS2Data$wt.group, "TGFb" = FigS2Data$tgfb.group))
FigS2aTab %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS2aTab1
pdf(file = file.path(dataDir, "AlluvialPlotS2a.pdf"), paper = "a4r")
alluvial(FigS2aTab1[,c(1:2)], freq = FigS2aTab1$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()


# alluvial plots - Figures S2b WT -> shH2AZ -------------------------------
tab1 <- data.table::as.data.table(table(FigS2Data$wt.group, FigS2Data$shZ.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2bTab <- data.table::as.data.table(table("WT" = FigS2Data$wt.group, "shH2AZ" = FigS2Data$shZ.group))
FigS2bTab %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> FigS2bTab1
pdf(file = file.path(dataDir, "AlluvialPlotS2b.pdf"), paper = "a4r")
alluvial(FigS2bTab1[,c(1:2)], freq = FigS2bTab1$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()


# alluvial plots - Figures S2a WT -> CA1A ---------------------------------
tab1 <- data.table::as.data.table(table(FigS2Data$wt.group, FigS2Data$CA1a.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2cTab <- data.table::as.data.table(table("WT" = FigS2Data$wt.group, "CA1a" = FigS2Data$CA1a.group))
FigS2cTab %>% group_by(WT, CA1a) %>% summarise(n = sum(N)) -> FigS2cTab1
pdf(file = file.path(dataDir, "AlluvialPlotS2c.pdf"), paper = "a4r")
alluvial(FigS2cTab1[,c(1:2)], freq = FigS2cTab1$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

#
apply(dt2[,c("left.MCF10A_WT.category", "right.MCF10A_TGFb.category", "right.MCF10A_shZ.category", "right.MCF10CA1A.category")], 2, function(x) table(x, useNA = "always"))

write.csv(table("WT" = dt2$wt.group, "TGFb" = dt2$tgfb.group), 
          file = file.path(dataDir, "Fig2_WT_TGFb_counts.csv"),
          row.names = T, col.names = T)
write.csv(table("WT" = dt2$wt.group, "shZ" = dt2$shZ.group),
          file = file.path(dataDir, "Fig2_WT_shZ_counts.csv"))
write.csv(table("WT" = dt2$wt.group, "CA1A" = dt2$CA1a.group),
          file = file.path(dataDir, "Fig2_WT_CA1A_counts.csv"))

# re-run analysis with EMT genes only - data preparation -----------
sigEMTCells <- data.table::fread("/Data/References/Annotations/GeneSets/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_cellLine.txt")
sigEMTTumor <- data.table::fread("/Data/References/Annotations/GeneSets/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_tumor.txt")
colnames(sigEMTCells)[1] <- "gene_symbol"
colnames(sigEMTTumor)[1] <- "gene_symbol"
setkey(sigEMTTumor, "gene_symbol")
sigEMT <- rbind(sigEMTCells,
                sigEMTTumor[!intersect(sigEMTCells$gene_symbol, sigEMTTumor$gene_symbol)])
setkey(dt2, "extGene")
emtData <- dt2[sigEMT$gene_symbol][!is.na(wt.group)]
emtData <- merge(emtData, 
                 sigEMT, 
                 by.x = "extGene", 
                 by.y = "gene_symbol", 
                 all.x = T, 
                 all.y = F)

# FigS2a WT -> TGFb all EMT genes ----------------------------------
tab1 <- data.table::as.data.table(table(emtData$wt.group, emtData$tgfb.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2aTabEMTAll <- data.table::as.data.table(table("WT" = emtData$wt.group, "TGFb" = emtData$tgfb.group))
FigS2aTabEMTAll %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS2aTabEMTAll1
pdf(file = file.path(dataDir, "FigS2aTabEMTAll.pdf"), paper = "a4r")
alluvial(FigS2aTabEMTAll1[,c(1:2)], freq = FigS2aTabEMTAll1$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()


# FigS2a WT -> TGFb epithelial genes --------------------------------------
tab1 <- data.table::as.data.table(table(emtData[epi_mes == "epi"]$wt.group, 
                                        emtData[epi_mes == "epi"]$tgfb.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2aEMTgenesEpi <- data.table::as.data.table(table("WT" = emtData[epi_mes == "epi"]$wt.group, 
                                                     "TGFb" = emtData[epi_mes == "epi"]$tgfb.group))
FigS2aEMTgenesEpi %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS2aEMTgenesEpi
pdf(file = file.path(dataDir, "FigS2aEMTgenesEpi.pdf"), paper = "a4r")
alluvial(FigS2aEMTgenesEpi[,c(1:2)], freq = FigS2aEMTgenesEpi$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2a WT -> TGFb down-regulated epithelial genes -----------------------
tab1 <- data.table::as.data.table(table(rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0]$wt.group, 
                                        rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0]$tgfb.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2aEMTgenesEpiDown <- data.table::as.data.table(table("WT" = rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0]$wt.group, 
                                                         "TGFb" = rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0]$tgfb.group))
FigS2aEMTgenesEpiDown %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS2aEMTgenesEpiDown
pdf(file = file.path(dataDir, "FigS2aEMTgenesEpiDown.pdf"), paper = "a4r")
alluvial(FigS2aEMTgenesEpiDown[,c(1:2)], freq = FigS2aEMTgenesEpiDown$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2a WT -> TGFb all mesenchymal genes ---------------------------------
tab1 <- data.table::as.data.table(table(emtData[epi_mes == "mes"]$wt.group, 
                                        emtData[epi_mes == "mes"]$tgfb.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2aEMTgenesMes <- data.table::as.data.table(table("WT" = emtData[epi_mes == "mes"]$wt.group, 
                                                         "TGFb" = emtData[epi_mes == "mes"]$tgfb.group))
FigS2aEMTgenesMes %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS2aEMTgenesMes
pdf(file = file.path(dataDir, "FigS2aEMTgenesMes.pdf"), paper = "a4r")
alluvial(FigS2aEMTgenesMes[,c(1:2)], freq = FigS2aEMTgenesMes$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2a WT -> TGFb up-regulated mesenchymal genes ---------------------------------
tab1 <- data.table::as.data.table(table(rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0]$wt.group, 
                                        rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0]$tgfb.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2aEMTgenesMesUp <- data.table::as.data.table(table("WT" = rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0]$wt.group, 
                                                     "TGFb" = rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0]$tgfb.group))
FigS2aEMTgenesMesUp %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS2aEMTgenesMesUp
pdf(file = file.path(dataDir, "FigS2aEMTgenesMesUp.pdf"), paper = "a4r")
alluvial(FigS2aEMTgenesMesUp[,c(1:2)], freq = FigS2aEMTgenesMesUp$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# CONTINUE HERE  - 2018-11-19
# FigS2b WT -> shH2AZ - all EMT genes -------------------------------------

tab1 <- data.table::as.data.table(table(emtData$wt.group, emtData$shZ.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2bTabEMTAll <- data.table::as.data.table(table("WT" = emtData$wt.group, "shH2AZ" = emtData$shZ.group))
FigS2bTabEMTAll %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> FigS2bTabEMTAll
pdf(file = file.path(dataDir, "FigS2bTabEMTAll.pdf"), paper = "a4r")
alluvial(FigS2bTabEMTAll[,c(1:2)], freq = FigS2bTabEMTAll$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2b WT -> shH2AZ - all epithelial genes -------------------------------------
tab1 <- data.table::as.data.table(table(emtData[epi_mes == "epi"]$wt.group, 
                                        emtData[epi_mes == "epi"]$shZ.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2bEMTgenesEpi <- data.table::as.data.table(table("WT" = emtData[epi_mes == "epi"]$wt.group, 
                                                     "shH2AZ" = emtData[epi_mes == "epi"]$shZ.group))
FigS2bEMTgenesEpi %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> FigS2bEMTgenesEpi
pdf(file = file.path(dataDir, "FigS2bEMTgenesEpi.pdf"), paper = "a4r")
alluvial(FigS2bEMTgenesEpi[,c(1:2)], freq = FigS2bEMTgenesEpi$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2b WT -> shH2AZ - down-regulated epithelial genes -------------------------------------
tab1 <- data.table::as.data.table(table(rT.shH2AZ[emtData[epi_mes == "epi"]][b < 0]$wt.group, 
                                        rT.shH2AZ[emtData[epi_mes == "epi"]][b < 0]$shZ.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2bEMTgenesEpiDown <- data.table::as.data.table(table("WT" = rT.shH2AZ[emtData[epi_mes == "epi"]][b < 0]$wt.group, 
                                                         "shH2AZ" = rT.shH2AZ[emtData[epi_mes == "epi"]][b < 0]$shZ.group))
FigS2bEMTgenesEpiDown %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> FigS2bEMTgenesEpiDown
pdf(file = file.path(dataDir, "FigS2bEMTgenesEpiDown.pdf"), paper = "a4r")
alluvial(FigS2bEMTgenesEpiDown[,c(1:2)], freq = FigS2bEMTgenesEpiDown$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2b WT -> shH2AZ - mesenchymal genes -------------------------------------
tab1 <- data.table::as.data.table(table(emtData[epi_mes == "mes"]$wt.group, 
                                        emtData[epi_mes == "mes"]$shZ.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2bEMTgenesMes <- data.table::as.data.table(table("WT" = emtData[epi_mes == "mes"]$wt.group, 
                                                     "shH2AZ" = emtData[epi_mes == "mes"]$shZ.group))
FigS2bEMTgenesMes %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> FigS2bEMTgenesMes
pdf(file = file.path(dataDir, "FigS2bEMTgenesMes.pdf"), paper = "a4r")
alluvial(FigS2bEMTgenesMes[,c(1:2)], freq = FigS2bEMTgenesMes$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2b WT -> shH2AZ - up-regulated mesenchymal genes -------------------------------------
tab1 <- data.table::as.data.table(table(rT.shH2AZ[emtData[epi_mes == "mes"]][b > 0]$wt.group, 
                                        rT.shH2AZ[emtData[epi_mes == "mes"]][b > 0]$shZ.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2bEMTgenesMesUp <- data.table::as.data.table(table("WT" = rT.shH2AZ[emtData[epi_mes == "mes"]][b > 0]$wt.group, 
                                                       "shH2AZ" = rT.shH2AZ[emtData[epi_mes == "mes"]][b > 0]$tgfb.group))
FigS2bEMTgenesMesUp %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> FigS2bEMTgenesMesUp
pdf(file = file.path(dataDir, "FigS2bEMTgenesMesUp.pdf"), paper = "a4r")
alluvial(FigS2bEMTgenesMesUp[,c(1:2)], freq = FigS2bEMTgenesMesUp$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2c WT -> Ca1a all EMT genes -----------------------------------------
tab1 <- data.table::as.data.table(table(emtData$wt.group, emtData$CA1a.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2cTabEMTAll <- data.table::as.data.table(table("WT" = emtData$wt.group, "CA1A" = emtData$CA1a.group))
FigS2cTabEMTAll %>% group_by(WT, CA1A) %>% summarise(n = sum(N)) -> FigS2cTabEMTAll
pdf(file = file.path(dataDir, "FigS2cTabEMTAll.pdf"), paper = "a4r")
alluvial(FigS2cTabEMTAll[,c(1:2)], freq = FigS2cTabEMTAll$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2c WT -> Ca1a epithelial EMT genes -----------------------------------------
tab1 <- data.table::as.data.table(table(emtData[epi_mes == "epi"]$wt.group, 
                                        emtData[epi_mes == "epi"]$CA1a.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2cEMTgenesEpi <- data.table::as.data.table(table("WT" = emtData[epi_mes == "epi"]$wt.group, 
                                                     "CA1A" = emtData[epi_mes == "epi"]$CA1a.group))
FigS2cEMTgenesEpi %>% group_by(WT, CA1A) %>% summarise(n = sum(N)) -> FigS2bEMTgenesEpi
pdf(file = file.path(dataDir, "FigS2bEMTgenesEpi.pdf"), paper = "a4r")
alluvial(FigS2bEMTgenesEpi[,c(1:2)], freq = FigS2bEMTgenesEpi$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2c WT -> Ca1a - down-regulated epithelial genes -------------------------------------
tab1 <- data.table::as.data.table(table(rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0]$wt.group, 
                                        rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0]$CA1a.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2cEMTgenesEpiDown <- data.table::as.data.table(table("WT" = rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0]$wt.group, 
                                                         "CA1A" = rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0]$CA1a.group))
FigS2cEMTgenesEpiDown %>% group_by(WT, CA1A) %>% summarise(n = sum(N)) -> FigS2cEMTgenesEpiDown
pdf(file = file.path(dataDir, "FigS2cEMTgenesEpiDown.pdf"), paper = "a4r")
alluvial(FigS2cEMTgenesEpiDown[,c(1:2)], freq = FigS2cEMTgenesEpiDown$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2c WT -> Ca1a mesenchymal EMT genes -----------------------------------------
tab1 <- data.table::as.data.table(table(emtData[epi_mes == "mes"]$wt.group, 
                                        emtData[epi_mes == "mes"]$CA1a.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2cEMTgenesMes <- data.table::as.data.table(table("WT" = emtData[epi_mes == "mes"]$wt.group, 
                                                     "CA1A" = emtData[epi_mes == "mes"]$CA1a.group))
FigS2cEMTgenesMes %>% group_by(WT, CA1A) %>% summarise(n = sum(N)) -> FigS2cEMTgenesMes
pdf(file = file.path(dataDir, "FigS2cEMTgenesMes.pdf"), paper = "a4r")
alluvial(FigS2cEMTgenesMes[,c(1:2)], freq = FigS2cEMTgenesMes$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS2c WT -> Ca1a - up-regulated mesenchymal genes -------------------------------------
tab1 <- data.table::as.data.table(table(rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0]$wt.group, 
                                        rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0]$CA1a.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2cEMTgenesMesUp <- data.table::as.data.table(table("WT" = rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0]$wt.group, 
                                                       "CA1A" = rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0]$CA1a.group))
FigS2cEMTgenesMesUp %>% group_by(WT, CA1A) %>% summarise(n = sum(N)) -> FigS2cEMTgenesMesUp
pdf(file = file.path(dataDir, "FigS2cEMTgenesMesUp.pdf"), paper = "a4r")
alluvial(FigS2cEMTgenesMesUp[,c(1:2)], freq = FigS2cEMTgenesMesUp$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# write tables for EMT genes
write.csv(table("WT" = emtData$wt.group, "TGFb" = emtData$tgfb.group), 
          file =  file.path(dataDir, "Fig2_WT_TGFb_counts_EMTall.csv"))
write.csv(table("WT" = emtData$wt.group, "shZ" = emtData$shZ.group),
          file =  file.path(dataDir, "Fig2_WT_shZ_counts_EMTall.csv"))
write.csv(table("WT" = emtData$wt.group, "CA1A" = emtData$CA1a.group),
          file =  file.path(dataDir, "Fig2_WT_CA1A_counts_EMTall.csv"))

write.csv(table("WT" = emtData[epi_mes == "epi"]$wt.group, "TGFb" = emtData[epi_mes == "epi"]$tgfb.group), 
          file =  file.path(dataDir, "Fig2_WT_TGFb_counts_EMTepi.csv"))
write.csv(table("WT" = emtData[epi_mes == "epi"]$wt.group, "shZ" = emtData[epi_mes == "epi"]$shZ.group),
          file =  file.path(dataDir, "Fig2_WT_shZ_counts_EMTepi.csv"))
write.csv(table("WT" = emtData[epi_mes == "epi"]$wt.group, "CA1A" = emtData[epi_mes == "epi"]$CA1a.group),
          file =  file.path(dataDir, "Fig2_WT_CA1A_counts_EMTepi.csv"))

write.csv(table("WT" = emtData[epi_mes == "mes"]$wt.group, "TGFb" = emtData[epi_mes == "mes"]$tgfb.group), 
          file =  file.path(dataDir, "Fig2_WT_TGFb_counts_EMTmes.csv"))
write.csv(table("WT" = emtData[epi_mes == "mes"]$wt.group, "shZ" = emtData[epi_mes == "mes"]$shZ.group),
          file =  file.path(dataDir, "Fig2_WT_shZ_counts_EMTmes.csv"))
write.csv(table("WT" = emtData[epi_mes == "mes"]$wt.group, "CA1A" = emtData[epi_mes == "mes"]$CA1a.group),
          file =  file.path(dataDir, "Fig2_WT_CA1A_counts_EMTmes.csv"))


