# Figure S6
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


# load ChIP-Seq cluster data ----------------------------------------------
dataDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/Figure_6/"
l6 <- lapply(list.files(path = dataDir, pattern=".tsv"), function(x){
  dt <- data.table::fread(paste(dataDir, x, sep = "/"))
  return(dt)
})
n6 <- gsub(x = list.files(path = dataDir, pattern=".tsv"), pattern = ".tsv", replacement =  "") 
names(l6) <- unlist(lapply(strsplit(n6, "_"), function(x) paste(x[2:3], collapse = "_")))

l6 <- lapply(l6, function(x) {
  ucscID <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[4]))
  extGene <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[6]))
  x$ucscID <- ucscID
  x$extGene <- extGene
  x$chr <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[1]))
  x$start <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[2]))
  x$end <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[3]))
  return(x)
})

lapply(names(l6), function(x){
  setkey(l6[[x]], "extGene")
})

dt6 <- merge(subset(l6$log2ratio_10A[!duplicated("extGene")], select = c("extGene", "group1")), 
             subset(l6$log2ratio_TGFb[!duplicated("extGene")], select = c("extGene", "group1")), 
             by.x = "extGene", by.y = "extGene")
colnames(dt6)[2:3] <- c("MCF10A_WT", "MCF10A_TGFb")
dt6 <- merge(dt6, 
             subset(l6$log2ratio_CA1a[!duplicated(extGene)], select = c("extGene", "group1")), 
             by.x = "extGene", by.y = "extGene")
colnames(dt6)[4] <- c("MCF10CA1A")

# For H2A.Z sensitivity
# 1. sensitive +2 ().
# 2. resistant TSS nucleosome ().
# 3. sensitive -1/+1 nucleosome ().
# 4. resistant fuzzy +1/+2 nucleosome ().
# 5. neither sensitive or resistant ()
# 6. Sensitive -2 ().
# 7. Resistant -2 ().
# 8. Sensitive +1 and +2 ().
# 9. Sensitive -1 ().
# 10. Sensitive -1 and -2 ()
# 11. Sensitive +1 ()


# Ca1a
# 1 - "Sensitive +2"
# 2 - "Resistant -2"
# 5 - "resistant TSS nucleosome"
# 6 - "resistant fuzzy +1/+2 nucleosome"
# 3 - "Sensitive -1 and -2"
# 4 - "neither sensitive or resistant"
# 7 - "Sensitive +1"

mcf10ca1aCategories <- data.table::data.table(cat = c(
                                                      "Sensitive +2",
                                                      "Resistant -2",
                                                      "Resistant TSS nucleosome",
                                                      "Resistant fuzzy +1/+2 nucleosome",
                                                      "Sensitive -1 and -2",
                                                      "Neither sensitive or resistant",
                                                      "Sensitive +1"
                                                    ),
                                              group = c(1,
                                                        2,
                                                        3,
                                                        4,
                                                        5,
                                                        6,
                                                        7))

mcf10ca1aCategories$cat <- factor(mcf10ca1aCategories$cat, levels = mcf10ca1aCategories$cat, ordered = T)
mcf10ca1aCategories <- merge(mcf10ca1aCategories, unique(l6$log2ratio_CA1a[, c("group1", "color")]), by.x = "group", by.y = "group1", all.x = F, all.y = F)

# TGFb
# 1 - "Sensitive +1 and +2"
# 2 - "Sensitive -1"
# 3 - "Resistant -2"
# 4 - "Resistant TSS nucleosome"
# 5 - "Resistant fuzzy +1/+2 nucleosome"
# 6 - "Sensitive -2"
# 7 - "Neither sensitive or resistant"
mcf10atgfbCategories <- data.table::data.table(cat = c(
                                                    "Sensitive +1 and +2",
                                                    "Sensitive -1",
                                                    "Resistant -2",
                                                    "Resistant TSS nucleosome",
                                                    "Resistant fuzzy +1/+2 nucleosome",
                                                    "Sensitive -2",
                                                    "Neither sensitive or resistant"
                                                  ),
                                               group = c(1,
                                                         2,
                                                         3,
                                                         4,
                                                         5,
                                                         6,
                                                         7))
mcf10atgfbCategories$cat <- factor(mcf10atgfbCategories$cat, levels = mcf10atgfbCategories$cat, ordered = T)
mcf10atgfbCategories <- merge(mcf10atgfbCategories, unique(l6$log2ratio_TGFb[, c("group1", "color")]), by.x = "group", by.y = "group1", all.x = F, all.y = F)

# For H2A.Z sensitivity
# 1 - "sensitive +2"
# 2 - "resistant TSS nucleosome"
# 3 - "sensitive -1/+1 nucleosome"
# 4 - "resistant fuzzy +1/+2 nucleosome"
# 5 - "neither sensitive or resistant"
# 6 - "Sensitive -2"
# 7 - "Resistant -2"
mcf10awtCategories <- data.table::data.table(cat = c("Sensitive +2",
                                                     "Resistant TSS nucleosome",
                                                     "Sensitive -1/+1 nucleosome",
                                                     "Resistant fuzzy +1/+2 nucleosome",
                                                     "Neither sensitive or resistant",
                                                     "Sensitive -2",
                                                     "Resistant -2"),
                                             group = c(1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6,
                                                       7))
mcf10awtCategories$cat <- factor(mcf10awtCategories$cat, levels = mcf10awtCategories$cat, ordered = T)
mcf10awtCategories <- merge(mcf10awtCategories, unique(l6$log2ratio_10A[, c("group1", "color")]), by.x = "group", by.y = "group1", all.x = F, all.y = F)

dt6$left.MCF10A_WT.category <- as.factor(mcf10awtCategories$cat[match(dt6$MCF10A_WT, mcf10awtCategories$group)])
dt6$right.MCF10A_TGFb.category <- as.factor(mcf10atgfbCategories$cat[match(dt6$MCF10A_TGFb, mcf10atgfbCategories$group)])
dt6$right.MCF10CA1A.category <- as.factor(mcf10ca1aCategories$cat[match(dt6$MCF10CA1A, mcf10ca1aCategories$group)])

FigS6Data <- dt6
rm(dt6)
save(FigS6Data, file = file.path(dataDir, "FigS6Data.rda"))
write.csv(FigS6Data, file = file.path(dataDir, "FigS6Data.csv"))

# FigS6a - alluvial plot WT -> TGFb --------------------------------
tab1 <- data.table::as.data.table(table(FigS6Data$MCF10A_WT, FigS6Data$MCF10A_TGFb))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6aTab <- data.table::as.data.table(table("WT" = FigS6Data$MCF10A_WT, "TGFb" = FigS6Data$MCF10A_TGFb))
FigS6aTab %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS6aTab1
pdf(file = file.path(dataDir, "AlluvialPlotS6a.pdf"), paper = "a4r")
alluvial(FigS6aTab[,c(1:2)], freq = FigS6aTab$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()


# FigS6b - alluvial plot WT -> CA1A ----------------------------------------
tab1 <- data.table::as.data.table(table(FigS6Data$MCF10A_WT, FigS6Data$MCF10CA1A))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6bTab <- data.table::as.data.table(table("WT" = FigS6Data$MCF10A_WT, "CA1A" = FigS6Data$MCF10CA1A))
FigS6bTab %>% group_by(WT, CA1A) %>% summarise(n = sum(N)) -> FigS6bTab
pdf(file = file.path(dataDir, "AlluvialPlotS6b.pdf"), paper = "a4r")
alluvial(FigS6bTab[,c(1:2)], freq = FigS6bTab$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# re-run analysis with EMT genes only - data preparation -----------
sigEMTCells <- data.table::fread("/Data/References/Annotations/GeneSets/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_cellLine.txt")
sigEMTTumor <- data.table::fread("/Data/References/Annotations/GeneSets/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_tumor.txt")
colnames(sigEMTCells)[1] <- "gene_symbol"
colnames(sigEMTTumor)[1] <- "gene_symbol"
setkey(sigEMTTumor, "gene_symbol")
sigEMT <- rbind(sigEMTCells,
                sigEMTTumor[!intersect(sigEMTCells$gene_symbol, sigEMTTumor$gene_symbol)])
setkey(FigS6Data, "extGene")
emtData <- FigS6Data[sigEMT$gene_symbol][!is.na(MCF10A_WT)]
emtData <- merge(emtData, 
                 sigEMT, 
                 by.x = "extGene", 
                 by.y = "gene_symbol", 
                 all.x = T, 
                 all.y = F)

# FigS6a WT -> TGFb all EMT genes ----------------------------------
tab1 <- data.table::as.data.table(table(emtData$MCF10A_WT, emtData$MCF10A_TGFb))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6aTabEMTAll <- data.table::as.data.table(table("WT" = emtData$MCF10A_WT, "TGFb" = emtData$MCF10A_TGFb))
FigS6aTabEMTAll %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS6aTabEMTAll1
pdf(file = file.path(dataDir, "FigS6aTabEMTAll.pdf"), paper = "a4r")
alluvial(FigS6aTabEMTAll1[,c(1:2)], freq = FigS6aTabEMTAll1$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()


# FigS6a WT -> TGFb epithelial genes --------------------------------------
tab1 <- data.table::as.data.table(table(emtData[epi_mes == "epi"]$MCF10A_WT, 
                                        emtData[epi_mes == "epi"]$MCF10A_TGFb))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6aEMTgenesEpi <- data.table::as.data.table(table("WT" = emtData[epi_mes == "epi"]$MCF10A_WT, 
                                                     "TGFb" = emtData[epi_mes == "epi"]$MCF10A_TGFb))
FigS6aEMTgenesEpi %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS6aEMTgenesEpi
pdf(file = file.path(dataDir, "FigS6aEMTgenesEpi.pdf"), paper = "a4r")
alluvial(FigS6aEMTgenesEpi[,c(1:2)], freq = FigS6aEMTgenesEpi$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS6a WT -> TGFb down-regulated epithelial genes -----------------------
tab1 <- data.table::as.data.table(table(rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0]$MCF10A_WT, 
                                        rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0]$MCF10A_TGFb))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6aEMTgenesEpiDown <- data.table::as.data.table(table("WT" = rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0]$MCF10A_WT, 
                                                         "TGFb" = rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0]$MCF10A_TGFb))
FigS6aEMTgenesEpiDown %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS6aEMTgenesEpiDown
pdf(file = file.path(dataDir, "FigS6aEMTgenesEpiDown.pdf"), paper = "a4r")
alluvial(FigS6aEMTgenesEpiDown[,c(1:2)], freq = FigS6aEMTgenesEpiDown$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS6a WT -> TGFb all mesenchymal genes ---------------------------------
tab1 <- data.table::as.data.table(table(emtData[epi_mes == "mes"]$MCF10A_WT, 
                                        emtData[epi_mes == "mes"]$MCF10A_TGFb))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6aEMTgenesMes <- data.table::as.data.table(table("WT" = emtData[epi_mes == "mes"]$MCF10A_WT, 
                                                     "TGFb" = emtData[epi_mes == "mes"]$MCF10A_TGFb))
FigS6aEMTgenesMes %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS6aEMTgenesMes
pdf(file = file.path(dataDir, "FigS6aEMTgenesMes.pdf"), paper = "a4r")
alluvial(FigS6aEMTgenesMes[,c(1:2)], freq = FigS6aEMTgenesMes$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS6a WT -> TGFb up-regulated mesenchymal genes ---------------------------------
tab1 <- data.table::as.data.table(table(rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0]$MCF10A_WT, 
                                        rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0]$MCF10A_TGFb))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6aEMTgenesMesUp <- data.table::as.data.table(table("WT" = rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0]$MCF10A_WT, 
                                                       "TGFb" = rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0]$MCF10A_TGFb))
FigS6aEMTgenesMesUp %>% group_by(WT, TGFb) %>% summarise(n = sum(N)) -> FigS6aEMTgenesMesUp
pdf(file = file.path(dataDir, "FigS6aEMTgenesMesUp.pdf"), paper = "a4r")
alluvial(FigS6aEMTgenesMesUp[,c(1:2)], freq = FigS6aEMTgenesMesUp$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS6b WT -> Ca1a all EMT genes -----------------------------------------
tab1 <- data.table::as.data.table(table(emtData$MCF10A_WT, emtData$MCF10CA1A))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6bTabEMTAll <- data.table::as.data.table(table("WT" = emtData$MCF10A_WT, "CA1A" = emtData$MCF10CA1A))
FigS6bTabEMTAll %>% group_by(WT, CA1A) %>% summarise(n = sum(N)) -> FigS6bTabEMTAll
pdf(file = file.path(dataDir, "FigS6bTabEMTAll.pdf"), paper = "a4r")
alluvial(FigS6bTabEMTAll[,c(1:2)], freq = FigS6bTabEMTAll$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS6b WT -> Ca1a epithelial EMT genes -----------------------------------------
tab1 <- data.table::as.data.table(table(emtData[epi_mes == "epi"]$MCF10A_WT, 
                                        emtData[epi_mes == "epi"]$MCF10CA1A))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6bEMTgenesEpi <- data.table::as.data.table(table("WT" = emtData[epi_mes == "epi"]$MCF10A_WT, 
                                                     "CA1A" = emtData[epi_mes == "epi"]$MCF10CA1A))
FigS6bEMTgenesEpi %>% group_by(WT, CA1A) %>% summarise(n = sum(N)) -> FigS2bEMTgenesEpi
pdf(file = file.path(dataDir, "FigS2bEMTgenesEpi.pdf"), paper = "a4r")
alluvial(FigS2bEMTgenesEpi[,c(1:2)], freq = FigS2bEMTgenesEpi$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS6b WT -> Ca1a - down-regulated epithelial genes -------------------------------------
tab1 <- data.table::as.data.table(table(rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0]$MCF10A_WT, 
                                        rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0]$MCF10CA1A))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6bEMTgenesEpiDown <- data.table::as.data.table(table("WT" = rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0]$MCF10A_WT, 
                                                         "CA1A" = rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0]$MCF10CA1A))
FigS6bEMTgenesEpiDown %>% group_by(WT, CA1A) %>% summarise(n = sum(N)) -> FigS6bEMTgenesEpiDown
pdf(file = file.path(dataDir, "FigS6bEMTgenesEpiDown.pdf"), paper = "a4r")
alluvial(FigS6bEMTgenesEpiDown[,c(1:2)], freq = FigS6bEMTgenesEpiDown$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS6b WT -> Ca1a mesenchymal EMT genes -----------------------------------------
tab1 <- data.table::as.data.table(table(emtData[epi_mes == "mes"]$MCF10A_WT, 
                                        emtData[epi_mes == "mes"]$MCF10CA1A))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6bEMTgenesMes <- data.table::as.data.table(table("WT" = emtData[epi_mes == "mes"]$MCF10A_WT, 
                                                     "CA1A" = emtData[epi_mes == "mes"]$MCF10CA1A))
FigS6bEMTgenesMes %>% group_by(WT, CA1A) %>% summarise(n = sum(N)) -> FigS6bEMTgenesMes
pdf(file = file.path(dataDir, "FigS6bEMTgenesMes.pdf"), paper = "a4r")
alluvial(FigS6bEMTgenesMes[,c(1:2)], freq = FigS6bEMTgenesMes$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# FigS6b WT -> Ca1a - up-regulated mesenchymal genes -------------------------------------
tab1 <- data.table::as.data.table(table(rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0]$MCF10A_WT, 
                                        rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0]$MCF10CA1A))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS6bEMTgenesMesUp <- data.table::as.data.table(table("WT" = rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0]$MCF10A_WT, 
                                                       "CA1A" = rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0]$MCF10CA1A))
FigS6bEMTgenesMesUp %>% group_by(WT, CA1A) %>% summarise(n = sum(N)) -> FigS6bEMTgenesMesUp
pdf(file = file.path(dataDir, "FigS6bEMTgenesMesUp.pdf"), paper = "a4r")
alluvial(FigS6bEMTgenesMesUp[,c(1:2)], freq = FigS6bEMTgenesMesUp$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
dev.off()

# write tables for EMT genes
write.csv(table("WT" = emtData$MCF10A_WT, "TGFb" = emtData$MCF10A_TGFb), 
          file =  file.path(dataDir, "Fig2_WT_TGFb_counts_EMTall.csv"))
write.csv(table("WT" = emtData$MCF10A_WT, "CA1A" = emtData$MCF10CA1A),
          file =  file.path(dataDir, "Fig2_WT_CA1A_counts_EMTall.csv"))

write.csv(table("WT" = emtData[epi_mes == "epi"]$MCF10A_WT, "TGFb" = emtData[epi_mes == "epi"]$MCF10A_TGFb), 
          file =  file.path(dataDir, "Fig2_WT_TGFb_counts_EMTepi.csv"))
write.csv(table("WT" = emtData[epi_mes == "epi"]$MCF10A_WT, "CA1A" = emtData[epi_mes == "epi"]$MCF10CA1A),
          file =  file.path(dataDir, "Fig2_WT_CA1A_counts_EMTepi.csv"))

write.csv(table("WT" = emtData[epi_mes == "mes"]$MCF10A_WT, "TGFb" = emtData[epi_mes == "mes"]$MCF10A_TGFb), 
          file =  file.path(dataDir, "Fig2_WT_TGFb_counts_EMTmes.csv"))
write.csv(table("WT" = emtData[epi_mes == "mes"]$MCF10A_WT, "CA1A" = emtData[epi_mes == "mes"]$MCF10CA1A),
          file =  file.path(dataDir, "Fig2_WT_CA1A_counts_EMTmes.csv"))

