## local functions
# Figure S2 (supplement for Figure 2)
# remove all variables from environment
rm(list = ls())


# load expression data ----------------------------------------------------
RDAs <- list.files("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /", 
                   pattern = "*.rda", 
                   full.names = T)
for (i in RDAs){
  load(i)
}
setkey(rT.MCF10Ca1a, "target_id")
setkey(rT.shH2AZ, "target_id")
setkey(rT.TGFbD6, "target_id")
setkey(kT1, "target_id")

# load parallel plot data -------------------------------------------------
dataDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/"
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
dt2 <- l2$log2ratio_10A[!duplicated(extGene), c("extGene", "left.MCF10A_WT.category")]
dt2 <- merge(dt2, l2$log2ratio_TGFb[!duplicated(extGene), c("extGene", "right.MCF10A_TGFb.category")], by.x = "extGene", by.y = "extGene")
dt2 <- merge(dt2, l2$log2ratio_shZ[!duplicated(extGene), c("extGene", "right.MCF10A_shZ.category")], by.x = "extGene", by.y = "extGene")
dt2 <- merge(dt2, l2$log2ratio_CA1a[!duplicated(extGene), c("extGene", "right.MCF10CA1A.category")], by.x = "extGene", by.y = "extGene")
colnames(dt2)

FigS2Data <- dt2
save(FigS2Data, file = file.path(dataDir, "FigS2Data.rda"))
write.csv(FigS2Data, file = file.path(dataDir, "FigS2Data.csv"))

#
apply(dt2[,c("left.MCF10A_WT.category", "right.MCF10A_TGFb.category", "right.MCF10A_shZ.category", "right.MCF10CA1A.category")], 2, function(x) table(x, useNA = "always"))


# Parallel plots ----------------------------------------------------------
FigS2a <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10A_TGFb.category"), #, "right.MCF10A_shZ.category", "right.MCF10CA1A.category"), 
                       data = dt2[!is.na(right.MCF10CA1A.category) & !is.na(right.MCF10A_TGFb.category) ],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F) +
  ggtitle("MCF10A WT -> TGFb (Figure 2)") +
  theme(legend.position = "none")

FigS2b <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10A_shZ.category"), #, , "right.MCF10CA1A.category"), 
                       data = dt2[!is.na(right.MCF10CA1A.category) & !is.na(right.MCF10A_shZ.category) ],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F) +
  ggtitle("MCF10A WT -> shH2AZ (Figure 2)") +
  theme(legend.position = "none")

FigS2c <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10CA1A.category"), #, , ), 
                       data = dt2[!is.na(right.MCF10CA1A.category) & !is.na(right.MCF10CA1A.category) ],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F, text.offset = 0.1) +
  ggtitle("MCF10A WT -> CA1A (Figure 2)") +
  theme(legend.position = "none")

lapply(list(FigS2a, FigS2b, FigS2c), function(x){
  x$layers[[3]]$data$labels <- ""
  x$layers[[4]]$aes_params$colour <- "black"
})

ggsave(plot = FigS2a, 
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/FigS2a_parallelPlots.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

ggsave(plot = FigS2b, 
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/FigS2b_parallelPlots.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

ggsave(plot = FigS2c, 
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/FigS2c_parallelPlots.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

write.csv(table("WT" = dt2$left.MCF10A_WT.category, "TGFb" = dt2$right.MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_TGFb_counts.csv")
write.csv(table("WT" = dt2$left.MCF10A_WT.category, "shZ" = dt2$right.MCF10A_shZ.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_shZ_counts.csv")
write.csv(table("WT" = dt2$left.MCF10A_WT.category, "CA1A" = dt2$right.MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_CA1A_counts.csv")



# re-run analysis with EMT genes only -------------------------------------
sigEMTCells <- data.table::fread("~/Data/References/Annotations/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_cellLine.txt")
setkey(dt2, "extGene")
emtData <- dt2[sigEMTCells$cellLine_sig][!is.na(left.MCF10A_WT.category)]
emtData <- merge(emtData, 
                 sigEMTCells, 
                 by.x = "extGene", 
                 by.y = "cellLine_sig", 
                 all.x = T, 
                 all.y = F)


# FigS2a WT -> TGFb all EMT genes ----------------------------------
FigS2aEMTgenesAll <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                 "right.MCF10A_TGFb.category"),
                                            data = emtData,
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) +
  ggtitle("MCF10A WT -> TGFb (Figure 2 - all EMT genes)") +
  theme(legend.position = "none")
FigS2aEMTgenesAll$layers[[3]]$data$labels <- ""
FigS2aEMTgenesAll$layers[[4]]$aes_params$colour <- "black"


# FigS2a WT -> TGFb epithelial genes --------------------------------------
FigS2aEMTgenesEpi <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                 "right.MCF10A_TGFb.category"),
                                            data = emtData[epi_mes == "epi"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure 2 - epithelial EMT genes)") +
  theme(legend.position = "none")
FigS2aEMTgenesEpi$layers[[3]]$data$labels <- ""
FigS2aEMTgenesEpi$layers[[4]]$aes_params$colour <- "black"


# FigS2a WT -> TGFb down-regulated epithelial genes -----------------------
FigS2aEMTgenesEpiDown <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                     "right.MCF10A_TGFb.category"),
                                              data = rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0],
                                              label.size = 4,
                                              text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure 2 - epithelial EMT genes, down-regulated)") +
  theme(legend.position = "none")
FigS2aEMTgenesEpiDown$layers[[3]]$data$labels <- ""
FigS2aEMTgenesEpiDown$layers[[4]]$aes_params$colour <- "black"


# FigS2a WT -> TGFb all mesenchymal genes ---------------------------------
FigS2aEMTgenesMes <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                 "right.MCF10A_TGFb.category"),
                                            data = emtData[epi_mes == "mes"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure 2 - mesenchymal EMT genes)") +
  theme(legend.position = "none")
FigS2aEMTgenesMes$layers[[3]]$data$labels <- ""
FigS2aEMTgenesMes$layers[[4]]$aes_params$colour <- "black"

# FigS2a WT -> TGFb up-regulated mesenchymal genes ---------------------------------
FigS2aEMTgenesMesUp <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                   "right.MCF10A_TGFb.category"),
                                              data = rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0],
                                              label.size = 4,
                                              text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure 2 - mesenchymal EMT genes, up-regulated)") +
  theme(legend.position = "none")
FigS2aEMTgenesMesUp$layers[[3]]$data$labels <- ""
FigS2aEMTgenesMesUp$layers[[4]]$aes_params$colour <- "black"

# FigS2a save DiffExp plots ---------------------------------
FigS2aEMTgenes <- gridExtra::grid.arrange(FigS2aEMTgenesEpi,
                                          FigS2aEMTgenesEpiDown,
                                          FigS2aEMTgenesMes,
                                          FigS2aEMTgenesMesUp,
                                          nrow = 2, ncol = 2)
ggsave(plot = FigS2aEMTgenes, 
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/FigS2a_parallelPlots_EMT_DiffExp.pdf", 
       device = "pdf",
       height = 450,
       width = 600,
       units = "mm")

# FigS2b WT -> shH2AZ - all EMT genes -------------------------------------
FigS2bEMTgenesAll <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                 "right.MCF10A_shZ.category"),
                                            data = emtData,
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) +
  ggtitle("MCF10A WT -> shH2AZ (Figure 2 - all EMT genes)") +
  theme(legend.position = "none")
FigS2bEMTgenesAll$layers[[3]]$data$labels <- ""
FigS2bEMTgenesAll$layers[[4]]$aes_params$colour <- "black"

# FigS2b WT -> shH2AZ - all epithelial genes -------------------------------------
FigS2bEMTgenesEpi <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                 "right.MCF10A_shZ.category"),
                                            data = emtData[epi_mes == "epi"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> shH2AZ (Figure 2 - epithelial EMT genes)") +
  theme(legend.position = "none")
FigS2bEMTgenesEpi$layers[[3]]$data$labels <- ""
FigS2bEMTgenesEpi$layers[[4]]$aes_params$colour <- "black"

# FigS2b WT -> shH2AZ - down-regulated epithelial genes -------------------------------------
FigS2bEMTgenesEpiDown <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                     "right.MCF10A_shZ.category"),
                                            data = rT.shH2AZ[emtData[epi_mes == "epi"]][b < 0],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> shH2AZ (Figure 2 - down-regulated epithelial EMT genes)") +
  theme(legend.position = "none")
FigS2bEMTgenesEpiDown$layers[[3]]$data$labels <- ""
FigS2bEMTgenesEpiDown$layers[[4]]$aes_params$colour <- "black"

# FigS2b WT -> shH2AZ - mesenchymal genes -------------------------------------
FigS2bEMTgenesMes <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                 "right.MCF10A_shZ.category"),
                                            data = emtData[epi_mes == "mes"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> shH2AZ (Figure 2 - mesenchymal EMT genes)") +
  theme(legend.position = "none")
FigS2bEMTgenesMes$layers[[3]]$data$labels <- ""
FigS2bEMTgenesMes$layers[[4]]$aes_params$colour <- "black"

# FigS2b WT -> shH2AZ - up-regulated mesenchymal genes -------------------------------------
FigS2bEMTgenesMesUp <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                   "right.MCF10A_shZ.category"),
                                              data = rT.shH2AZ[emtData[epi_mes == "mes"]][b > 0],
                                              label.size = 4,
                                              text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> shH2AZ (Figure 2 - mesenchymal EMT genes, up-regulated)") +
  theme(legend.position = "none")
FigS2bEMTgenesMesUp$layers[[3]]$data$labels <- ""
FigS2bEMTgenesMesUp$layers[[4]]$aes_params$colour <- "black"


FigS2bEMTgenes <- gridExtra::grid.arrange(FigS2bEMTgenesEpi,
                                          FigS2bEMTgenesEpiDown,
                                          FigS2bEMTgenesMes,
                                          FigS2bEMTgenesMesUp,
                                          nrow = 2, ncol = 2)
ggsave(plot = FigS2bEMTgenes, 
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/FigS2b_parallelPlots_EMT_DiffEx.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")


# FigS2c WT -> Ca1a all EMT genes -----------------------------------------
FigS2cEMTgenesAll <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                 "right.MCF10CA1A.category"),
                                            data = emtData,
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) +
  ggtitle("MCF10A WT -> MCF10Ca1a (Figure 2 - all EMT genes)") +
  theme(legend.position = "none")
FigS2cEMTgenesAll$layers[[3]]$data$labels <- ""
FigS2cEMTgenesAll$layers[[4]]$aes_params$colour <- "black"

# FigS2c WT -> Ca1a epithelial EMT genes -----------------------------------------
FigS2cEMTgenesEpi <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                 "right.MCF10CA1A.category"),
                                            data = emtData[epi_mes == "epi"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> MCF10Ca1a (Figure 2 - epithelial EMT genes)") +
  theme(legend.position = "none")
FigS2cEMTgenesEpi$layers[[3]]$data$labels <- ""
FigS2cEMTgenesEpi$layers[[4]]$aes_params$colour <- "black"

# FigS2c WT -> Ca1a mesenchymal EMT genes -----------------------------------------
FigS2cEMTgenesMes <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                 "right.MCF10CA1A.category"),
                                            data = emtData[epi_mes == "mes"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> MCF10Ca1a (Figure 2 - mesenchymal EMT genes)") +
  theme(legend.position = "none")
FigS2cEMTgenesMes$layers[[3]]$data$labels <- ""
FigS2cEMTgenesMes$layers[[4]]$aes_params$colour <- "black"

# FigS2c WT -> Ca1a - up-regulated mesenchymal genes -------------------------------------
FigS2cEMTgenesMesUp <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                 "right.MCF10CA1A.category"),
                                              data = rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0],
                                              label.size = 4,
                                              text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> MCF10Ca1a (Figure 2 - mesenchymal EMT genes, up-regulated)") +
  theme(legend.position = "none")
FigS2cEMTgenesMesUp$layers[[3]]$data$labels <- ""
FigS2cEMTgenesMesUp$layers[[4]]$aes_params$colour <- "black"

# FigS2c WT -> Ca1a - down-regulated epithelial genes -------------------------------------
FigS2cEMTgenesEpiDown <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                   "right.MCF10CA1A.category"),
                                              data = rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0],
                                              label.size = 4,
                                              text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> MCF10Ca1a (Figure 2 - epithelial EMT genes, down-regulated)") +
  theme(legend.position = "none")
FigS2cEMTgenesEpiDown$layers[[3]]$data$labels <- ""
FigS2cEMTgenesEpiDown$layers[[4]]$aes_params$colour <- "black"

FigS2cEMTgenes <- gridExtra::grid.arrange(FigS2cEMTgenesEpi,
                                          FigS2cEMTgenesEpiDown,
                                          FigS2cEMTgenesMes,
                                          FigS2cEMTgenesMesUp,
                                          nrow = 2, ncol = 2)
ggsave(plot = FigS2cEMTgenes, 
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/FigS2c_parallelPlots_EMT_DiffExp.pdf", 
       device = "pdf",
       height = 450,
       width = 600,
       units = "mm")

# write tables for EMT genes
write.csv(table("WT" = emtData$left.MCF10A_WT.category, "TGFb" = emtData$right.MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_TGFb_counts_EMTall.csv")
write.csv(table("WT" = emtData$left.MCF10A_WT.category, "shZ" = emtData$right.MCF10A_shZ.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_shZ_counts_EMTall.csv")
write.csv(table("WT" = emtData$left.MCF10A_WT.category, "CA1A" = emtData$right.MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_CA1A_counts_EMTall.csv")

write.csv(table("WT" = emtData[epi_mes == "epi"]$left.MCF10A_WT.category, "TGFb" = emtData[epi_mes == "epi"]$right.MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_TGFb_counts_EMTepi.csv")
write.csv(table("WT" = emtData[epi_mes == "epi"]$left.MCF10A_WT.category, "shZ" = emtData[epi_mes == "epi"]$right.MCF10A_shZ.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_shZ_counts_EMTepi.csv")
write.csv(table("WT" = emtData[epi_mes == "epi"]$left.MCF10A_WT.category, "CA1A" = emtData[epi_mes == "epi"]$right.MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_CA1A_counts_EMTepi.csv")

write.csv(table("WT" = emtData[epi_mes == "mes"]$left.MCF10A_WT.category, "TGFb" = emtData[epi_mes == "mes"]$right.MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_TGFb_counts_EMTmes.csv")
write.csv(table("WT" = emtData[epi_mes == "mes"]$left.MCF10A_WT.category, "shZ" = emtData[epi_mes == "mes"]$right.MCF10A_shZ.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_shZ_counts_EMTmes.csv")
write.csv(table("WT" = emtData[epi_mes == "mes"]$left.MCF10A_WT.category, "CA1A" = emtData[epi_mes == "mes"]$right.MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 2/Fig2_WT_CA1A_counts_EMTmes.csv")


# add GO analysis of MCF10A WT clusters -----------------------------------
library("clusterProfiler")
library("org.Hs.eg.db")

keyType <- "SYMBOL"
ontologies <- c("MF", "BP", "CC")
pAdjust <- "fdr"
pCutoff <- 0.25
analysis <- "clusterProfiler"

geneUniverse <- dt2$extGene

GO.MCF10A_WT.clusters <- lapply(levels(dt2$left.MCF10A_WT.category), function(x){
  geneList <- dt2[left.MCF10A_WT.category == x]$extGene
  MF <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "MF",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  BP <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "BP",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  CC <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "CC",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  return(list(MF = MF, CC = CC, BP = BP))
})

names(GO.MCF10A_WT.clusters) <- levels(dt2$left.MCF10A_WT.category)
