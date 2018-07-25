# Figure S6
rm(list = ls())
# load expression data
RDAs <- list.files("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /", 
                   pattern = "*.rda", 
                   full.names = T)
for (i in RDAs){
  load(i)
}

lapply(list(rT.MCF10Ca1a, rT.shH2AZ, rT.TGFbD6, kT1), function(x) {
  setkey(x, "target_id")
})

# load ChIP-Seq cluster data
dataDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/"
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


dt6$left.MCF10A_WT.category <- as.factor(mcf10awtCategories$cat[match(dt6$MCF10A_WT, mcf10awtCategories$group)])
dt6$right.MCF10A_TGFb.category <- as.factor(mcf10atgfbCategories$cat[match(dt6$MCF10A_TGFb, mcf10atgfbCategories$group)])
dt6$right.MCF10CA1A.category <- as.factor(mcf10ca1aCategories$cat[match(dt6$MCF10CA1A, mcf10ca1aCategories$group)])

FigS6Data <- dt6
save(FigS6Data, file = file.path(dataDir, "FigS6Data.rda"))
write.csv(FigS6Data, file = file.path(dataDir, "FigS6Data.csv"))


# FigS6a plotting ---------------------------------------------------------
FigS6a <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10A_TGFb.category"),
                             data = dt6[!is.na(left.MCF10A_WT.category) & !is.na(right.MCF10A_TGFb.category)],
                             label.size = 4,
                             text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure 6 - H2A.Z sensitivity)") +
  theme(legend.position = "none")

FigS6b <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10CA1A.category"),
                             data = dt6[!is.na(left.MCF10A_WT.category) & !is.na(right.MCF10CA1A.category)],
                             label.size = 4,
                             text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> CA1A (Figure 6 - H2A.Z sensitivity)") +
  theme(legend.position = "none")

ggsave(FigS6a,
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6a_parallelPlots.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

ggsave(FigS6b,
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6b_parallelPlots.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

write.csv(table("WT" = dt6$left.MCF10A_WT.category, "TGFb" = dt6$right.MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/Fig6_WT_TGFb_counts.csv")
write.csv(table("WT" = dt6$left.MCF10A_WT.category, "CA1A" = dt6$right.MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/Fig6_WT_CA1A_counts.csv")

# re-run analysis with EMT genes only -------------------------------------
sigEMTCells <- data.table::fread("~/Data/References/Annotations/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_cellLine.txt")
setkey(dt6, "extGene")
emtData <- dt6[sigEMTCells$cellLine_sig][!is.na(left.MCF10A_WT.category)]
emtData <- merge(emtData, 
                 sigEMTCells, 
                 by.x = "extGene", 
                 by.y = "cellLine_sig", 
                 all.x = T, 
                 all.y = F)

# for FigS6a
FigS6aEMTgenesAll <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10A_TGFb.category"),
                                            data = emtData,
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) +
  ggtitle("MCF10A WT -> TGFb (Figure 6 - all EMT genes)") +
  theme(legend.position = "none")
FigS6aEMTgenesAll$layers[[3]]$data$labels <- ""
FigS6aEMTgenesAll$layers[[4]]$aes_params$colour <- "black"

FigS6aEMTgenesEpi <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10A_TGFb.category"),
                                            data = emtData[epi_mes == "epi"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure 6 - epithelial EMT genes)") +
  theme(legend.position = "none")
FigS6aEMTgenesEpi$layers[[3]]$data$labels <- ""
FigS6aEMTgenesEpi$layers[[4]]$aes_params$colour <- "black"

FigS6aEMTgenesEpiDown <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                     "right.MCF10A_TGFb.category"),
                                            data = rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure 6 - epithelial EMT genes, down-regulated)") +
  theme(legend.position = "none")
FigS6aEMTgenesEpiDown$layers[[3]]$data$labels <- ""
FigS6aEMTgenesEpiDown$layers[[4]]$aes_params$colour <- "black"

# Mesenchymal genes
FigS6aEMTgenesMes <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10A_TGFb.category"),
                                            data = emtData[epi_mes == "mes"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure 6 - mesenchymal EMT genes)") +
  theme(legend.position = "none")
FigS6aEMTgenesMes$layers[[3]]$data$labels <- ""
FigS6aEMTgenesMes$layers[[4]]$aes_params$colour <- "black"

# up-regulated mesenchymal genes
FigS6aEMTgenesMesUp <- ggparallel::ggparallel(list("left.MCF10A_WT.category", 
                                                   "right.MCF10A_TGFb.category"),
                                            data =  rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure 6 - mesenchymal EMT genes, up-regulated)") +
  theme(legend.position = "none")
FigS6aEMTgenesMesUp$layers[[3]]$data$labels <- ""
FigS6aEMTgenesMesUp$layers[[4]]$aes_params$colour <- "black"

FigS6aEMTgenes <- gridExtra::grid.arrange(FigS6aEMTgenesAll, FigS6aEMTgenesEpi, FigS6aEMTgenesMes, nrow = 3, ncol = 1)

FigS6aEMTgenesDE <- gridExtra::grid.arrange(FigS6aEMTgenesEpi, 
                                            FigS6aEMTgenesEpiDown,
                                            FigS6aEMTgenesMes,
                                            FigS6aEMTgenesMesUp,
                                            nrow = 2, ncol = 2)

ggsave(plot = FigS6aEMTgenes, 
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6a_parallelPlots_EMT.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

ggsave(plot = FigS6aEMTgenesDE, 
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6a_parallelPlots_EMT_DiffExp.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")


# for FigS6b ------------------------------------------------------------
setkey(rT.MCF10Ca1a)

# all EMT genes
FigS6bEMTgenesAll <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10CA1A.category"),
                                            data = emtData,
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) +
  ggtitle("MCF10A WT -> TGFb (Figure 6 - all EMT genes)") +
  theme(legend.position = "none")
FigS6bEMTgenesAll$layers[[3]]$data$labels <- ""
FigS6bEMTgenesAll$layers[[4]]$aes_params$colour <- "black"

# all epithelial genes
FigS6bEMTgenesEpi <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10CA1A.category"),
                                            data = emtData[epi_mes == "epi"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> MCF10Ca1a (Figure 6 - all EMT genes)") +
  theme(legend.position = "none")
FigS6bEMTgenesEpi$layers[[3]]$data$labels <- ""
FigS6bEMTgenesEpi$layers[[4]]$aes_params$colour <- "black"

# down-regulated epithelial genes
FigS6bEMTgenesEpiDown <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10CA1A.category"),
                                            data = rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> MCF10Ca1a (Figure 6 - epithelial EMT genes, down-regulated)") +
  theme(legend.position = "none")
FigS6bEMTgenesEpiDown$layers[[3]]$data$labels <- ""
FigS6bEMTgenesEpiDown$layers[[4]]$aes_params$colour <- "black"

# all mesenchymal genes
FigS6bEMTgenesMes <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10CA1A.category"),
                                            data = emtData[epi_mes == "mes"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> MCF10Ca1a (Figure 6 - mesenchymal EMT genes)") +
  theme(legend.position = "none")
FigS6bEMTgenesMes$layers[[3]]$data$labels <- ""
FigS6bEMTgenesMes$layers[[4]]$aes_params$colour <- "black"

# up-regulated mesenchymal genes
FigS6bEMTgenesMesUp <- ggparallel::ggparallel(list("left.MCF10A_WT.category", "right.MCF10CA1A.category"),
                                            data = rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> MCF10Ca1a (Figure 6 - mesenchymal EMT genes, up-regulated)") +
  theme(legend.position = "none")
FigS6bEMTgenesMesUp$layers[[3]]$data$labels <- ""
FigS6bEMTgenesMesUp$layers[[4]]$aes_params$colour <- "black"

# saving
FigS6bEMTgenes <- gridExtra::grid.arrange(FigS6bEMTgenesAll, 
                                          FigS6bEMTgenesEpi, 
                                          FigS6bEMTgenesMes, nrow = 3, ncol = 1)
ggsave(plot = FigS6bEMTgenes, 
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6b_parallelPlots_EMT.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

FigS6bEMTgenesDE <- gridExtra::grid.arrange(FigS6bEMTgenesEpi, 
                                            FigS6bEMTgenesEpiDown,
                                            FigS6bEMTgenesMes, 
                                            FigS6bEMTgenesMesUp,
                                            nrow = 2, ncol = 2)
ggsave(plot = FigS6bEMTgenesDE, 
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6b_parallelPlots_EMT_DiffExp.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

# write tables for EMT genes
write.csv(table("WT" = emtData$left.MCF10A_WT.category, "TGFb" = emtData$right.MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_TGFb_counts_EMTall.csv")
write.csv(table("WT" = emtData$left.MCF10A_WT.category, "CA1A" = emtData$right.MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_CA1A_counts_EMTall.csv")

write.csv(table("WT" = emtData[epi_mes == "epi"]$left.MCF10A_WT.category, "TGFb" = emtData[epi_mes == "epi"]$right.MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_TGFb_counts_EMTepi.csv")
write.csv(table("WT" = emtData[epi_mes == "epi"]$left.MCF10A_WT.category, "CA1A" = emtData[epi_mes == "epi"]$right.MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_CA1A_counts_EMTepi.csv")

write.csv(table("WT" = emtData[epi_mes == "mes"]$left.MCF10A_WT.category, "TGFb" = emtData[epi_mes == "mes"]$right.MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_TGFb_counts_EMTmes.csv")
write.csv(table("WT" = emtData[epi_mes == "mes"]$left.MCF10A_WT.category, "CA1A" = emtData[epi_mes == "mes"]$right.MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_CA1A_counts_EMTmes.csv")

write.csv(table("WT" = emtData[epi_mes == "epi"]$left.MCF10A_WT.category, "TGFb" = emtData[epi_mes == "epi"]$right.MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_TGFb_counts_EMTepi.csv")
write.csv(table("WT" = emtData[epi_mes == "epi"]$left.MCF10A_WT.category, "CA1A" = emtData[epi_mes == "epi"]$right.MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_CA1A_counts_EMTepi.csv")

write.csv(table("WT" = rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0]$left.MCF10A_WT.category, "TGFb" = rT.TGFbD6[emtData[epi_mes == "mes"]][b > 0]$right.MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_TGFb_counts_EMTmes_Up.csv")
write.csv(table("WT" = rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0]$left.MCF10A_WT.category, "CA1A" =  rT.MCF10Ca1a[emtData[epi_mes == "mes"]][b > 0]$right.MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_CA1A_counts_EMTmes_Up.csv")

write.csv(table("WT" = rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0]$left.MCF10A_WT.category, "TGFb" = rT.TGFbD6[emtData[epi_mes == "epi"]][b < 0]$right.MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_TGFb_counts_EMTepi_Down.csv")
write.csv(table("WT" = rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0]$left.MCF10A_WT.category, "CA1A" = rT.MCF10Ca1a[emtData[epi_mes == "epi"]][b < 0]$right.MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 6/FigS6_WT_CA1A_counts_EMTepi_Down.csv")

