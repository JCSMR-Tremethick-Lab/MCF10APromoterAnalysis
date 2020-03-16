# just to make sure that no remnants of previous plot remain
# Figure S1 (supplement to Figure 1)
rm(list = ls())
## local functions
dataDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 1/"
l1 <- lapply(list.files(path = dataDir, pattern=".tsv"), function(x){
  dt <- data.table::fread(paste(dataDir, x, sep = "/"))
  return(dt)
})
n1 <- gsub(x = list.files(path = dataDir, pattern=".tsv"), pattern = ".tsv", replacement =  "") 
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

dt1 <- merge(subset(l1$Total_10A[!duplicated(extGene)], select = c("extGene", "group1")), subset(l1$Total_TGFb[!duplicated(extGene)], select = c("extGene", "group1")), by.x = "extGene", by.y = "extGene")
colnames(dt1)[2:3] <- c("MCF10A_WT", "MCF10A_TGFb")
dt1 <- merge(dt1, subset(l1$Total_shZ[!duplicated(extGene)], select = c("extGene", "group1")), by.x = "extGene", by.y = "extGene")
dt1 <- merge(dt1, subset(l1$Total_CA1a[!duplicated(extGene)], select = c("extGene", "group1")), by.x = "extGene", by.y = "extGene")
colnames(dt1)[4:5] <- c("MCF10A_shZ", "MCF10CA1A")

# For input (figure 1)
# We have decided on the following classes 
# 1. Strong +2 (1 for 10A, 1 for TGFb, 1 for shH2AZ, 1 for CAC1)
# 2. Strong +1 (2 for 10A, 2 for TGFb, 2 for shH2AZ, 2 for CAC1)
# 3. weak upstream (3 for 10A, 5 for TGFb, 5 for shH2AZ, 5 for CAC1)
# 4. weak  array across promoter (4 for 10A, 7 for TGFb, 7 for shH2AZ, 7 for CAC1)
# 5. Strong -2 (5 for 10A, 6 for TGFb, 6 for shH2AZ, 6 for CAC1)
# 6. Strong -1 (6 for 10A, 4 for TGFb, 4 for shH2AZ, 4 for CAC1)
# 7. Overall depleted  (7 for 10A, 3 for TGFb, 3 for shH2AZ, 3 for CAC1)

mcf10awtCategories <- data.table::data.table(cat = c("Strong +2", 
                                                     "Strong +1", 
                                                     "weak upstream", 
                                                     "weak array across promoter",
                                                     "Strong -2", 
                                                     "Strong -1", 
                                                     "Overall depleted"),
                                             group = c(1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6,
                                                       7))

mcf10awtCategories$cat <- factor(mcf10awtCategories$cat, levels = mcf10awtCategories$cat, ordered = T)

# For input (figure 1)
# We have decided on the following classes 
# 1. Strong +2 (1 for 10A, 1 for TGFb, 1 for shH2AZ, 1 for CAC1)
# 2. Strong +1 (2 for 10A, 2 for TGFb, 2 for shH2AZ, 2 for CAC1)
# 3. weak upstream (3 for 10A, 5 for TGFb, 5 for shH2AZ, 5 for CAC1)
# 4. weak  array across promoter (4 for 10A, 7 for TGFb, 7 for shH2AZ, 7 for CAC1)
# 5. Strong -2 (5 for 10A, 6 for TGFb, 6 for shH2AZ, 6 for CAC1)
# 6. Strong -1 (6 for 10A, 4 for TGFb, 4 for shH2AZ, 4 for CAC1)
# 7. Overall depleted  (7 for 10A, 3 for TGFb, 3 for shH2AZ, 3 for CAC1)

mcf10atgfbCategories <- data.table::data.table(cat = c(
                                                       "Strong +2", #1
                                                       "Strong +1", #2
                                                       "Overall depleted", #3
                                                       "Strong -1", #4
                                                       "weak upstream", #5 
                                                       "Strong -2", #6
                                                       "weak array across promoter" #7
                                                       ),
                                               group = c(1,
                                                         2,
                                                         3,
                                                         4,
                                                         5,
                                                         6,
                                                         7))
mcf10atgfbCategories$cat <- factor(mcf10atgfbCategories$cat, levels = mcf10atgfbCategories$cat, ordered = T)

mcf10ashh2azCategories <- data.table::data.table(cat = c(
                                                         "Strong +2", #1 
                                                         "Strong +1", #2
                                                         "Overall depleted", #3
                                                         "Strong -1", #4
                                                         "weak upstream", #5
                                                         "Strong -2", #6
                                                         "weak array across promoter"#7
                                                         ),
                                                 group = c(1,
                                                           2,
                                                           3,
                                                           4,
                                                           5,
                                                           6,
                                                           7))
mcf10ashh2azCategories$cat <- factor(mcf10ashh2azCategories$cat, levels = mcf10ashh2azCategories$cat, ordered = T)

mcf10ca1aCategories <- data.table::data.table(cat = c(
                                                      "Strong +2", #1 
                                                      "Strong +1", #2
                                                      "Overall depleted", #3
                                                      "Strong -1", #4
                                                      "weak upstream", #5
                                                      "Strong -2", #6
                                                      "weak array across promoter" #7
                                                      ),
                                              group = c(1,
                                                        2,
                                                        3,
                                                        4,
                                                        5,
                                                        6,
                                                        7))
mcf10ca1aCategories$cat <- factor(mcf10ca1aCategories$cat, levels = mcf10ca1aCategories$cat, ordered = T)

dt1$MCF10A_WT.category <- as.factor(mcf10awtCategories$cat[match(dt1$MCF10A_WT, mcf10awtCategories$group)])
dt1$MCF10A_TGFb.category <- as.factor(mcf10atgfbCategories$cat[match(dt1$MCF10A_TGFb, mcf10atgfbCategories$group)])
dt1$MCF10A_shZ.category <- as.factor(mcf10ashh2azCategories$cat[match(dt1$MCF10A_shZ, mcf10ashh2azCategories$group)])
dt1$MCF10CA1A.category <- as.factor(mcf10ca1aCategories$cat[match(dt1$MCF10CA1A, mcf10ca1aCategories$group)])

apply(dt1[,c("MCF10A_WT.category", "MCF10A_TGFb.category", "MCF10A_shZ.category", "MCF10CA1A.category")], 2, function(x) table(x, useNA = "always"))

FigS1Data <- dt1
save(FigS1Data, file = file.path(dataDir, "FigS1Data.rda"))
write.csv(FigS1Data, file = file.path(dataDir, "FigS1Data.csv"))


# FigS1a plotting ---------------------------------------------------------


FigS1a <- ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_TGFb.category"),
                       data = dt1[!is.na(MCF10A_WT.category) & !is.na(MCF10A_TGFb.category)],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure1)") 

FigS1b <- ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_shZ.category"),
                       data = dt1[!is.na(MCF10A_WT.category) & !is.na(MCF10A_shZ.category)],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> shH2AZ (Figure1") 
FigS1c <- ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10CA1A.category"),
                       data = dt1[!is.na(MCF10A_WT.category) & !is.na(MCF10CA1A.category)],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> CA1A (Figure 1)") 

ggsave(FigS1a,
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 1/FigS1a_parallelPlots.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

ggsave(FigS1b,
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 1/FigS1b_parallelPlots.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

ggsave(FigS1c,
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 1/FigS1c_parallelPlots.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

write.csv(table("WT" = dt1$MCF10A_WT.category, "TGFb" = dt1$MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 1/Fig1_WT_TGFb_counts.csv")
write.csv(table("WT" = dt1$MCF10A_WT.category, "shH2AZ" = dt1$MCF10A_shZ.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 1/Fig1_WT_shH2AZ_counts.csv")
write.csv(table("WT" = dt1$MCF10A_WT.category, "CA1A" = dt1$MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 1/Fig1_WT_CA1A_counts.csv")


# re-run analysis with EMT genes only -------------------------------------
sigEMTCells <- data.table::fread("~/Data/References/Annotations/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_cellLine.txt")
setkey(dt1, "extGene")
emtData <- dt1[sigEMTCells$cellLine_sig][!is.na(MCF10A_WT.category)]
emtData <- merge(emtData, 
                 sigEMTCells, 
                 by.x = "extGene", 
                 by.y = "cellLine_sig", 
                 all.x = T, 
                 all.y = F)

FigS1aEMTgenesAll <- ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_TGFb.category"),
                                 data = emtData,
                                 label.size = 4,
                                 text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure1 - all EMT genes)") 

FigS1aEMTgenesEpi <- ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_TGFb.category"),
                                            data = emtData[epi_mes == "epi"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure1 - epithelial EMT genes)") 

FigS1aEMTgenesMes <- ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_TGFb.category"),
                                            data = emtData[epi_mes == "mes"],
                                            label.size = 4,
                                            text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure1 - mesenchymal EMT genes)") 
