rm(list = ls())
dataDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 5/"

l5 <- lapply(list.files(path = dataDir, pattern=".tsv"), function(x){
  dt <- data.table::fread(paste(dataDir, x, sep = "/"))
  return(dt)
})
n5 <- gsub(x = list.files(path = dataDir, pattern=".tsv"), pattern = ".tsv", replacement =  "") 
names(l5) <- unlist(lapply(strsplit(n5, "_"), function(x) paste(x[2:3], collapse = "_")))

l5 <- lapply(l5, function(x) {
  ucscID <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[4]))
  extGene <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[6]))
  x$ucscID <- ucscID
  x$extGene <- extGene
  x$chr <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[1]))
  x$start <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[2]))
  x$end <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[3]))
  return(x)
})
l5

rm(dt5)
dt5 <- merge(subset(l5$TOTAL_10A[!duplicated(extGene)], select = c("extGene", "group1")), subset(l5$TOTAL_TGFb[!duplicated(extGene)], select = c("extGene", "group1")), by.x = "extGene", by.y = "extGene")
colnames(dt5)[2:3] <- c("MCF10A_WT", "MCF10A_TGFb")
dt5 <- merge(dt5, subset(l5$TOTAL_CA1a[!duplicated(extGene)], select = c("extGene", "group1")), by.x = "extGene", by.y = "extGene")
colnames(dt5)[4] <- c("MCF10CA1A")
colnames(dt5)
dt5

# For Figure 5. H2A.Z 
# 1.  Weak -1 with weak upstream and downstream (red cluster 1 in 10A, -not seen in other clusters).
# 2. No positioned nucleosomes (yellow cluster 2 in 10A, 5 TGFb, cluster 3 in CA1A).
# 3. Strong -1b (cluster 3 in 10A, cluster 3 in TGFb, cluster 1 CA1A)
# 4. Strong +1 (cluster 4 in 10A, cluster 2 in TGFb, cluster 7 CA1A).
# 5. Weak upstream -1 and -2 (cluster 5 10A, cluster 7 TGFb, cluster 2 CA1A).
# 6. Strong -1a (cluster 6 in 10A, cluster 4A TGFb, cluster 6 and 4 in CAC1) 
# 7. Strong -2 (cluster 7 10A, cluster 6 TGFb, cluster 5 TGFb).

# WT
# 1 - Weak -1 with weak upstream and downstream
# 2 - No positioned nucleosomes
# 3 - Strong -1b
# 4 - Strong +1
# 5 - Weak upstream -1 & -2
# 6 - Strong -1a
# 7 - Strong -2
mcf10awtCategories <- data.table::data.table(cat = c("Weak -1 with weak upstream and downstream", 
                                                     "No positioned nucleosomes", 
                                                     "Strong -1b", 
                                                     "Strong +1", 
                                                     "Weak upstream -1 & -2", 
                                                     "Strong -1a",
                                                     "Strong -2"),
                                             group = c(1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6,
                                                       7))
mcf10awtCategories$cat <- factor(mcf10awtCategories$cat, levels = mcf10awtCategories$cat, ordered = T)

# TGFb
# 1 - Weak/fuzzy downstream
# 2 - Strong +1
# 3 - Strong -1b
# 4 - Strong -1a
# 5 - No positioned nucleosomes
# 6 - Strong -2
# 7 - Weak upstream -1 & -2
mcf10atgfbCategories <- data.table::data.table(cat = c("Weak/fuzzy downstream",
                                                       "Strong +1",
                                                       "Strong -1b",
                                                       "Strong -1a",
                                                       "No positioned nucleosomes",
                                                       "Strong -2",
                                                       "Weak upstream -1 & -2"),
                                               group = c(1,
                                                         2,
                                                         3,
                                                         4,
                                                         5,
                                                         6,
                                                         7))
mcf10atgfbCategories$cat <- factor(mcf10atgfbCategories$cat, levels = mcf10atgfbCategories$cat, ordered = T)

# Ca1a
# 1 - Strong -1b
# 2 - Weak upstream -1 & -2
# 3 - No positioned nucleosomes
# 4 - Strong -1a
# 5 - Strong -2
# 6 - Strong -1c
# 7 - Strong +1

mcf10ca1aCategories <- data.table::data.table(cat = c("Strong -1b",
                                                      "Weak upstream -1 & -2",
                                                      "No positioned nucleosomes",
                                                      "Strong -1a",
                                                      "Strong -2",
                                                      "Strong -1c",
                                                      "Strong +1"),
                                              group = c(1,
                                                        2,
                                                        3,
                                                        4,
                                                        5,
                                                        6,
                                                        7))
mcf10ca1aCategories$cat <- factor(mcf10ca1aCategories$cat, levels = mcf10ca1aCategories$cat, ordered = T)

dt5$MCF10A_WT.category <- as.factor(mcf10awtCategories$cat[match(dt5$MCF10A_WT, mcf10awtCategories$group)])
dt5$MCF10A_TGFb.category <- as.factor(mcf10atgfbCategories$cat[match(dt5$MCF10A_TGFb, mcf10atgfbCategories$group)])
dt5$MCF10CA1A.category <- as.factor(mcf10ca1aCategories$cat[match(dt5$MCF10CA1A, mcf10ca1aCategories$group)])

sapply(dt5[, -c("extGene")], function(x) table(x, useNA = "always"))

FigS5Data <- dt5
save(FigS5Data, file = file.path(dataDir, "FigS5Data.rda"))
write.csv(FigS5Data, file = file.path(dataDir, "FigS5Data.csv"))

# FigS5a plotting ---------------------------------------------------------
FigS5a <- ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_TGFb.category"),
                       data = dt5[!is.na(MCF10A_WT.category) & !is.na(MCF10A_TGFb.category)],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> TGFb (Figure4/5)") +
  theme(legend.position = "none")

FigS5b <- ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10CA1A.category"), #, , ), 
                       data = dt5[!is.na(MCF10A_WT.category) & !is.na(MCF10A_TGFb.category)],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F) + 
  ggtitle("MCF10A WT -> CA1A (Figure4/5)") +
  theme(legend.position = "none")

ggsave(FigS5a,
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 5/FigS5a_parallelPlots.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

ggsave(FigS5b,
       filename = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 5/FigS5b_parallelPlots.pdf", 
       device = "pdf",
       height = 300,
       width = 400,
       units = "mm")

write.csv(table("WT" = dt5$MCF10A_WT.category, "TGFb" = dt5$MCF10A_TGFb.category), 
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 5/Fig5_WT_TGFb_counts.csv")
write.csv(table("WT" = dt5$MCF10A_WT.category, "CA1A" = dt5$MCF10CA1A.category),
          file = "~/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 5/Fig5_WT_CA1A_counts.csv")

