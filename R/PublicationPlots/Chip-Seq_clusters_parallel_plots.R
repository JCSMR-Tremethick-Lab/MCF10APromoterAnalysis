## local functions
median.quartile <- function(x){
  out <- quantile(x, probs = c(0.25,0.5,0.75))
  names(out) <- c("ymin","y","ymax")
  return(out) 
}


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
gr.promoters <- GenomicRanges::makeGRangesFromDataFrame(l1$Total_10A, ignore.strand = T, keep.extra.columns = T)

rm(dt1)
dt1 <- merge(subset(l1$Total_10A, select = c("extGene", "group1")), subset(l1$Total_TGFb, select = c("extGene", "group1")), by.x = "extGene", by.y = "extGene")
colnames(dt1)[2:3] <- c("MCF10A_WT", "MCF10A_TGFb")
dt1 <- merge(dt1, subset(l1$Total_shZ, select = c("extGene", "group1")), by.x = "extGene", by.y = "extGene")
dt1 <- merge(dt1, subset(l1$Total_CA1a, select = c("extGene", "group1")), by.x = "extGene", by.y = "extGene")
colnames(dt1)[4:5] <- c("MCF10A_shZ", "MCF10CA1A")
colnames(dt1)


# For input (figure 1)
# We have decided on the following classes 
# 1. Strong +2 (1 for 10A, 1 for TGFb, 1 for shH2AZ, 1 for CAC1)
# 2. Strong +1 (2 for 10A, 2 for TGFb, 2 for shH2AZ, 2 for CAC1)
# 3. weak upstream (3 for 10A, 5 for TGFb, 5 for shH2AZ, 5 for CAC1)
# 4. weak  array across promoter (4 for 10A, 7 for TGFb, 7 for shH2AZ, 7 for CAC1)
# 5. Strong -2 (5 for 10A, 6 for TGFb, 6 for shH2AZ, 6 for CAC1)
# 6. Strong -1 (6 for 10A, 4 for TGFb, 4 for shH2AZ, 4 for CAC1)
# 7. Overall depleted  (7 for 10A, 3 for TGFb, 3 for shH2AZ, 3 for CAC1)

mcf10awtCategories <- data.table::data.table(cat = c("Strong +2", "Strong +1", "weak upstream", "weak array across promoter", "Strong -2", "Strong -1", "Overall depleted"),
                               group = c(1,2,3,4,5,6,7))

mcf10atgfbCategories <- data.table::data.table(cat = c("Strong +2", "Strong +1", "weak upstream", "weak array across promoter", "Strong -2", "Strong -1", "Overall depleted"),
                               group = c(1,2,5,7,6,4,3))

mcf10ashh2azCategories <- data.table::data.table(cat = c("Strong +2", "Strong +1", "weak upstream", "weak array across promoter", "Strong -2", "Strong -1", "Overall depleted"),
                                   group = c(1,2,5,7,6,4,3))

mcf10cac1Categories <- data.table::data.table(cat = c("Strong +2", "Strong +1", "weak upstream", "weak array across promoter", "Strong -2", "Strong -1", "Overall depleted"),
                                   group = c(1,2,5,7,6,4,3))

dt1$MCF10A_WT.category <- as.factor(mcf10awtCategories$cat[match(dt1$MCF10A_WT, mcf10awtCategories$group)])
dt1$MCF10A_TGFb.category <- as.factor(mcf10atgfbCategories$cat[match(dt1$MCF10A_TGFb, mcf10atgfbCategories$group)])
dt1$MCF10A_shZ.category <- as.factor(mcf10ashh2azCategories$cat[match(dt1$MCF10A_shZ, mcf10ashh2azCategories$group)])
dt1$MCF10CA1A.category <- as.factor(mcf10cac1Categories$cat[match(dt1$MCF10CA1A, mcf10cac1Categories$group)])

ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_TGFb.category"), #, "MCF10A_shZ.category", "MCF10CA1A.category"), 
                       data = dt1,
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F)

ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_shZ.category"), #, , "MCF10CA1A.category"), 
                       data = dt1,
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F)

ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10CA1A.category"), #, , ), 
                       data = dt1,
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F)

factor(dt1$MCF10A_TGFb.category)

table(dt1$MCF10A_WT.category) 
table(dt1$MCF10A_TGFb.category)
table("WT" = dt1$MCF10A_WT.category, "TGFb" = dt1$MCF10A_TGFb.category)
table("WT" = dt1$MCF10A_WT.category, "shZ" = dt1$MCF10A_shZ.category)
table("WT" = dt1$MCF10A_WT.category, "WT" = dt1$MCF10A_WT.category)
table("TGFb" = dt1$MCF10A_TGFb.category, "TGFb" = dt1$MCF10A_TGFb.category)

levels(dt1$MCF10A_shZ.category)

## gather table
Fig2Categories <- dt1[, c("extGene", "MCF10A_WT.category", "MCF10A_TGFb.category", "MCF10A_shZ.category", "MCF10CA1A.category")]
Fig2Categories <- data.table(gather(Fig2Categories, key = "condition", value = "category", c("MCF10A_WT.category", "MCF10A_TGFb.category", "MCF10A_shZ.category", "MCF10CA1A.category")))
Fig2Categories$condition <- unlist(lapply(strsplit(Fig2Categories$condition, "\\."), function(x) x[1]))
Fig2Categories$category <- factor(Fig2Categories$category, levels = levels(as.factor(Fig2Categories$category))[c(5,4,7,6,3,2,1)])

deGenesTGFbIDs <- deGenesTGFb[deGenesTGFb$qval < 0.1]$target_id
deGenesTGFbUpIDs <- deGenesTGFb[deGenesTGFb$qval < 0.1 & deGenesTGFb$b > 0]$target_id
deGenesTGFbDownIDs <- deGenesTGFb[deGenesTGFb$qval < 0.1 & deGenesTGFb$b < 0]$target_id

left <- "MCF10A_WT"
right <- "MCF10A_TGFb"
left.categories <- llist
left.categories <- c("Overall depleted", "weak array across promoter")
df1 <- data.table(Fig2Categories[condition == left]$extGene,  Fig2Categories[condition == left]$category, Fig2Categories[condition == right]$category)
colnames(df1) <- c("extGene", left, right)
#df1 <- df1[extGene %in% deGenesTGFbUpIDs]
#df1 <- df1[extGene %in% deGenesTGFbDownIDs]
df1 <- df1[extGene %in% sigEMTCells$cellLine_sig]

df1 <- df2["epi"]
data <- data.frame(df1[ df1[[left]] %in% left.categories , ])
ggparallel::ggparallel(vars = list(left, right), data = data, same.level = T, label.size = 3,
                       text.angle = 0, label.colour= "black") + ggtitle("Epithelial")

geneIDs <- data$extGene
s1 <- which(kt1D6$target_id %in% geneIDs)
s2 <- which(deGenesTGFb$target_id %in% sigEMTCells$cellLine_sig)
s2
plot(deGenesTGFb[s2]$b, -log10(deGenesTGFb[s2]$qval))

bp1 <- ggplot(data = kt1D6[s1],
              aes(x= condition, y = log2(tpm + 1), colour = condition)) + geom_boxplot()
bp1
bp2 <- ggplot(data =m3[m3[geneIDs][condition == "MCF10A_wt"][category %in% "Strong +1"]$target_id],
              aes(x= condition, y = log2(tpm + 1), colour = condition)) + geom_boxplot() + facet_wrap(~category)
  
bp2
vp1 <- ggplot(data = kt1D6[s1],
              aes(x= condition, y = log2(tpm + 1)), fill = condition) + geom_violin()
vp1
h1 <- ggplot(data = kt1D6[s1],
             aes(log2(tpm + 1), ..density.., fill = condition)) + geom_histogram(binwidth = 0.11, alpha = 0.5)
h1

###
data <- as.data.frame(df1)
vars <- unlist(list(left, right))
llist <- NULL
same.level = T
for (i in unique(vars)) {
  if (!same.level) levels(data[,i]) <- paste(i, levels(data[,i]), sep=":")
  llist <- unique(c(llist, levels(data[,i])))
}
(llist)
## integrate expression data
kt1 <- results$full$kallisto_table_genes
#kt1 <- results$full$kallisto_table_genes
setkey(kt1, "target_id")

kt1D6 <- kt1[grep("D6", kt1$sample)]
kt1D6.wide <- tidyr::spread(kt1D6[, c("target_id", "sample", "tpm")], sample, tpm)
kt1D6 <- kt1D6[, list("tpm" = mean(tpm)), by = list(condition, target_id)]
s1 <- dt1[(dt1$MCF10A_WT.category == "Strong +1" & dt1$MCF10A_TGFb.category == "Strong -1")]$extGene
s2 <- which(kt1D6.wide$target_id %in% s1)
s1 <- which(kt1D6$target_id %in% s1)
length(s1)

gplots::heatmap.2(scale(log2(as.matrix(subset(kt1D6.wide[s2], select = c("MCF10AD6_1",
                                                              "MCF10AD6_2",
                                                              #"MCF10AD6_3",
                                                              "MCF10ATGFbD6_1",
                                                              "MCF10ATGFbD6_2",
                                                              "MCF10ATGFbD6_3"))) + 1)), trace = "none")
m1 <- scale(log2(as.matrix(subset(kt1D6.wide[s2], select = c("MCF10AD6_1",
                                                             "MCF10AD6_2",
                                                             #"MCF10AD6_3",
                                                             "MCF10ATGFbD6_1",
                                                             "MCF10ATGFbD6_2",
                                                             "MCF10ATGFbD6_3"))) + 1))
m2 <- log2(as.matrix(subset(kt1D6.wide[s2], select = c("MCF10AD6_1",
                                                       "MCF10AD6_2",
                                                       #"MCF10AD6_3",
                                                       "MCF10ATGFbD6_1",
                                                       "MCF10ATGFbD6_2",
                                                       "MCF10ATGFbD6_3"))) + 1)

rownames(m1) <- kt1D6.wide[s2]$target_id
rownames(m2) <- kt1D6.wide[s2]$target_id


head(m1)
heatmaply::heatmaply(m1,#[sd1 > 0.3,], 
                     trace = "none",
                     scale = "none",
                     protocol = "ggplotly")
library(vegan)
km1 <- vegan::cascadeKM(m2, 2, 11, criterion = "calinski")
plot(km1)
data.table(km1$partition)
g1 <- km1$partition[,5]
names(g1)  <- rownames(m2)
g1
setkey(kt1D6, "target_id")
w1 <-which(kt1D6$target_id %in% names(g1))
kt1D6.g1 <- kt1D6[w1]
kmeansClusters <- data.table("extGene" = names(g1), "group" = g1)
kt1D6.g1 <- merge(kt1D6.g1, kmeansClusters, by.x = "target_id", by.y = "extGene", all.y = T)

bp2 <- ggplot(data = kt1D6.g1,
              aes(x= condition, y = log2(tpm + 1), colour = condition)) +
  geom_boxplot() + 
  facet_wrap(~group)
bp2

vp2 <- ggplot(data = kt1D6.g1,
              aes(x= condition, y = log2(tpm + 1), colour = condition)) +
  geom_violin() + 
  facet_wrap(~group) +
  stat_summary(fun.y=median.quartile,geom='point')
vp2

h2 <- ggplot(data = kt1D6.g1,
             aes(log2(tpm + 1), fill = condition)) + geom_histogram(bins = 30, alpha = 0.5,)
h2

f2 <- ggplot(data = kt1D6.g1,
             aes(log2(tpm + 1), colour = condition)) + geom_freqpoly(bins = 30, alpha = 0.5)
f2

median(log2(kt1D6[s1][condition == "MCF10A_wt"]$tpm + 1))
median(log2(kt1D6[s1][condition == "MCF10A_TGFb"]$tpm + 1))

median(log2(kt1D6[condition == "MCF10A_wt"]$tpm + 1))
median(log2(kt1D6[condition == "MCF10A_TGFb"]$tpm + 1))


# EMT genes-specific parallel plots ---------------------------------------
ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_shZ.category"), #, , "MCF10CA1A.category"), 
                       data = dt1[unique(ucscTranscriptsSigEMTCells[epi_mes == "epi"]$target_id)],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F) + ggtitle("shZ - EPI")
"ParallelPlot_WT-shZ_EPI_Genes.pdf"


ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_shZ.category"), #, , "MCF10CA1A.category"), 
                       data = dt1[unique(ucscTranscriptsSigEMTCells[epi_mes == "mes"]$target_id)],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F) + ggtitle("shZ - MES")
"ParallelPlot_WT-shZ_MES_Genes.pdf"

#----- TGFb
ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_TGFb.category"), #, , "MCF10CA1A.category"), 
                       data = dt1[unique(ucscTranscriptsSigEMTCells$target_id)],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F)

ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_TGFb.category"), #, , "MCF10CA1A.category"), 
                       data = dt1[unique(ucscTranscriptsSigEMTCells[epi_mes == "epi"]$target_id)],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F) + ggtitle("TGFb - EPI")
"ParallelPlot_WT-TGFb_EPI_Genes.pdf"

ggparallel::ggparallel(list("MCF10A_WT.category", "MCF10A_TGFb.category"), #, , "MCF10CA1A.category"), 
                       data = dt1[unique(ucscTranscriptsSigEMTCells[epi_mes == "mes"]$target_id)],
                       label.size = 4,
                       text.angle = 0, label.colour= "black", same.level = F) + ggtitle("TGFb - MES")
"ParallelPlot_WT-TGFb_MES_Genes.pdf"

