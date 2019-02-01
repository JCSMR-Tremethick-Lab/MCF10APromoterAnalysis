# MEME pipeline post-processing
library(data.table)
library(XML)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("Homo.sapiens")
library(gtools)

# load reference data -----------------------------------------------------
Hsap <- Homo.sapiens
genome <- BSgenome.Hsapiens.UCSC.hg19
TxDbUCSC <- TxDb.Hsapiens.UCSC.hg19.knownGene
knownGenes <- genes(Hsap, columns = "SYMBOL")
knownExons <- exonsBy(TxDbUCSC, by = "gene")


# basic routine for loading and extracting MEME CE XML data ------------------
memeResultsFile <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/ce/TOTALcombined_shH2AZ_Inp_000-125/meme.xml"
memeResultsCE <- XML::xmlToList(memeResultsFile)
l1 <- lapply(memeResultsCE$motifs, function(x) {
  r <- data.frame(lapply(x$.attrs, type.convert), stringsAsFactors=FALSE) # converts vector to DF
  return(r)
})
l2 <- lapply(memeResultsCE$training_set[which(names(memeResultsCE$training_set) == "sequence")], function(x) {
  r <- data.frame(lapply(x, type.convert), stringsAsFactors=FALSE) # converts vector to DF
  return(r)
})
l3 <- lapply(memeResultsCE$motifs, function(x) {
  motifName <- x$.attrs["name"]
  y <- x$contributing_sites
  l4 <- lapply(y, function(z){
    r <- as.data.table(lapply(z$`.attrs`, type.convert), stringsAsFactors = F)
    return(r)
  })
  tab1 <- do.call(rbind, lapply(l4, as.data.table))
  tab1$motifName <- motifName
  return(tab1)
})

memeResultsCETab <- do.call(rbind, lapply(l1, as.data.table))
setkey(memeResultsCETab, name)
memeResultsCESequences <- do.call(rbind, lapply(l2, as.data.table))
setkey(memeResultsCESequences, "name")
memeResultsCESites <- do.call(rbind, lapply(l3, as.data.table))
setkey(memeResultsCESites, "sequence_id")
rm(l1, l2, l3)

# basic routine for loading and extracting MEME CD XML data ------------------
memeResultsFile <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/meme/cd/TOTALcombined_shH2AZ_Inp_000-125/meme.xml"
memeResultsCD <- XML::xmlToList(memeResultsFile)
l1 <- lapply(memeResultsCD$motifs, function(x) {
  r <- data.frame(lapply(x$.attrs, type.convert), stringsAsFactors=FALSE) # converts vector to DF
  return(r)
})
l2 <- lapply(memeResultsCD$training_set[which(names(memeResultsCD$training_set) == "sequence")], function(x) {
  r <- data.frame(lapply(x, type.convert), stringsAsFactors=FALSE) # converts vector to DF
  return(r)
})
l3 <- lapply(memeResultsCD$motifs, function(x) {
  motifName <- x$.attrs["name"]
  y <- x$contributing_sites
  l4 <- lapply(y, function(z){
    r <- as.data.table(lapply(z$`.attrs`, type.convert), stringsAsFactors = F)
    return(r)
  })
  tab1 <- do.call(rbind, lapply(l4, as.data.table))
  tab1$motifName <- motifName
  return(tab1)
})

memeResultsCDTab <- do.call(rbind, lapply(l1, as.data.table))
setkey(memeResultsCDTab, name)
memeResultsCDSequences <- do.call(rbind, lapply(l2, as.data.table))
setkey(memeResultsCDSequences, "name")
memeResultsCDSites <- do.call(rbind, lapply(l3, as.data.table))
setkey(memeResultsCDSites, "sequence_id")
rm(l1, l2, l3)

memeResults <- rbind(memeResultsCETab,
                     memeResultsCDTab[!intersect(memeResultsCDTab$name, memeResultsCETab$name)])
# fixing column type in place
memeResults[, name := as.character(name)]
memeResults[, id := as.character(id)]
memeResults[, alt := as.character(alt)]

memeResultsCDTab[, name := as.character(name)]
memeResultsCDTab[, id := as.character(id)]
memeResultsCDTab[, alt := as.character(alt)]
memeResultsCDTab <- memeResultsCDTab[mixedorder(alt, decreasing = T)]

# load FIMO results -------------------------------------------------------
fimoResultsFile <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/fimo/ce/TOTALcombined_shH2AZ_Inp_000-125_first_run/fimo.tsv"
fimoResultsCE <- data.table::fread(fimoResultsFile)
fimoResultsFile <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/fimo/cd/TOTALcombined_shH2AZ_Inp_000-125_first_run/fimo.tsv"
fimoResultsCD <- data.table::fread(fimoResultsFile)
setkey(fimoResultsCE, motif_alt_id)
setkey(fimoResultsCD, motif_alt_id)




# load TOMTOM results -----------------------------------------------------
tomtomResultsFile <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/tomtom/ce/TOTALcombined_shH2AZ_Inp_000-125/tomtom.tsv"
tomtomResults <- data.table::fread(tomtomResultsFile)
setkey(tomtomResults, Target_ID)

# load hocomoco annotation ------------------------------------------------
hocomocoAnnotationFile <- "http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv"
hocomocoAnnotation <- data.table::fread(hocomocoAnnotationFile)
setkey(hocomocoAnnotation, Model)


# load ChIP-Seq data ------------------------------------------------------
topDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/"
clusterDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 9/"

peaks <- lapply(list.files(topDir, recursive = T, full.names = T, pattern = "125_peaks.xls"), function(x){
  return(data.table::fread(x))
})
names(peaks) <- unlist(lapply(strsplit(list.files(topDir, recursive = T, full.names = F, pattern = "125_peaks.xls"), "/"), function(x) x[1]))
peaks <- lapply(peaks, function(x) {
  x$peak_alt_name <- unlist(lapply(strsplit((x$name), "_"), function(x) paste(x[c(2,3,5,6)], collapse = "_")))
  return(x)
})
lapply(peaks, function(x) setkey(x, peak_alt_name))
key(peaks$`TOTALcombined_shH2AZ_Inp_000-125`)

# load cluster information from Figure 9 ----------------------------------
Fig9Clusters <- lapply(list.files(clusterDir, full.names = T), function (x){
  data.table::fread(x)
})
names(Fig9Clusters) <- unlist(lapply(strsplit(list.files(clusterDir, full.names = F), "_"), function(x) paste(x[1:4], collapse = "_")))

Fig9Clusters <- lapply(Fig9Clusters, function(x) {
  ucscID <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[4]))
  extGene <- unlist(lapply(strsplit(x$gene, ";"), function(x) x[6]))
  x$ucscID <- ucscID
  x$extGene <- extGene
  x$chr <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[1]))
  x$start <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[2]))
  x$end <-  unlist(lapply(strsplit(x$gene, ";"), function(x) x[3]))
  return(x)
})
Fig9Clusters

grlFig9Clusters <- lapply(Fig9Clusters, function(x){
  gr <- GenomicRanges::makeGRangesFromDataFrame(x, 
                                                seqinfo = seqinfo(genome),
                                                keep.extra.columns = T,
                                                starts.in.df.are.0based = T)
  return(gr)
})
grlFig9Clusters <- GRangesList(grlFig9Clusters)
names(grlFig9Clusters)
# Figure 9A is most relevant
smallFragments <- grlFig9Clusters$Fig9A_SFs_shZ_max
table(smallFragments$group1, smallFragments$color)
plot(c(1,2,3,4), rep(1, 4), type = "p", col = c("#FF0000FF","#80FF00FF", "#00FFFFFF", "#8000FFFF"))

grPromoters <- GenomicRanges::makeGRangesFromDataFrame(Fig9Clusters$Fig9A_SFs_10A_k4, 
                                                       seqinfo = seqinfo(genome),
                                                       keep.extra.columns = T, 
                                                       starts.in.df.are.0based = T)
grPromoters <- grPromoters[, c("extGene", "ucscID")]
names(grPromoters) <- grPromoters$ucscID


# load sleuth processed expression data -----------------------------------
load("~/Data/Tremethick/Breast/RNA-Seq/combined/R_Analysis/sleuthResults_GRCh37_hg19_UCSC_V1.rda")
kT <- results$kallisto_table_genes
rT <- results$sleuth_results_genes$conditionMCF10AshZD8
rm(results)
setkey(rT, target_id)
kT <- merge(kT, rT[,c("gene_id", "target_id")], all.x = T, all.y = F)
setkey(kT, gene_id)

kT1 <- kT[grep("MCF10A_wt|MCF10AshZD8", condition)]

# diagnostic plots (expression) -------------------------------------------
plotData <- kT1[hocomocoAnnotation[tomtomResults[Query_ID == "YSATTGGCT"]$Target_ID]$EntrezGene]
plotData <- plotData[!is.na(target_id)]
MOTIF1 <- ggplot(data = plotData, 
              aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
MOTIF1

# integrative analysis ----------------------------------------------------
# meme -> promoters
grPeakshH2AZ <- GenomicRanges::makeGRangesFromDataFrame(peaks$`TOTALcombined_shH2AZ_Inp_000-125`,
                                                        keep.extra.columns = T)

names(grPeakshH2AZ) <- grPeakshH2AZ$peak_alt_name
hocomocoAnnotation[tomtomResults[Query_ID %in% memeResults$name][,.(unique(Target_ID))]]

# looking at MEME-1 -------------------------------------------------------
m <- "MEME-1"
meme1Targets <- IRanges::subsetByOverlaps(grPromoters,
                                          grPeakshH2AZ[fimoResults[memeResults[m]][start > 230 & start < 300]$sequence_name])$extGene
meme1TFs <- hocomocoAnnotation[tomtomResults[Query_ID == unique(fimoResults[memeResults[m]]$motif_id)]$Target_ID]$EntrezGene
rT[gene_id %in% meme1TFs][b > 0]
table(rT[meme1Targets]$b > 0)
table(rT[gene_id %in% meme1TFs]$b >0)

# differential expression MEME-8 target genes-----------------------------
vp1 <- ggplot(data = rT, 
              aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)
vp1 <- vp1 + geom_point(data = rT[meme1Targets], color = "black")
vp1 <- vp1 + geom_point(data = rT[gene_id %in% meme1TFs], color = "lightblue")
vp1

# looking at MEME-8 -------------------------------------------------------
m <- "MEME-8"
meme8Targets <- IRanges::subsetByOverlaps(grPromoters,
                                          grPeakshH2AZ[fimoResults[memeResults[m]][start > 230 & start < 300]$sequence_name])$extGene
meme8TFs <- hocomocoAnnotation[tomtomResults[Query_ID == unique(fimoResults[memeResults[m]]$motif_id)]$Target_ID]$EntrezGene
rT[gene_id %in% meme8TFs][b > 0]
table(rT[meme8Targets]$b > 0)

# differential expression MEME-8 target genes-----------------------------
vp8 <- ggplot(data = rT, 
              aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)
vp8 <- vp8 + geom_point(data = rT[meme8Targets], color = "black")
vp8 <- vp8 + geom_point(data = rT[gene_id %in% meme8TFs], color = "lightblue")
vp8





