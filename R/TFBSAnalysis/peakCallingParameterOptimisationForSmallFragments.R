# inspect peak calling parameters
library("rtracklayer")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("BSgenome.Hsapiens.UCSC.hg19")
library("data.table")
library("GenomicRanges")


# load reference data -----------------------------------------------------
genome <- BSgenome.Hsapiens.UCSC.hg19
TxDbUCSC <- TxDb.Hsapiens.UCSC.hg19.knownGene
knownGenes <- genes(TxDbUCSC)
knownExons <- exonsBy(TxDbUCSC, by = "gene")


# set working and data directories ----------------------------------------
topDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/"
clusterDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 9/"

# load peaks first --------------------------------------------------------
peaks <- lapply(list.files(topDir, recursive = T, full.names = T, pattern = "125_peaks.xls"), function(x){
  return(data.table::fread(x))
})
names(peaks) <- unlist(lapply(strsplit(list.files(topDir, recursive = T, full.names = F, pattern = "125_peaks.xls"), "/"), function(x) x[1]))

lapply(peaks, function(x){
  summary(x$fold_enrichment)
})

lapply(peaks, function(x){
  summary(x$`-log10(qvalue)`)
})

lapply(peaks, function(x){
  table(x$`-log10(qvalue)` >= 2)
})

lapply(peaks, function(x){
  summary(x$length)
})

lapply(peaks, function(x){
  x[which.max(x$fold_enrichment)]
})


# load summits ------------------------------------------------------------
summits <- lapply(list.files(topDir, recursive = T, full.names = T, pattern = "125_summits.bed"), function(x){
  y <- rtracklayer::import(x)
  names(y) <- y$name
  return(y)
})
names(summits) <- unlist(lapply(strsplit(list.files(topDir, recursive = T, full.names = F, pattern = "125_summits.bed"), "/"), function(x) x[1]))

lapply(peaks, function(x){
  summary(x$fold_enrichment)
  nrow(x)
  table(x$fold_enrichment > 2)
  summary(x$pileup)
})


# extract DNA sequences around peaks --------------------------------------
summitSeqs <- lapply(summits, function(x){
  return(getSeq(genome, resize(x, 500, fix = "center")))
})
names(summitSeqs) <- names(summits)


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
                                                keep.extra.columns = T)
  return(gr)
})
grlFig9Clusters <- GRangesList(grlFig9Clusters)


grPromoters <- GenomicRanges::makeGRangesFromDataFrame(Fig9Clusters$Fig9A_SFs_10A_k4, 
                                                       seqinfo = seqinfo(genome),
                                                       keep.extra.columns = T, 
                                                       starts.in.df.are.0based = T)
grPromoters <- grPromoters[, c("extGene", "ucscID")]
names(grPromoters) <- grPromoters$ucscID

# exploratory analysis of the peak calling data to determine thres --------
peaks$`TOTALcombined_A_Inp_000-125`[(peaks$`TOTALcombined_A_Inp_000-125`$fold_enrichment > 3 & peaks$`TOTALcombined_A_Inp_000-125`$`-log10(qvalue)` > 14)]

# have to split data into ChIP and Input (MNase-Seq)
peaksMNase <- list(WT = peaks$`TOTALcombined_A_Inp_000-125`, 
                   TGFb = peaks$`TOTALcombined_A_TGFb_Inp_000-125`, 
                   Ca1a = peaks$`TOTALcombined_CA1a_Inp_000-125`,
                   shH2AZ = peaks$`TOTALcombined_shH2AZ_Inp_000-125`)
peaksH2AZ <- list(WT = peaks$`TOTALcombined_A_H2AZ_000-125`,
                  TGFb = peaks$`TOTALcombined_A_TGFb_H2AZ_000-125`,
                  shH2AZ = peaks$`TOTALcombined_CA1a_H2AZ_000-125`)

grlPeaksMNase <- lapply(peaksMNase, function(x){
  gr <- GenomicRanges::makeGRangesFromDataFrame(x, 
                                                seqinfo = seqinfo(genome),
                                                keep.extra.columns = T)
  return(gr)
})
grlPeaksMNase <- GenomicRanges::GRangesList(grlPeaksMNase)

lapply(names(peaksMNase), function(x){
  print(x)
  print(summary(peaksMNase[[x]]$fold_enrichment))
  print(nrow(peaksMNase[[x]]))
  print(table(peaksMNase[[x]]$fold_enrichment > 2))
  print(summary(peaksMNase[[x]]$`-log10(qvalue)`))
  print(table(peaksMNase[[x]]$`-log10(qvalue)` > 2))
  print(summary(peaksMNase[[x]]$pileup))
})

lapply(names(peaksH2AZ), function(x){
  print(x)
  print(summary(peaksH2AZ[[x]]$fold_enrichment))
  print(nrow(peaksH2AZ[[x]]))
  print(table(peaksH2AZ[[x]]$fold_enrichment > 2))
  print(summary(peaksH2AZ[[x]]$`-log10(qvalue)`))
  print(table(peaksH2AZ[[x]]$`-log10(qvalue)` > 2))
  print(summary(peaksH2AZ[[x]]$pileup))
})


grlSigPeaksMNase <- lapply(grlPeaksMNase, function(x){
  x <- x[which(x$pileup > 4 & x$`-log10(qvalue)` > 2)]
  x <- sort(x)
  return(x)
})
grlSigPeaksMNase <- GRangesList(grlSigPeaksMNase)
grlSigPeaksMNase <- subsetByOverlaps(grlSigPeaksMNase, grPromoters)

# testing with shH2AZ MNase-Seq data --------------------------------------
TOTALcombined_shH2AZ_Inp_PeakSeqs <- getSeq(genome,
                                            resize(summits$`TOTALcombined_shH2AZ_Inp_000-125`[grlSigPeaksMNase$shH2AZ$name],
                                                   500, 
                                                   fix = "center"))

names(TOTALcombined_shH2AZ_Inp_PeakSeqs) <- unlist(lapply(strsplit((names(TOTALcombined_shH2AZ_Inp_PeakSeqs)), "_"), function(x) paste(x[c(2,3,5,6)], collapse = "_")))

rtracklayer::export(TOTALcombined_shH2AZ_Inp_PeakSeqs, 
                    con = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling//TOTALcombined_shH2AZ_Inp_000-125/summitSequences.fa",
                    format = "fasta")
promoterBackground <- getSeq(genome, grPromoters)
rtracklayer::export(promoterBackground,
                    con = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling/promoterSequences.fa",
                    format = "fasta")


# use inverse of peakSeqs from promoter sequences as negative ctrl --------
f1 <- findOverlaps(reduce(grlSigPeaksMNase$shH2AZ), grPromoters)
seq1 <- getSeq(genome, grPromoters[to(f1)])

maskList <- lapply(1:length(seq1), function(x) {
  s1 <- seq1[[x]]
  start <- start(grlSigPeaksMNase$shH2AZ[from(f1[x])]) - start(grPromoters[to(f1[x])])
  start <- unlist(sapply(start, function(x) {if (x < 1){x <- 1}; return (x)}))
  end <- start + width(grlSigPeaksMNase$shH2AZ[from(f1[x])])
  m1 <- mask(s1, start = start, end = end)
  return(m1)
})

shH2AZNegativeControls1 <- getSeq(genome, grPromoters[-to(f1)])
names(shH2AZNegativeControls1)
rtracklayer::export(shH2AZNegativeControls1,
                    con = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling//TOTALcombined_shH2AZ_Inp_000-125/shH2AZNegativeControlPromoters.fa",
                    format = "fasta")

shH2AZNegativeControls2 <- DNAStringSetList(lapply(maskList, function(x) as(x, "DNAStringSet")))
rtracklayer::export(unlist(shH2AZNegativeControls2), 
                    con = "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling//TOTALcombined_shH2AZ_Inp_000-125/shH2AZNegativeControlPeakFreePromoterSequences.fa",
                    format = "fasta")
l
