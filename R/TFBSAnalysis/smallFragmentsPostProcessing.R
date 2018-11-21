# run after peakCallingParameterOptimisation.R

library(GenomicRanges)

# load smallFragments -----------------------------------------------------
smallFragmentsDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling"
smallFragmentPeaks <- list.files(smallFragmentsDir, recursive = T, full.names = T, pattern = "peaks.xls")
smallFragmentPeaks <- lapply(smallFragmentPeaks, function(x) {dT <- data.table::fread(x); return(dT)})
names(smallFragmentPeaks) <- unlist(lapply(strsplit(list.files(smallFragmentsDir, recursive = T, full.names = F, pattern = "peaks.xls"), "/"), function(x) x[1]))
grlSallFragmentPeaks <- lapply(smallFragmentPeaks, function(x) gr <- makeGRangesFromDataFrame(x, keep.extra.columns = T))

sfp <- smallFragmentPeaks$`TOTALcombined_shH2AZ_Inp_000-125`[(fold_enrichment > 3 & `-log10(qvalue)` > 2)]
grSfp <- makeGRangesFromDataFrame(sfp, keep.extra.columns = T)
names(grSfp) <- unlist(lapply(strsplit(grSfp$name, "_"), function(x) paste(x[c(2,3,5,6)], collapse = "_")))


# use results of IDR analysis ---------------------------------------------
grNPIDR <- data.table::fread("/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/GRCh37_hg19_UCSC/idr_noMacs2Cutoff/H2AZ_H_idr_values.txt")
colnames(grNPIDR) <- c("chrom", "chromStart", "chromEnd", "name", "score", 
                       "strand", "signalValue", "p-value", "q-value", "summit",
                       "localIDR", "globalIDR", 
                       "rep1_chromStart", "rep1_chromEnd", "rep1_signalValue", "rep1_summit",
                       "rep2_chromStart", "rep2_chromEnd", "rep2_signalValue", "rep2_summit")

grNPIDR <- makeGRangesFromDataFrame(grNPIDR, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd")
hist(grNPIDR$score)
summary(grNPIDR$score)


smallFragmentsSummits <- resize(summits$`TOTALcombined_shH2AZ_Inp_000-125`[grlSigPeaksMNase$shH2AZ$name],
       124, 
       fix = "center")
names(smallFragmentsSummits) <- gsub("TOTALcombined_", "", names(smallFragmentsSummits))
names(smallFragmentsSummits) <- gsub("000-125_", "", names(smallFragmentsSummits))

grNPIDRFiltered <- grNPIDR[which(grNPIDR$score > 415)]
Hits <- findOverlaps(smallFragmentsSummits, grNPIDRFiltered, minoverlap = 1)
grNPIDR[subjectHits(Hits)]
summary(grNPIDRFiltered[subjectHits(Hits)]$score)
smallFragmentsSummits[queryHits(Hits)]
length(unique(subjectHits(Hits)))
length(unique(queryHits(Hits)))

fimoResultsFiltered <- fimoResults[names(smallFragmentsSummits[queryHits(Hits)])]
fimoResultsFiltered <- fimoResultsFiltered[!is.na(fimoResultsFiltered$motif_alt_id)]

library(dplyr)
t1 <- fimoResultsFiltered %>% 
  group_by(motif_alt_id, motif_id) %>% 
  summarise(counts = length(sequence_name))
t1$motif_id

# approach:
# 1. identify promoters with H2A.Z peaks that give rise to small fragments in shH2A.Z
# 2. record changes in gene expression
# 3. correlate actual peak positions within promoters of changing genes and small fragment peaks in shH2AZ
knownGenesPromoters <- promoters(knownGenes, upstream = 1000, downstream = 1000)
promoterHits <- findOverlaps(4, knownGenesPromoters)
knownGenesPromoters[unique(subjectHits(promoterHits))]
knownGenesPromotersExpression <- rTshH2AZ[unlist(knownGenesPromoters[unique(subjectHits(promoterHits))]$SYMBOL)]
knownGenesPromotersExpression <- knownGenesPromotersExpression[!is.na(knownGenesPromotersExpression$pval)]

vpKGP <- ggplot(data = rTshH2AZ[qval < 0.1], 
              aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", size = 0.8) +
  #  geom_point(aes(color = target_id)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")
vpKGP + geom_point(data = rTshH2AZ[unlist(knownGenesPromoters[unique(subjectHits(promoterHits))]$SYMBOL)], aes(color = target_id), size = 2)
