# run after peakCallingParameterOptimisation.R

library(GenomicRanges)

# load smallFragments -----------------------------------------------------
smallFragmentsDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/macs2PeakCalling"
smallFragmentPeaks <- list.files(smallFragmentsDir, recursive = T, full.names = T, pattern = "peaks.xls")
smallFragmentPeaks <- lapply(smallFragmentPeaks, function(x) {dT <- data.table::fread(x); return(dT)})
names(smallFragmentPeaks) <- unlist(lapply(strsplit(list.files(smallFragmentsDir, recursive = T, full.names = F, pattern = "peaks.xls"), "/"), function(x) x[1]))
grlSallFragmentPeaks <- lapply(smallFragmentPeaks, function(x) gr <- makeGRangesFromDataFrame(x, keep.extra.columns = T))

sfp <- smallFragmentPeaks$`TOTALcombined_shH2AZ_Inp_000-125`[(fold_enrichment > 4 & `-log10(qvalue)` > 2)]
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
grNPIDRFiltered[subjectHits(Hits)]
summary(grNPIDRFiltered[subjectHits(Hits)]$score)
smallFragmentsSummits[queryHits(Hits)]
length(unique(subjectHits(Hits)))
length(unique(queryHits(Hits)))


# analysis of expression of small fragment promoter genes -----------------
# approach:
# 1. identify promoters with H2A.Z peaks that give rise to small fragments in shH2A.Z
# 2. record changes in gene expression
# 3. correlate actual peak positions within promoters of changing genes and small fragment peaks in shH2AZ

# load pre-processed expression data --------------------------------------
load("~/Data/Tremethick/Breast/RNA-Seq_run2/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/R_Analysis/sleuthResults_GRCh37_hg19_ensembl75_ERCC_V2.rda")
rTshH2AZ <- as.data.table(results$MCF10A_vs_shZ$sleuth_results_genes$conditionMCF10A_shZ)
setkey(rTshH2AZ, target_id)
rTshH2AZ <- rTshH2AZ[!is.na(pval)]
kTGenes <- results$MCF10A_vs_shZ$kallisto_table_genes
setkey(kTGenes, target_id)
rm(results)

# prepare promoter regions and cross with H2AZ peaks giving rise t --------
knownGenesPromoters <- promoters(knownGenes, upstream = 1000, downstream = 1000)
# V this filters for promoters with H2A.Z peaks associated with small fragments
promoterHits <- IRanges::findOverlaps(grNPIDR[subjectHits(Hits)], knownGenesPromoters) 
knownGenesPromoters[unique(subjectHits(promoterHits))]
knownGenesPromotersExpression <- rTshH2AZ[unlist(knownGenesPromoters[unique(subjectHits(promoterHits))]$SYMBOL)]
knownGenesPromotersExpression <- knownGenesPromotersExpression[!is.na(knownGenesPromotersExpression$pval)]

ph2 <- findOverlaps(grNPIDR, knownGenesPromoters)
knownGenesPromoters[unique(subjectHits(ph2))]
rTshH2AZ[unlist(knownGenesPromoters[unique(subjectHits(ph2))]$SYMBOL)]
i1 <- intersect(unlist(knownGenesPromoters[unique(subjectHits(ph2))]$SYMBOL), knownGenesPromotersExpression$target_id)


# looking at some summary statistics --------------------------------------
summary(kTGenes[condition == "MCF10A_wt"]$tpm)
summary(kTGenes[condition == "MCF10A_wt"][unlist(knownGenesPromoters[unique(subjectHits(ph2))]$SYMBOL)][!target_id %in% i1]$tpm)
summary(kTGenes[condition == "MCF10A_wt"][unlist(knownGenesPromoters[unique(subjectHits(promoterHits))]$SYMBOL)]$tpm)

summary(kTGenes[condition == "MCF10A_shZ"]$tpm)
summary(kTGenes[condition == "MCF10A_shZ"][unlist(knownGenesPromoters[unique(subjectHits(ph2))]$SYMBOL)][!target_id %in% i1]$tpm)
summary(kTGenes[condition == "MCF10A_shZ"][unlist(knownGenesPromoters[unique(subjectHits(promoterHits))]$SYMBOL)]$tpm)

summary(rTshH2AZ$var_obs)
summary(rTshH2AZ[unlist(knownGenesPromoters[unique(subjectHits(ph2))]$SYMBOL)][!target_id %in% i1]$var_obs)
summary(knownGenesPromotersExpression$var_obs)

summary(rTshH2AZ$b)
summary(rTshH2AZ[unlist(knownGenesPromoters[unique(subjectHits(ph2))]$SYMBOL)][!target_id %in% i1][!is.na(b)]$b)
summary(knownGenesPromotersExpression$b)

ggplot(data = rTshH2AZ, aes(x = b)) + geom_density() +
  geom_density(data = rTshH2AZ[unlist(knownGenesPromoters[unique(subjectHits(ph2))]$SYMBOL)][!target_id %in% i1][!is.na(b)], aes(x = b), colour = "green") +
  geom_density(data = knownGenesPromotersExpression, aes(x = b), colour = "red")

ggplot(data = rTshH2AZ, aes(x = qval)) + geom_density(colour = "grey", fill = "grey", alpha = 0.3) +
  geom_density(data = rTshH2AZ[unlist(knownGenesPromoters[unique(subjectHits(ph2))]$SYMBOL)][!target_id %in% i1][!is.na(b)], aes(x = qval), colour = "green", fill = "green", alpha = 0.3) +
  geom_density(data = knownGenesPromotersExpression, aes(x = qval), colour = "red", fill = "red", alpha = 0.3)

# volcano plot with small fragment genes highlighted ----------------------
library(ggrepel)
vpKGP <- ggplot(data = rTshH2AZ, 
              aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", size = 0.6) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")
vpKGP + geom_point(data = knownGenesPromotersExpression, color = "purple", size = 0.7) +
        ggrepel::geom_text_repel(data = knownGenesPromotersExpression[qval < 0.0001][abs(b) > 0.5], aes(label = target_id), color = "blue", size = 3)




