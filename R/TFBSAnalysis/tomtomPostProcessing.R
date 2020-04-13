library(data.table)
tomtom <- data.table::fread("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/tomtom/cd/TOTALcombined_shH2AZ_Inp_000-125/tomtom.tsv")
tomtom <- rbind(tomtom, data.table::fread("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments/tomtom/ce/TOTALcombined_shH2AZ_Inp_000-125/tomtom.tsv"))
load("/home/sebastian/Data/Tremethick/Breast/RNA-Seq/combined/R_Analysis/sleuthResults_GRCh37_hg19_UCSC_V1.rda")
rTshH2AZ <- results$sleuth_results_genes$conditionMCF10A_shZ
setkey(rTshH2AZ, target_id)

summary(tomtom$`q-value`)
summary(tomtom$`E-value`)

tomtom[`E-value` < 0.1][`q-value` < 0.1]
length(unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$Target_ID))
tomtom$extGene <- unlist(lapply(strsplit(tomtom$Target_ID, "_"), function(x) x[1]))
setkey(tomtom, extGene)

rTshH2AZ[which(rTshH2AZ$target_id %in% unique(tomtom[`E-value` < 0.1]$extGene))]
rTshH2AZ[which(rTshH2AZ$target_id %in% unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$extGene))]
table(rTshH2AZ[which(rTshH2AZ$target_id %in% unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$extGene))]$b > 0)

# 1. find query IDs which had hits in tomtom search & are expressed --------
extGenes <- rTshH2AZ[unique(tomtom[`E-value` < 0.1]$extGene)][!is.na(pval)]$target_id 
# returns 103 unique genes

extGenesUp <- rTshH2AZ[unique(tomtom[`E-value` < 0.1]$extGene)][!is.na(pval)][b > 0]$target_id 
# 61 up regulated

extGenesDown <- rTshH2AZ[unique(tomtom[`E-value` < 0.1]$extGene)][!is.na(pval)][b < 0]$target_id 
# 42 down regulated

# a substantial number of query motifs are redundant between the up/down regulated TFs
# differential expression analysis thus only makes sense to run on the unique sets
setkey(tomtom, "extGene")
queryIDs <- unique(tomtom[extGenes][`E-value` < 0.1]$Query_ID)
queryIDsU <- unique(tomtom[extGenesUp][`E-value` < 0.1]$Query_ID)
queryIDsD <- unique(tomtom[extGenesDown][`E-value` < 0.1]$Query_ID)

# mutually exclusiv query motifs
queryIDsUp <- queryIDsU[-which(queryIDsU %in% intersect(queryIDsU, queryIDsD))]
queryIDsDown <- queryIDsD[-which(queryIDsD %in% intersect(queryIDsU, queryIDsD))]

# 2. get sequencIDs of those sites ----------------------------------------
setkey(memeResultsCDSites, "motifName")
setkey(memeResultsCESites, "motifName")

sequenceIDs <- as.character(memeResultsCDSites[queryIDs]$sequence_id)
sequenceIDs <- c(sequenceIDs, 
                 as.character(memeResultsCESites[intersect(memeResultsCESites$motifName, queryIDs)]$sequence_id))
length(unique(sequenceIDs))
length((sequenceIDs))

sequenceIDsUp <- as.character(memeResultsCDSites[queryIDsUp]$sequence_id)
sequenceIDsUp <- c(sequenceIDsUp, 
                 as.character(memeResultsCESites[intersect(memeResultsCESites$motifName, queryIDsUp)]$sequence_id))
length(unique(sequenceIDsUp))
length((sequenceIDsUp))

sequenceIDsDown <- as.character(memeResultsCDSites[queryIDsDown]$sequence_id)
sequenceIDsDown <- c(sequenceIDsDown, 
                   as.character(memeResultsCESites[intersect(memeResultsCESites$motifName, queryIDsDown)]$sequence_id))
length(unique(sequenceIDsDown))
length((sequenceIDsDown))

length(intersect(unique(sequenceIDsDown), 
          unique(sequenceIDsUp)))
# 3. get peakIDs ----------------------------------------------------------
setkey(memeResultsCDSequences, "id")
setkey(memeResultsCESequences, "id")
peakIDs <- as.character(memeResultsCDSequences[unique(sequenceIDs)]$name)


peakIDs <- peakIDs[!is.na(peakIDs)]

# 4. now we can get the actual locations of these peaks -------------------
smallFragmentsSummits[peakIDs]

# 5. and cross them with H2A.Z peaks --------------------------------------
H2AZHits <- findOverlaps(smallFragmentsSummits[peakIDs], grNPIDRFiltered, minoverlap = 1)

# 6. and check in which gene promoters they are located -------------------
expressedTFPromoters <- findOverlaps(grNPIDRFiltered[subjectHits(H2AZHits)], grPromoters)
grPromoters[subjectHits(expressedTFPromoters)]


# todo --------------------------------------------------------------------
# determine how many TFBS in small fragment sites overlapping H2A.Z peaks are recognized by EXPRESSED TFs
setkey(tomtom, "extGene")

expressedTFBSmotifs <- unique(tomtom[extGenes]$Query_ID)
memeResultsSites <- rbind(memeResultsCDSites, memeResultsCESites)
memeResultsSites[, sequence_id := as.character(sequence_id)]
setkey(memeResultsSites, motifName)

memeResultsSequences <- rbind(memeResultsCDSequences, memeResultsCESequences)
memeResultsSequences[, id := as.character(id)]
memeResultsSequences[, name := as.character(name)]
setkey(memeResultsSequences, id)

memeResultsSequences[unique(memeResultsSites[expressedTFBSmotifs]$sequence_id)]
TFBSpeakIDs <- unique(memeResultsSequences[unique(memeResultsSites[expressedTFBSmotifs]$sequence_id)]$name)
smallFragmentsSummits[TFBSpeakIDs]

TFBSH2AZHits <- findOverlaps(smallFragmentsSummits[peakIDs], grNPIDRFiltered, minoverlap = 100)



# other stuff -------------------------------------------------------------
expressedTFTargets <- ggplot(data = rTshH2AZ, 
                          aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", size = 0.6) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")
expressedTFTargets + geom_point(data = rTshH2AZ[grPromoters[subjectHits(expressedTFPromoters)]$extGene][!is.na(pval)], color = "blue", size = 0.9)

table(rTshH2AZ[grPromoters[subjectHits(expressedTFPromoters)]$extGene][qval < 0.1]$b > 0)
kTGenes[grPromoters[subjectHits(H2AZHits)]$extGene]




vp2  <- ggplot(data = rTshH2AZ[rTshH2AZ[which(rTshH2AZ$target_id %in% tomtom[`E-value` < 0.1]$extGene)]], 
               aes(x = b, y = -(log10(qval)))) +
  xlim(c(-6.5,6.5)) +
  labs(x = "log2 fold-change") +
  geom_point(aes(color = target_id), size = 3) +
  geom_hline(yintercept = 2, linetype = "dotted") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")
vp2

table(rTshH2AZ[which(rTshH2AZ$target_id %in% unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$extGene))][qval < 0.1]$b > 0)
table(rTshH2AZ[which(rTshH2AZ$target_id %in% unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$extGene))]$qval < 0.1)
length(unique(expressedTFs))

# dot plots of single genes -----------------------------------------------
motif1 <- ggplot(data = kTGenes[tomtom[Query_ID == "YSATTGGC"]$extGene][grep("MCF10AshZD8|MCF10AD8", sample)], 
                     aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
motif1

motif2 <- ggplot(data = kTGenes[tomtom[Query_ID == "GCCAATCA"]$extGene][grep("MCF10AshZD8|MCF10AD8", sample)], 
                 aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
motif2

motif3 <- ggplot(data = kTGenes[tomtom[Query_ID == "GGAGGCGG"]$extGene][grep("MCF10AshZD8|MCF10AD8", sample)], 
                 aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
motif3

motif4 <- ggplot(data = kTGenes[tomtom[Query_ID == "TCCCAGCA"]$extGene][grep("MCF10AshZD8|MCF10AD8", sample)], 
                 aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
motif4

motif5 <- ggplot(data = kTGenes[tomtom[Query_ID == "TGASTCAB"]$extGene][grep("MCF10AshZD8|MCF10AD8", sample)], 
                 aes(x = condition, y = log2(tpm + 1))) + 
  geom_jitter(aes(color = sample), width = 0.05) + 
  facet_wrap(~target_id) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, size = 14))
motif5

# against background of all genes
vp1 <- ggplot(data = rTshH2AZ[qval < 0.1], 
              aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey") +
  #  geom_point(aes(color = target_id)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")
vp1 + geom_point(data = rTshH2AZ[rTshH2AZ[which(rTshH2AZ$target_id %in% tomtom[`E-value` < 0.1]$extGene)]], aes(color = target_id), size = 4)


# select TFs which are expressed ------------------------------------------
expressedTFs <- intersect(unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$extGene), rTshH2AZ$target_id)
length(expressedTFs)

expressedTFsUp <- intersect(unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$extGene), rTshH2AZ[b > 0 & qval < 0.1]$target_id)
expressedTFsDown <- intersect(unique(tomtom[`E-value` < 0.1][`q-value` < 0.1]$extGene), rTshH2AZ[b < 0 & qval < 0.1]$target_id)

vp1 <- ggplot(data = rTshH2AZ[qval < 0.1], 
              aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", size = 0.76) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")
vp1 <- vp1 + geom_point(data = rTshH2AZ[expressedTFsDown], color = "red", size = 1.4) 
vp1 + geom_point(data = rTshH2AZ[expressedTFsUp], color = "green", size = 1.4)

setkey(tomtom, extGene)
expTFsMotifs <- unique(tomtom[expressedTFs]$Query_ID)
expTFsUpMotifs <- unique(tomtom[expressedTFsUp]$Query_ID)
expTFsDownMotifs <- unique(tomtom[expressedTFsDown]$Query_ID)
length(expTFsMotifs)

fimoResults <- rbind(fimoResultsCD, fimoResultsCE)
setkey(fimoResults, motif_id)
fimoResultsSum <- fimoResults[`q-value` < 0.1][ , length(sequence_name), by = list(motif_alt_id, motif_id)]
fimoResultsSum[order(V1, decreasing = T)]
fimoResults <- fimoResults[`q-value` < 0.1]
fimoResults[expTFsMotifs][ , length(sequence_name), by = list(motif_alt_id, motif_id)][order(V1, decreasing = T)]
fimoResults[expTFsUpMotifs][ , length(sequence_name), by = list(motif_alt_id, motif_id)][order(motif_id)][!is.na(motif_alt_id)]
fimoResults[expTFsDownMotifs][ , length(sequence_name), by = list(motif_alt_id, motif_id)][order(motif_id)][!is.na(motif_alt_id)]

fimoResults[expTFsMotifs][!is.na(motif_alt_id)]

i1 <- intersect(names(smallFragmentsSummits), fimoResults[expTFsMotifs][!is.na(motif_alt_id)]$sequence_name)
upSummits <- unique(intersect(names(smallFragmentsSummits), fimoResults[expTFsUpMotifs][!is.na(motif_alt_id)]$sequence_name))
downSummits <- intersect(names(smallFragmentsSummits), fimoResults[expTFsDownMotifs][!is.na(motif_alt_id)]$sequence_name)
length(i1)

fO1 <- findOverlaps(smallFragmentsSummits[i1], grNPIDR[which(grNPIDR$score > 415)])

sFSUp <- smallFragmentsSummits[i1][queryHits(fO1)]
sFSUp[names(sFSUp) %in% upSummits]

grNPIDR[which(grNPIDR$score > 415)][subjectHits(fO1)]




fOUp <- findOverlaps(smallFragmentsSummits[queryHits(fO1)][upSummits], grNPIDR[which(grNPIDR$score > 415)])
fODown <- findOverlaps(smallFragmentsSummits[downSummits], grNPIDR[which(grNPIDR$score > 415)])

grSfp[queryHits(fO1)]
grNPIDR[which(grNPIDR$score > 415)][subjectHits(fO1)]


h2azPromoterHits <- findOverlaps(knownGenesPromoters, grNPIDR[which(grNPIDR$score > 415)][subjectHits(fO1)])
knownGenesPromoters[queryHits(h2azPromoterHits)]$SYMBOL
knownGenesPromotersExpression <- rTshH2AZ[unlist(knownGenesPromoters[queryHits(h2azPromoterHits)]$SYMBOL)]
knownGenesPromotersExpression <- knownGenesPromotersExpression[!is.na(knownGenesPromotersExpression$pval)]

TFUpTargets <- findOverlaps(smallFragmentsSummits[upSummits][queryHits(fOUp)], knownGenesPromoters)
TFDownTargets <- findOverlaps(smallFragmentsSummits[downSummits], knownGenesPromoters)



vpKGP <- ggplot(data = rTshH2AZ[qval < 0.1], 
                aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", size = 0.8) +
  #  geom_point(aes(color = target_id)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")
vpKGP + geom_point(data = knownGenesPromotersExpression, aes(color = "red"), size = 1.4)


