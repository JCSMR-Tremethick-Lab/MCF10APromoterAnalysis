
# forward search ----------------------------------------------------------
# 1. peakIDs from smallFragmentsSummits X H2A.Z peaks
peakIDs <- names(smallFragmentsSummits[queryHits(Hits)])

# 2. get sequenceIDs - we only need to use one table as sequences are identical for CD/CE
setkey(memeResultsCDSequences, "name")
sequenceIDs <- as.character(memeResultsCDSequences[peakIDs]$id)

# 3. get motifNames
setkey(memeResultsCDSites, "sequence_id")
setkey(memeResultsCESites, "sequence_id")
memeResultsCESites
CDmotifNames <- as.character(memeResultsCDSites[sequenceIDs]$motifName)
CEmotifNames <- as.character(memeResultsCESites[sequenceIDs]$motifName)
CEmotifNames <- CEmotifNames[!is.na(CEmotifNames)]

length(unique(CDmotifNames))
length(unique(CEmotifNames))
intersect(CDmotifNames, CEmotifNames)

motifNames <- c(CDmotifNames, CEmotifNames[-which(CEmotifNames %in% intersect(CDmotifNames, CEmotifNames))])
motifTable <- table(motifNames)
motifName <- sort(motifTable, decreasing = T)
# 4. get ext_genes
setkey(tomtom, "Query_ID")
tomtom[names(motifNames)[6]]

# reverse search (from TF gene to sites) ----------------------------------
# 1. get query IDs that returned FOS sites
setkey(tomtom, "extGene")
fosQueries <- tomtom["FOS"]$Query_ID
junQueries <- tomtom["JUN"]$Query_ID
intersect(fosQueries, junQueries)
mycQueries <- tomtom["MYC"]$Query_ID

# 2. get sequencIDs of FOS/JUN sites
setkey(memeResultsCDSites, "motifName")
junFosSequenceIDs <- as.character(memeResultsCDSites[intersect(fosQueries, junQueries)]$sequence_id)
junFosSequenceIDs <- junFosSequenceIDs[!is.na(junFosSequenceIDs)]

# 3. get peakIDs of FOS/JUN sites
setkey(memeResultsCDSequences, "id")
junFosPeakIDs <- memeResultsCDSequences[junFosSequenceIDs]$name

# 4. now we can get the actual locations of these peaks
smallFragmentsSummits[junFosPeakIDs]

# 5. and cross them with H2A.Z peaks
junFosHits <- findOverlaps(smallFragmentsSummits[junFosPeakIDs], grNPIDRFiltered, minoverlap = 1)

# 6. and check in which gene promoters they are located
junFosPromoters <- findOverlaps(grNPIDRFiltered[subjectHits(junFosHits)], grPromoters)
kTGenes[grPromoters[subjectHits(junFosPromoters)]$extGene]

# 7. check what the differential expression looks like
library(ggrepel)
vpJunFosTargets <- ggplot(data = rTshH2AZ, 
                aes(x = b, y = -(log10(qval)))) +
  geom_point(color = "grey", size = 0.6) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")
vpJunFosTargets + geom_point(data = rTshH2AZ[grPromoters[subjectHits(junFosPromoters)]$extGene], color = "blue", size = 0.9)

table(rTshH2AZ[grPromoters[subjectHits(junFosPromoters)]$extGene]$b > 0)





