# script to extract short and long (nucleosomal) fragments 
# from alignments of PromoterSeqCap

library("data.table")
library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg19")

genome <- BSgenome.Hsapiens.UCSC.hg19
bedPEDir <- "/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/GRCh37_hg19_UCSC/bedtools/bedpe/"
fastaDir <- "/home/sebastian/Data/Tremethick/Breast/PromoterSeqCap/processed_data/fragments"
list.files(bedPEDir)

fileNames <- list.files(bedPEDir)
sampleNames <- unlist(lapply(strsplit(fileNames, "\\."), function(x) x[1]))
bedFiles <- list.files(bedPEDir, full.names = T)
lapply(1:length(bedFiles), function(x) {
  dT1 <- data.table::fread(bedFiles[x])
  gr1 <- GenomicRanges::makeGRangesFromDataFrame(dT1, seqnames.field = c("V1"), start.field = "V2", end.field = "V6", strand.field = "V9")
  n1 <- openssl::md5(dT1$V7)
  names(gr1) <- n1
  shortFrags <- getSeq(genome, gr1[width(gr1) < 121])
  nucleosomeFrags <- getSeq(genome, gr1[width(gr1) > 120 & width(gr1) < 200])
  BSgenome::export(shortFrags, con = file.path(fastaDir, paste(sampleNames[x], "short", "fa", sep = ".")), format = "fasta")
  BSgenome::export(nucleosomeFrags, con = file.path(fastaDir, paste(sampleNames[x], "long", "fa", sep = ".")), format = "fasta")
})
