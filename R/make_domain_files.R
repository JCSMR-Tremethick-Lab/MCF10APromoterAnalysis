library(GenomicRanges)
library(rtracklayer)

files <- list.files(pattern = "consensus.txt")
genomeSize <- "~/Data/References/Annotations/Homo_sapiens/hg19/UCSC/hg19.chrom.sizes.sorted"
chromInfo <- read.csv(genomeSize, sep = "\t", header = F)

lapply(files, function(infile) {
  outfile <- paste(paste(unlist(strsplit(infile, "\\."))[c(1:5)], collapse = "."), "domains", sep = ".")
  tab1 <- read.csv(infile, sep = "\t", header = F)
  gr1 <- GenomicRanges::GRanges(seqnames = tab1$V1,
                                ranges = IRanges(tab1$V2, tab1$V3),
                                strand = "*")
  seqlevels(gr1) <- paste("chr", seqlevels(gr1), sep = "")
  seqlengths(gr1) <- chromInfo[chromInfo$V1 == seqlevels(gr1),]$V2
  gr1 <- trim(gr1)
  gr1 <- reduce(gr1)
  gr1$annotation <- "domain"
  gr1.gaps <- gaps(gr1)[-c(1:2)]
  gr1.gaps$annotation <- "gap"
  gr1 <- c(gr1, gr1.gaps)
  gr1 <- sort(gr1)
  df1 <- as(gr1, "data.frame")
  df1 <- df1[,c("seqnames", "start", "end", "annotation")]
  readr::write_delim(df1, outfile, delim = "\t", col_names = F)
})
