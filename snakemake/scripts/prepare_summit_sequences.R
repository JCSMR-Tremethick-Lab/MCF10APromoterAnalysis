# external snakemake script to extract summit sequences prior to meme processing
prepare_summit_sequences <- function(summitsFile, peaksFile, summitsSeqFile, summitsSeqWidth, peaksMinPileUp, peaksMinQval) {
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  summits <- rtracklayer::import(summitsFile)
  summits <- GenomeInfoDb::sortSeqlevels(summits, X.is.sexchrom = T)
  GenomeInfoDb::seqlevels(summits) <- GenomeInfoDb::seqlevels(genome)
  GenomeInfoDb::seqinfo(summits) <- GenomeInfoDb::seqinfo(genome)
  names(summits) <- unlist(lapply(strsplit(summits$name, "_"), function(x) paste(x[c(2,5:6)], collapse = "_")))
  peaks <- data.table::fread(peaksFile)
  peaks$name <- unlist(lapply(strsplit(peaks$name, "_"), function(x) paste(x[c(2,5:6)], collapse = "_")))
  significantPeaks <-peaks[which(peaks$pileup > peaksMinPileUp & peaks$`-log10(qvalue)` > peaksMinQval)]
  summitsSeqs <- BSgenome::getSeq(genome, resize(summits[significantPeaks$name], summitsSeqWidth, fix = "center"))
  rtracklayer::export(summitsSeqs, con = summitsSeqFile, format = "fasta")
}

prepare_summit_sequences(summitsFile = snakemake@input[["summits"]], 
                         peaksFile = snakemake@input[["peaks"]],
                         summitsSeqFile = snakemake@output[[1]], 
                         summitSeqWidth = snakemake@config[["summitSeqWidth"]],
                         peakMinPileUp = snakemake@config[["peakMinPileUp"]],
                         peakMinQval = snakemake@config[["peakMinQval"]])

