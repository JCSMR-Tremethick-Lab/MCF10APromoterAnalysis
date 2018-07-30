# external snakemake script to extract summit sequences prior to meme processing
prepare_summit_sequences <- function(summitsFile, 
                                     genomePackage, 
                                     peaksFile, 
                                     summitsSeqFile, 
                                     summitsSeqWidth, 
                                     peaksMinPileUp, 
                                     peaksMinQval) {
    if (genomePackage == "BSgenome.Hsapiens.UCSC.hg19"){
        genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    } else {
      stop("wrong genome annotation. Please use BSgenome.Hsapiens.UCSC.hg19!")
    }
    summits <- rtracklayer::import(summitsFile)
    summits <- GenomeInfoDb::sortSeqlevels(summits, X.is.sexchrom = T)
    GenomeInfoDb::seqlevels(summits) <- GenomeInfoDb::seqlevels(genome)
    GenomeInfoDb::seqinfo(summits) <- GenomeInfoDb::seqinfo(genome)
    names(summits) <- gsub("TOTALcombined_", "", summits$name)
    names(summits) <- gsub("_000-125", "", names(summits)) 
    peaks <- data.table::fread(peaksFile)
    peaks$name <- gsub("TOTALcombined_", "", peaks$name)
    peaks$name <- gsub("_000-125", "", peaks$name)
    significantPeaks <-peaks[which(peaks$pileup > peaksMinPileUp & peaks$`-log10(qvalue)` > peaksMinQval)]
    summitsSeqs <- BSgenome::getSeq(genome, 
                                    GenomicRanges::resize(summits[significantPeaks$name], 
                                                          summitsSeqWidth, 
                                                          fix = "center"))
    rtracklayer::export(summitsSeqs, 
                        con = summitsSeqFile, 
                        format = "fasta")
}
save.image()
prepare_summit_sequences(summitsFile = snakemake@input[["summits"]],
                         peaksFile = snakemake@input[["peaks"]],
                         summitsSeqFile = snakemake@output[[1]],
                         genomePackage = snakemake@params[["BSgenome"]],
                         summitsSeqWidth = snakemake@params[["summitsSeqWidth"]],
                         peaksMinPileUp = snakemake@params[["peaksMinPileUp"]],
                         peaksMinQval = snakemake@params[["peaksMinQval"]])
