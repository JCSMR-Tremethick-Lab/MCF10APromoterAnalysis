library("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg19")

genome <- BSgenome.Hsapiens.UCSC.hg19

smallFragmentsDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/SmallFragments"

smallFragments <- lapply(list.files(smallFragmentsDir, pattern = ".bed$", full.names = T), function(x) {
  y <- rtracklayer::import(x)
})
names(smallFragments) <- unlist(lapply(strsplit(list.files(smallFragmentsDir, pattern = ".bed$"), "\\."), function(x) x[1]))

smallFragmentsSeq <- lapply(smallFragments, function(x) {
  return(getSeq(genome, resize(x, 500, fix = "center")))
})

lapply(names(smallFragmentsSeq), function(x){
  export(smallFragmentsSeq[[x]], con = file.path(smallFragmentsDir, paste(x, "_500bp.fa", sep = "")), format = "fasta")
})
