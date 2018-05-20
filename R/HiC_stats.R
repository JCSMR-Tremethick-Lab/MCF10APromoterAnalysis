library(data.table)
library(tidyr)

setwd("/home/sebastian/mount/gduserv/Data/Tremethick/Breast/HiC/HiCPro_output_run3.1/hic_results/data")
l1 <- lapply(list.files(recursive = T, pattern = "RSstat"), function(x){
  dT <- data.table::fread(x)
  id <- unlist(strsplit(unlist(strsplit(x, split = "_hg38_ensembl84"))[1], "/"))[1]
  return(list(data = dT, id = id))
})

names(l1) <- unlist(lapply(strsplit(unlist(lapply(strsplit(list.files(recursive = T, pattern = "RSstat"), split = "_hg38_ensembl84"), function(x) x[1])), "/"), function(x) x[2]))

HiCStats <- data.table(sample = unlist(lapply(l1, function(x) x$id)), replicate = names(l1), do.call("rbind", lapply(l1, function(x) spread(data = x$data, key = "V1", value = "V2")))) 
subset(HiCStats, select = c(1, 3:12))[, lapply(.SD, sum), by = sample]

# sinlge sample single end statistics
setwd("/home/sebastian/mount/gduserv/Data/Tremethick/Breast/HiC/HiCPro_output_run3/logs/")
l2 <- lapply(list.files(recursive = T, pattern = "multiqc_bowtie2.txt"), function (x) {
  dT <- data.table::fread(x)
})
bowtieSE2Stats <- do.call("rbind", l2)
bowtieSE2Stats <- bowtieSE2Stats[!grep("local", bowtieSE2Stats$Sample)]

# single sample paired end statistics
setwd("/home/sebastian/mount/gduserv/Data/Tremethick/Breast/HiC/HiCPro_output_run3.1/bowtie_results/")
l3 <- lapply(list.files(recursive = T, pattern = "bwt2pairs.pairstat"), function(x){
  dTAbs <- data.table::fread(x, select = c(1,2))
  dTAbs$V1 <- factor(dTAbs$V1, levels =  dTAbs$V1)
  dTPerc <- data.table::fread(x, select = c(1,3))
  dTPerc$V1 <- factor(dTPerc$V1, levels =  dTPerc$V1)
  group <- unlist(lapply(strsplit(unlist(lapply(strsplit(x, split = "_hg38_ensembl84"), function(x) x[1])), "/"), function(x) x[2]))
  l <- list(absolute = dTAbs, percent = dTPerc, group = group)
  return(l)
})

names(l3) <- unlist(lapply(strsplit(unlist(lapply(strsplit(list.files(recursive = T, pattern = "bwt2pairs.pairstat"), split = "_hg38_ensembl84"), function(x) x[1])), "/"), function(x) x[3]))
bowtie2PEStats <- data.table::data.table(sample = unlist(lapply(l3, function(x) x$group)), replicate = names(l3), do.call("rbind", lapply(l3, function(x) spread(data = x$absolute, key = "V1", value = "V2"))))
subset(bowtie2PEStats, select = c(1, 3:12))[, lapply(.SD, sum), by = sample]

m1 <- merge(bowtie2PEStats, subset(HiCStats, select = c(2:12)), by.x = "replicate", by.y = "replicate")
write.csv(x = subset(m1, select = c(2:22))[, lapply(.SD, sum), by = sample], file = "HiC-Pro_stats_per_sample.csv")
write.csv(x = m1, file = "HiC-Pro_stats_per_replicate.csv")
