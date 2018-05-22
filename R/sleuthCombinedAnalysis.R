# combined analysis
library("sleuth")
library("data.table")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("Homo.sapiens")

Hsap <- Homo.sapiens
t2gHsap <- data.table(select(Hsap, 
                      keys(Hsap, "TXNAME"),
                      columns = c("TXNAME", "SYMBOL", "ENSEMBL", "GENEID"),
                      keytype = "TXNAME"))
t2gHsap <- t2gHsap[!is.na(t2gHsap$SYMBOL)]
t2gHsap <- dplyr::rename(t2gHsap, target_id = TXNAME, gene_id = GENEID, external_gene_name = SYMBOL, ensembl_gene_id = ENSEMBL)
t2gHsap <- t2gHsap[!duplicated(t2gHsap$target_id)]
t2gHsap[is.na(external_gene_name)]$external_gene_name <- t2gHsap[is.na(external_gene_name)]$target_id

TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
t2gTxDb <- data.table(select(TxDb, 
                             keys(TxDb, "TXNAME"), 
                             columns = c("TXNAME", "GENEID"), 
                             keytype = "TXNAME"))
t2gTxDb[is.na(t2gTxDb$GENEID)]$GENEID <- t2gTxDb[is.na(t2gTxDb$GENEID)]$TXNAME
t2gTxDb <- dplyr::rename(t2gTxDb, target_id = TXNAME, gene_id = GENEID)


# get snakemake run configuration -----------------------------------------
runConfig <- jsonlite::fromJSON("~/Development/JCSMR-Tremethick-Lab/Breast/snakemake/configs/config_RNA-Seq.json")
runNo <- names(runConfig$samples)[1]
refVersion <- "hg19"
annotationVersion <- runConfig$references[[refVersion]]$version
annotationVersion <- annotationVersion[3]
runID <- "combined" # first run

dataDir <- file.path("/home/sebastian/Data/Tremethick/Breast/RNA-Seq", runID, "processed_data/kallisto", annotationVersion)
devPath <- "~/Development"
annotationDataPath <- file.path("~/Data/Tremethick/Breast", runNo, runID ,"R_Analysis")
getwd()
# global variables --------------------------------------------------------
if (refVersion == "hg38"){
  ensemblHost <- "mar2016.archive.ensembl.org"
} else if (refVersion == "hg19"){
  ensemblHost <- "grch37.ensembl.org"
}

dataset <- "hsapiens_gene_ensembl"
biomart <- "ensembl"
colors <- RColorBrewer::brewer.pal(3, "Set2")
biotype <- "protein_coding"

if (!dir.exists(annotationDataPath)){
  dir.create(annotationDataPath)
  setwd(annotationDataPath)
} else {
  setwd(annotationDataPath)
}
getwd()

# file names for data output  ---------------------------------------------
analysis_version <- "1"
sleuth_results_output <- paste("sleuthResults_", annotationVersion, "_V", analysis_version, ".rda", sep = "")
sleuth_resultsCompressed_file <- paste("sleuthResultsCompressed_", annotationVersion, "_V", analysis_version, ".rda", sep = "")

# preparing annotation data from Ensembl ----------------------------------
if (length(grep("UCSC", annotationVersion)) > 0) {
  t2g_file <- paste(annotationDataPath, "t2g_", annotationVersion, ".rda", sep = "")
  mart <- biomaRt::useEnsembl(biomart = biomart, dataset = dataset, host = ensemblHost)
  annotationDataPath <- paste("~/Data/References/Annotations/Homo_sapiens/", annotationVersion, "/", sep = "")
  ucscTranscripts_file <- paste(annotationDataPath, "ucscTranscripts_", annotationVersion, ".rda", sep = "")
  load(ucscTranscripts_file)
  ucscTranscripts <- data.table::data.table(ucscTranscripts)
  ucscTranscripts <- ucscTranscripts[!is.na(GENEID)]
  ucscGenes <- data.table::data.table(biomaRt::getBM(attributes = c("external_gene_name",
                                             "description",
                                             "entrezgene"),
                              # filters = "entrezgene",
                              # values = ucscTranscripts$GENEID,
                              mart = mart))
  ucscGenes <- ucscGenes[!is.na(entrezgene)]
  ucscGenes$entrezgene <- paste("entrez:", ucscGenes$entrezgene, sep = "")
  t2g <- ucscTranscripts
  t2g$GENEID <- paste("entrez:", t2g$GENEID, sep = "")
  t2g <- dplyr::rename(t2g, target_id = TXNAME, gene_id = GENEID)
  t2g <- t2g[!is.na(t2g$gene_id)]
  t2g <- merge(t2g, ucscGenes, by.x = "gene_id", by.y = "entrezgene", all.x = T, all.y = F, allow.cartesian = F)
  setcolorder(t2g, c(2,1,3:8))
  save(t2g, file = t2g_file)
} else {
  stop()
}

# load kallisto data with tximport and inspect via PCA -------------------------
base_dir <- paste(pathPrefix, 
                  "Data/Tremethick/Breast",
                  runNo, 
                  runID, 
                  "processed_data",
                  annotationVersion, 
                  "kallisto",
                  sep = "/")

sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
# analysis of first MCF10A sequencing run:
# evidently shZ and WT samples were swapped. have to check which one shows H2A.Z knockdown
# looks like:
# MCF10A_shZ_rep2 is WT
# MCF10A_wt_rep2 is shZ
sample_id[grep("MCF10A_shZ_rep2|MCF10A_wt_rep2", sample_id)] <- c("MCF10A_wt_rep2", "MCF10A_shZ_rep2")
condition <- unlist(lapply(strsplit(sample_id, "_"), function(x) paste(x[1:2], collapse = "_")))
names(condition) <- sample_id
files <- paste(kal_dirs, "abundance.h5", sep = "/")
names(files) <- sample_id

txi <- tximport::tximport(files, 
                          type = "kallisto",
                          geneIdCol = gene_id,
                          txIdCol = target_id,
                          tx2gene = t2gTxDb[,c(1,2)],
                          ignoreTxVersion = F)

sd1 <- apply(txi$abundance, 1, sd)
summary(sd1)
pca1 <- ade4::dudi.pca(t(txi$abundance[sd1 > 3, ]), scannf = F, nf = 6)
ade4::s.arrow(pca1$li, boxes = F)

sample_id <- sample_id[grep(paste(names(files), collapse = "|"), sample_id)]
condition <- condition[grep(paste(names(files), collapse = "|"), names(condition))]
condition <- gsub("_1", "", condition)
condition <- gsub("_2", "", condition)
condition <- gsub("_3", "", condition)
ade4::s.class(pca1$li, fac = as.factor(condition))

s2c <- data.table::data.table(sample = sample_id, condition = condition)
s2c <- data.table::data.table(dplyr::mutate(s2c, path = kal_dirs))
s2c[grep("AD6|AD8", s2c$condition),]$condition <- "MCF10A_wt"
s2c$condition <- as.factor(s2c$condition)
s2c$condition <- relevel(s2c$condition, ref = "MCF10A_wt")
s2c <- s2c[condition != "MCF10Ca1a_shZ"]
s2c$condition <- droplevels(s2c$condition)
table(s2c$condition)

filter_function <- function(row, min_reads = 2, min_prop = 0.47) {
  mean(row >= min_reads) >= min_prop
}
log2_transform <- function(x) {log2(x + 0.5)}


if(!file.exists(sleuth_results_output)){
  design <- model.matrix(~ condition, data = s2c)
  #-----------------------------------------------------------------------------
  # transcript-level DE
  so <- sleuth::sleuth_prep(s2c, 
                            ~ condition,
                            target_mapping = t2g,
                            filter_fun = filter_function,
                            transformation_function = log2_transform,
                            max_bootstrap = 30, 
                            read_bootstrap_tpm = T, 
                            extra_bootstrap_summary = T)
  so <- sleuth::sleuth_fit(so, formula = design)
  so <- sleuth::sleuth_fit(so, ~1, "reduced")
  so <- sleuth::sleuth_lrt(so, "reduced", "full")
  for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
    so <- sleuth::sleuth_wt(so, i)  
  }
  rt.list <- lapply(colnames(design)[grep("Intercept", colnames(design), invert = T)], function(x){
    rt <- data.table::data.table(sleuth::sleuth_results(so, x))
    rt <- rt[order(rt$qval),]
  })
  names(rt.list) <- colnames(design)[grep("Intercept", colnames(design), invert = T)]
  kt <- data.table::data.table(sleuth::kallisto_table(so, normalized = T, include_covariates = T))
  kt_wide <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
  kt_wide <- data.table::as.data.table(kt_wide)
  data.table::setkey(kt_wide, "target_id")
  so.gene <- sleuth::sleuth_prep(s2c, ~ condition,
                                 target_mapping = t2gHsap, 
                                 aggregation_column = "external_gene_name",
                                 filter_fun = filter_function,
                                 transformation_function = log2_transform, 
                                 max_bootstrap = 30, 
                                 read_bootstrap_tpm = T, 
                                 extra_bootstrap_summary = T)
  so.gene <- sleuth::sleuth_fit(so.gene, formula = design)
  so.gene <- sleuth::sleuth_fit(so.gene, ~1, "reduced")
  so.gene <- sleuth::sleuth_lrt(so.gene, "reduced", "full")
  for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
    so.gene <- sleuth::sleuth_wt(so.gene, i)  
  }
  rt.gene.list <- lapply(colnames(design)[grep("Intercept", colnames(design), invert = T)], function(x){
    rt.gene <- sleuth::sleuth_results(so.gene, x)
    rt.gene <- rt.gene[order(rt.gene$qval),]
    rt.gene <- data.table::data.table(rt.gene)
  })
  names(rt.gene.list) <- colnames(design)[grep("Intercept", colnames(design), invert = T)]
  # gene-level expression is summed from transcript level data (sum(TPM))
  target_mapping <- data.table::as.data.table(so$target_mapping)
  kt.gene <- sleuth::kallisto_table(so.gene, use_filtered = T, normalized = T)
  kt.gene <- data.table::as.data.table(kt.gene)
  kt_wide.gene <- tidyr::spread(kt.gene[, c("target_id", "sample", "tpm")], sample, tpm)
  kt_wide.gene <- data.table::as.data.table(merge(kt_wide.gene, subset(t2g, select = c("target_id", "ens_gene")), all.x = TRUE, all.y = FALSE))
  cols.chosen <- as.character(s2c.list[[x]]$sample)
  kt_wide.gene <- kt_wide.gene[,lapply(.SD,sum),by=ens_gene, .SDcols = cols.chosen]
  data.table::setkey(kt_wide.gene, "ens_gene")
  kt_wide.gene <- data.table::as.data.table(merge(kt_wide.gene, unique(subset(ensGenes, select =  c("ensembl_gene_id", "external_gene_name", "description")), by = "ensembl_gene_id"), by.x = "ens_gene", by.y = "ensembl_gene_id", all.x = TRUE, all.y = FALSE))
  return(list(sleuth_object = so,
              sleuth_object_genes = so.gene,
              sleuth_results = rt.list,
              sleuth_results_genes = rt.gene.list,
              kallisto_table = kt,
              kallisto_table_wide = kt_wide,
              kallisto_table_genes = kt.gene,
              kallisto_table_genes_wide = kt_wide.gene))
}
  names(results) <- names(s2c.list)
  save(results, file = sleuth_results_output)
} else {
  load(sleuth_results_output)
}


# EMT signature analysis --------------------------------------------------
load(paste(annotationDataPath, "emtGenes_", annotationVersion, ".rda", sep = ""))
ucscTranscriptsSigEMTCells
ucscTranscriptsSigEMTCells <- dplyr::rename(ucscTranscriptsSigEMTCells, target_id = external_gene_name)
ucscGenesSigEMTCells <- ucscTranscriptsSigEMTCells[!duplicated(target_id), c("target_id", "epi_mes")]
setkey(ucscTranscriptsSigEMTCells, "target_id")
rt.gene.list$conditionMCF10A_TGFb[ucscGenesSigEMTCells$target_id]

key(rt.gene.list$conditionMCF10A_TGFb)
EMT.TGFb <- na.omit(rt.gene.list$conditionMCF10A_TGFb[ucscTranscriptsSigEMTCells[!duplicated(target_id),c("target_id", "epi_mes")], 
                                                      c("target_id", "qval", "b", "epi_mes")], 
                    cols = "qval")
EMT.TGFb$experiment <- "TGFb"
EMT.shZ <- na.omit(rt.gene.list$conditionMCF10A_shZ[ucscTranscriptsSigEMTCells[!duplicated(target_id), c("target_id", "epi_mes")], 
                                                    c("target_id", "qval", "b", "epi_mes")], 
                   cols = "qval")
EMT.shZ$experiment <- "H2A.Z KD"

volcanoData <- rbind(EMT.TGFb, EMT.shZ)
volcanoData$qval <- -log10(volcanoData$qval)
volcanoData <- dplyr::rename(volcanoData, "-log10qval" = qval, logFC = b)
pTGFb <- ggplot(volcanoData[experiment == "TGFb"], aes(x = logFC, y = `-log10qval`, color = epi_mes)) + geom_point() + ggtitle("TGFb vs WT")
pshZ <- ggplot(volcanoData[experiment == "H2A.Z KD"], aes(x = logFC, y = `-log10qval`, color = epi_mes)) + geom_point() + ggtitle("shZ vs WT")
table(volcanoData[experiment == "H2A.Z KD" & logFC > 0 & `-log10qval` > 1.5]$epi_mes)
table(volcanoData[experiment == "H2A.Z KD" & logFC < 0 & `-log10qval` > 1.5]$epi_mes)
table(volcanoData[experiment == "TGFb" & logFC > 0 & `-log10qval` > 1.5]$epi_mes)
table(volcanoData[experiment == "TGFb" & logFC < 0 & `-log10qval` > 1.5]$epi_mes)
table(volcanoData[experiment == "H2A.Z KD"]$epi_mes)

plot(EMT.TGFb$b, -log10(EMT.TGFb$qval))
plot(EMT.shZ$b, -log10(EMT.shZ$qval))



# subgroup analysis -------------------------------------------------------
# second rund D8 samples shZ vs WT
s2c.shZD8 <- s2c[grep("D8", s2c$path)]
s2c.shZD8$condition <- droplevels(s2c.shZD8$condition)
so.gene.shZD8 <- sleuth::sleuth_prep(s2c.shZD8, ~ condition,
                                     target_mapping = t2gHsap, 
                                     aggregation_column = "external_gene_name",
                                     filter_fun = filter_function,
                                     transformation_function = log2_transform, 
                                     max_bootstrap = 30, 
                                     read_bootstrap_tpm = T, 
                                     extra_bootstrap_summary = T)
design.shZD8 <- model.matrix(~ condition, data = s2c.shZD8)
so.gene.shZD8 <- sleuth::sleuth_fit(so.gene.shZD8, formula = design.shZD8)
so.gene.shZD8 <- sleuth::sleuth_fit(so.gene.shZD8, ~1, "reduced")
so.gene.shZD8 <- sleuth::sleuth_lrt(so.gene.shZD8, "reduced", "full")
for (i in colnames(design.shZD8)[grep("Intercept", colnames(design.shZD8), invert = T)]){
  so.gene.shZD8 <- sleuth::sleuth_wt(so.gene.shZD8, i)  
}
sleuth_live(so.gene.shZD8)
rt.gene.shZD8 <- data.table::data.table(sleuth::sleuth_results(so.gene.shZD8, "conditionMCF10AshZD8"))
setkey(rt.gene.shZD8, "target_id")
rt.gene.shZD8.EMT <- na.omit(rt.gene.shZD8[ucscGenesSigEMTCells], cols = "pval")
setkey(rt.gene.shZD8.EMT, "target_id")
# first run shZ vs WT
s2c.shZ <- s2c[grep("MCF10A_shZ_rep|MCF10A_wt_rep", s2c$sample)]
s2c.shZ$condition <- droplevels(s2c.shZ$condition)
so.gene.shZ <- sleuth::sleuth_prep(s2c.shZ, ~ condition,
                                     target_mapping = t2gHsap, 
                                     aggregation_column = "external_gene_name",
                                     filter_fun = filter_function,
                                     transformation_function = log2_transform, 
                                     max_bootstrap = 30, 
                                     read_bootstrap_tpm = T, 
                                     extra_bootstrap_summary = T)
design.shZ <- model.matrix(~ condition, data = s2c.shZ)
so.gene.shZ <- sleuth::sleuth_fit(so.gene.shZ, formula = design.shZ)
so.gene.shZ <- sleuth::sleuth_fit(so.gene.shZ, ~1, "reduced")
so.gene.shZ <- sleuth::sleuth_lrt(so.gene.shZ, "reduced", "full")
for (i in colnames(design.shZ)[grep("Intercept", colnames(design.shZ), invert = T)]){
  so.gene.shZ <- sleuth::sleuth_wt(so.gene.shZ, i)  
}
sleuth_live(so.gene.shZ)
rt.gene.shZ <- data.table::data.table(sleuth::sleuth_results(so.gene.shZ, "conditionMCF10A_shZ"))
setkey(rt.gene.shZ, "target_id")
rt.gene.shZ.EMT <- na.omit(rt.gene.shZ.EMT[ucscGenesSigEMTCells], cols = "pval")
setkey(rt.gene.shZ.EMT, "target_id")
table(rt.gene.shZ.EMT[b < 0 & qval < 0.1]$epi_mes)

######################################
# COMMENT: shZD8 clearly did not work
#####################################


# compare the two knock-down experiments ----------------------------------
i1 <- intersect(rt.gene.shZD8.EMT$target_id, rt.gene.shZ.EMT$target_id)
cor(rt.gene.shZD8.EMT[i1]$b, rt.gene.shZ.EMT[i1]$b)


# plot all histone genes --------------------------------------------------
histoneGenes <- data.table::fread("~/Development/JCSMR-Tremethick-Lab/adhoc/data.csv")
histoneGenes <- data.table::data.table(histoneGenes[,c("hgnc_symbol", "Group")])
tmp <- data.table::data.table(biomaRt::getBM(attributes = c("external_gene_name", "description", "ensembl_gene_id", "hgnc_symbol"),
                                      filters = "hgnc_symbol",
                                      values = histoneGenes$hgnc_symbol,
                                      mart = mart))
setkey(tmp, "hgnc_symbol")
setkey(histoneGenes, "hgnc_symbol")
histoneGenes <- histoneGenes[tmp]
key(histoneGenes)
setkey(histoneGenes, "external_gene_name")

plotData <- na.omit(kt.gene[histoneGenes])[na.omit(kt.gene[histoneGenes])$tpm > 5]
ggplot(data = plotData, aes(x = target_id, y = tpm)) + 
  geom_boxplot(aes(colour = condition)) +
  facet_wrap(~Group + condition, scales = "free_x", ncol = 6) +
  theme(axis.text = element_text(angle = 45), legend.position = "none")

