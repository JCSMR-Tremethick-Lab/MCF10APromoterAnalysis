require(deepToolsUtils)
require(biomaRt)
require(tidyr)
require(rtracklayer)
require(RColorBrewer)
require(tibble)
require(snowfall)
require(sleuth)
require(tibble)


# external functions ------------------------------------------------------
source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")

# local functions ---------------------------------------------------------
lDir <- function(x, y){
  paste(x, y, sep = "/")
}


# global variables --------------------------------------------------------
ensemblHost <- "mar2016.archive.ensembl.org"
dataset <- "hsapiens_gene_ensembl"
biomart <- "ensembl"
colors <- RColorBrewer::brewer.pal(3, "Set2")

if (amILocal("JCSMR027564ML")){
  pathPrefix <- "~/mount/gduserv"
  cpus <- 8
} else {
  pathPrefix <- "~"
  cpus <- 16
}
options(mc.cores = cpus)
setwd(lDir(pathPrefix, "Data/Tremethick/Breast/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/R_Analysis/"))
dataPath <- lDir(pathPrefix, "Data/Tremethick/Breast/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/processed_data/GRCh38_ensembl84_ERCC/HTSeq/count/")

files <- list.files(path = dataPath, full.names = T)
names(files) <- list.files(path = dataPath, full.names = F)

# get snakemake run configuration -----------------------------------------
runConfig <- jsonlite::fromJSON("~/Development/JCSMR-Tremethick-Lab/Breast/snakemake/configs/config_RNA-Seq.json")
runConfig$references$version

# preparing annotation data from Ensembl ----------------------------------
if (!file.exists("ensGenes.rda")){
  mart <- biomaRt::useEnsembl(biomart = biomart, dataset = dataset, host = ensemblHost)
                              attribs <- biomaRt::listAttributes(mart)
                              ensGenes <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                              "external_gene_name",
                              "chromosome_name",
                              "start_position",
                              "end_position",
                              "strand",
                              "band",
                              "description",
                              "percentage_gc_content",
                              "gene_biotype"),
                              mart = mart)
  save(ensGenes, file = "ensGenes.rda")

  # get Ensembl transcripts
  ensTranscripts <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                                  "ensembl_gene_id",
                                                  "transcript_length",
                                                  "version", 
                                                  "transcript_version",
                                                  "external_gene_name"),
                                   mart = mart,
                                   filter = "ensembl_gene_id",
                                   values = ensGenes$ensembl_gene_id)
  save(ensTranscripts, file = "ensTranscripts.rda")
  
  # create t2g object
  t2g <- ensTranscripts[, c("ensembl_transcript_id", 
                            "ensembl_gene_id", 
                            "external_gene_name", 
                            "version", 
                            "transcript_version")]
  t2g$ensembl_transcript_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep = ".")
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  save(t2g, file = "t2g.rda")
  
  mylength <- sapply(ensGenes$ensembl_gene_id, function(x){
    y <- ensTranscripts[which(ensTranscripts$ensembl_gene_id == x), ]
    y <- y[which.max(y$transcript_length), ]$transcript_length})
  save(mylength, file = "mylength.rda")
  mygc <- ensGenes$percentage_gc_content
  names(mygc) <- ensGenes$ensembl_gene_id
  save(mygc, file = "mygc.rda")
  mybiotypes <- ensGenes$gene_biotype
  names(mybiotypes) <- ensGenes$ensembl_gene_id
  save(mybiotypes, file = "mybiotypes.rda")
  mychroms <- data.frame(Chr = ensGenes$chromosome_name, GeneStart = ensGenes$start_position, GeneEnd = ensGenes$end_position)
  save(mychroms, file = "mychroms.rda")
  } else {
  load("ensGenes.rda")
  load("ensTranscripts.rda")
  load("mylength.rda")
  load("mygc.rda")
  load("mybiotypes.rda")
  load("mychroms.rda")
  load("t2g.rda")
}

# load kallisto data and normalize with sleuth ------------
base_dir <- paste(pathPrefix, "Data/Tremethick/Breast/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/processed_data", runConfig$references$version, "kallisto", sep = "/")
sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
sample_id[c(1,2,7,8)] <- unlist(lapply(strsplit(sample_id[c(1,2,7,8)], "_"), function(x) paste(x[1], "wt", x[2], sep = "_")))
condition <- unlist(lapply(strsplit(sample_id, "_"), function(x) paste(x[1:2], collapse = "_")))
################################################################################
# IMPORTANT!!!!!
# PCA shows that two samples were mislabelled
# inspection of the s.class plot reveals that most likely:
# MCF10A_shZ_rep2 is MCF10A_wt_rep2
condition <- as.character(condition[c(1,4,3,2,5:10)])
condition <- as.factor(condition)
condition <- factor(as.character(condition), levels = levels(condition)[c(3,2,1,5,4)])
s2c <- data.frame(sample = sample_id, condition = condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
################################################################################

s2c.TGFb_vs_WT <- rbind(s2c[grep("MCF10A_wt", s2c$condition),], s2c[grep("MCF10A_TGFb", s2c$condition),])
s2c.shZ_vs_WT <- rbind(s2c[grep("MCF10A_wt", s2c$condition),], s2c[grep("MCF10A_shZ", s2c$condition),])
s2c.Ca1aWT_vs_WT <- rbind(s2c[grep("MCF10A_wt", s2c$condition),], s2c[grep("MCF10Ca1a_wt", s2c$condition),])
s2c.Ca1AshZ_vs_WT <- rbind(s2c[grep("MCF10Ca1a_wt", s2c$condition),], s2c[grep("MCF10Ca1a_shZ", s2c$condition),])

s2c.list <- list(MCF10A_TGFb = s2c.TGFb_vs_WT,
                 MCF10A_shZ = s2c.shZ_vs_WT,
                 MCF10Ca1a_wt = s2c.Ca1aWT_vs_WT,
                 MCF10Ca1a_shZ = s2c.Ca1AshZ_vs_WT)

s2c.list <- S4Vectors::endoapply(s2c.list, function(x){
  x$sample <- as.character(x$sample)
  x$condition <- as.character(x$condition)
  x$condition <- factor(x$condition, levels = levels(as.factor(x$condition))[c(2,1)])
  return(x)
})

results <- lapply(names(s2c.list), function(x){
  print(paste("Processing ", x, sep = ""))
  if (grep(x, levels(s2c.list[[x]]$condition)) == 1) {
    cond <- as.character(s2c.list[[x]]$condition)
    s2c.list[[x]]$condition <- factor(cond, levels = levels(as.factor(cond))) 
  }
  design <- model.matrix(~ condition, data = s2c.list[[x]])
  #-----------------------------------------------------------------------------
  # transcript-level DE
  so <- sleuth::sleuth_prep(s2c.list[[x]], ~ condition, target_mapping = t2g)
  so <- sleuth::sleuth_fit(so, formula = design)
  so <- sleuth::sleuth_wt(so, paste("condition", x, sep = ""))
  so <- sleuth::sleuth_fit(so, ~1, "reduced")
  so <- sleuth::sleuth_lrt(so, "reduced", "full")
  rt <- sleuth::sleuth_results(so, paste("condition", x, sep = ""))
  rt <- rt[order(rt$qval),]
  kt <- sleuth::kallisto_table(so, include_covariates = T)
  kt <- kt[order(kt$target_id, kt$condition),]
  kt_wide <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
  rownames(kt_wide) <- kt_wide[,1]
  kt_wide <- kt_wide[,-1]
  # kt_pca <- ade4::dudi.pca(t(kt_wide), scannf = F, nf = 6)
  #-----------------------------------------------------------------------------
  # gene-level DE  
  so.gene <- sleuth::sleuth_prep(s2c.list[[x]], ~ condition, target_mapping = t2g, aggregation_column = "ens_gene")
  so.gene <- sleuth::sleuth_fit(so.gene, formula = design)
  so.gene <- sleuth::sleuth_wt(so.gene, paste("condition", x, sep = ""))
  so.gene <- sleuth::sleuth_fit(so.gene, ~1, "reduced")
  so.gene <- sleuth::sleuth_lrt(so.gene, "reduced", "full")
  rt.gene <- sleuth::sleuth_results(so.gene, paste("condition", x, sep = ""))
  rt.gene <- rt.gene[order(rt.gene$qval),]
  # gene-level expression is summed from transcript level data (sum(TPM))
  sfInit(parallel = T, cpus = cpus)
  sfExport("so")
  l1 <- sfLapply(unique(so$target_mapping$ens_gene), function(x) {
    s <- apply(kt_wide[so$target_mapping[so$target_mapping$ens_gene == x, "target_id"], ], 2, sum)
  })
  sfStop()
  names(l1) <- unique(so$target_mapping$ens_gene)
  kt_genes <- as.data.frame(do.call("rbind", l1))
  kt_genes$ensembl_gene_id <- rownames(kt_genes)
  kt_genes <- tibble::as_tibble(merge(kt_genes, ensGenes, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id"))
  kt_wide <- tibble::as_tibble(tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm))
  return(list(sleuth_object = so,
              sleuth_results = rt,
              kallisto_table = kt,
              kallisto_table_wide = kt_wide,
              #kallisto_pca = kt_pca,
              sleuth_results.gene = rt.gene,
              kallisto_table_genes = kt_genes))
})

names(results) <- names(s2c.list)
resultsCompressed <- lapply(names(results), function(x){
  results[[x]][grep("sleuth_object", names(results[[x]]), invert = T)]
})
names(resultsCompressed) <- names(results)

resultsCompressed <- lapply(names(resultsCompressed), function(x){
  resultsCompressed[[x]][grep("kallisto_pca", names(resultsCompressed[[x]]), invert = T)]
})
names(resultsCompressed) <- names(results)
save(resultsCompressed, file = "resultsCompressed.rda")

load("resultsCompressed.rda")
resultsCompressed <- resultsCompressedBU

resultsCompressed <- lapply(names(resultsCompressed), function(x) {
  resultsCompressed[[x]]$kallisto_table_wide <- resultsCompressed[[x]]$kallisto_table_wide[, c("target_id", s2c.list[[x]]$sample)]
  return(resultsCompressed[[x]])
})

names(resultsCompressed) <- names(s2c.list)

cor(kt_wide)
pdf("Pairwise_correlation_transcript_level.pdf")
pairs(kt_wide)
dev.off()


# save all results tables into one image file -----------------------------
save(list = ls()[grep("^rt.", ls())], file = "sleuth_results_tables.rda")


# get list of EMT genes (from qPCR array) ---------------------------------
qPCRGeneList <- readLines(lDir(pathPrefix, "Data/Tremethick/EMT/ChIP-Seq/MDCK qPCR data/genelist.txt"))
mart <- human <- useEnsembl(biomart = biomart, host = ensemblHost, dataset = dataset)
hsap.qPCRGenesTab <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = qPCRGeneList, mart = mart)


