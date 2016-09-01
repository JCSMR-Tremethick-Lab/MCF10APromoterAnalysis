require(deepToolsUtils)
require(biomaRt)
require(tidyr)
require(RColorBrewer)
require(tibble)
require(snowfall)
require(sleuth)
require(tibble)
require(tximport)
require(readr)

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
  options(width = 137)
}
options(mc.cores = cpus)
setwd(lDir(pathPrefix, 
           "Data/Tremethick/Breast/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/R_Analysis/"))
dataPath <- lDir(pathPrefix, 
                 "Data/Tremethick/Breast/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/processed_data/GRCh38_ensembl84_ERCC/HTSeq/count/")
devPath <- "~/Development"

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

# load kallisto data with tximport and inspect via PCA -------------------------
base_dir <- paste(pathPrefix, 
                  "Data/Tremethick/Breast/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/processed_data", 
                  runConfig$references$version, 
                  "kallisto",
                  sep = "/")
sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
sample_id[c(1,2,7,8)] <- unlist(lapply(strsplit(sample_id[c(1,2,7,8)], "_"), function(x) paste(x[1], "wt", x[2], sep = "_")))
condition <- unlist(lapply(strsplit(sample_id, "_"), function(x) paste(x[1:2], collapse = "_")))
files <- paste(kal_dirs, "abundance.tsv", sep = "/")
names(files) <- sample_id


# import kallisto data ----------------------------------------------------
txi <- tximport::tximport(files, 
                          type = "kallisto",
                          geneIdCol = "ens_gene",
                          txIdCol = "target_id",
                          tx2gene = t2g,
                          reader = read_tsv)
# perform PCA for first inspection of data --------------------------------
pca1 <- ade4::dudi.pca(t(txi$abundance), scannf = F, nf = 6)
ade4::s.arrow(pca1$li)
ade4::s.class(pca1$li, fac = as.factor(condition))

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
s2c.mcf10a <- s2c[grep("MCF10A_|MCF10Ca1a_wt", s2c$condition),]
s2c.mcf10Ca1a <- s2c[grep("MCF10Ca1a", s2c$condition),]

s2c.list <- list(MCF10A = s2c.mcf10a,
                 MCF10Ca1a = s2c.mcf10Ca1a)

s2c.list <- lapply(names(s2c.list), function(x) {
  df <- s2c.list[[x]]
  df$sample <- as.character(df$sample)
  df$condition <- as.factor(as.character(df$condition))
  df$condition <- relevel(df$condition, paste(x, "wt", sep = "_"))
  df <- df[order(df$condition), ]
  return(df)
})
names(s2c.list) <- c("MCF10A", "MCF10Ca1a")

################################################################################
# actual processing using sleuth------------------------------------------------
sleuth_results_output <- "sleuthResults.rda"
if(!file.exists(sleuth_results_output)){
  results <- lapply(names(s2c.list), function(x){
    print(paste("Processing ", x, sep = ""))
    design <- model.matrix(~ condition, data = s2c.list[[x]])
    #-----------------------------------------------------------------------------
    # transcript-level DE
    so <- sleuth::sleuth_prep(s2c.list[[x]], ~ condition, target_mapping = t2g)
    so <- sleuth::sleuth_fit(so, formula = design)
    so <- sleuth::sleuth_fit(so, ~1, "reduced")
    so <- sleuth::sleuth_lrt(so, "reduced", "full")
    for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
      so <- sleuth::sleuth_wt(so, i)  
    }
    rt.list <- lapply(colnames(design)[grep("Intercept", colnames(design), invert = T)], function(x){
      rt <- sleuth::sleuth_results(so, x)
      rt <- rt[order(rt$qval),]
    })
    names(rt.list) <- colnames(design)[grep("Intercept", colnames(design), invert = T)]
    kt <- sleuth::kallisto_table(so, include_covariates = T)
    kt <- sleuth::kallisto_table(so, normalized = T, include_covariates = T)
    kt_wide <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
    rownames(kt_wide) <- kt_wide[,1]
    kt_wide <- kt_wide[,-1]
    #-----------------------------------------------------------------------------
    # gene-level DE  
    so.gene <- sleuth::sleuth_prep(s2c.list[[x]], ~ condition, target_mapping = t2g, aggregation_column = "ens_gene")
    so.gene <- sleuth::sleuth_fit(so.gene, formula = design)
    so.gene <- sleuth::sleuth_fit(so.gene, ~1, "reduced")
    so.gene <- sleuth::sleuth_lrt(so.gene, "reduced", "full")
    for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
      so.gene <- sleuth::sleuth_wt(so.gene, i)  
    }
    rt.gene.list <- lapply(colnames(design)[grep("Intercept", colnames(design), invert = T)], function(x){
      rt.gene <- sleuth::sleuth_results(so.gene, x)
      rt.gene <- rt.gene[order(rt.gene$qval),]
    })
    names(rt.gene.list) <- colnames(design)[grep("Intercept", colnames(design), invert = T)]
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
    kt_genes <- kt_genes[, c("ensembl_gene_id", s2c.list[[x]]$sample)]
    kt_genes <- tibble::as_tibble(merge(kt_genes, ensGenes, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id"))
    kt_wide <- tibble::as_tibble(tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm))
    return(list(sleuth_object = so,
                sleuth_results = rt.list,
                kallisto_table = kt,
                kallisto_table_wide = kt_wide,
                sleuth_results.gene = rt.gene.list,
                kallisto_table_genes = kt_genes))
  })
  save(results, file = sleuth_results_output)
} else {
  load(sleuth_results_output)
}

# re-formatting of list object --------------------------------------------
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
resultsCompressedBU <- resultsCompressed

resultsCompressed <- lapply(names(resultsCompressed), function(x) {
  resultsCompressed[[x]]$kallisto_table_wide <- resultsCompressed[[x]]$kallisto_table_wide[, c("target_id", s2c.list[[x]]$sample)]
  return(resultsCompressed[[x]])
})

names(resultsCompressed) <- names(s2c.list)

# prepare table for output ------------------------------------------------
tab1 <- resultsCompressed[[1]]$kallisto_table_genes
#########################################
# rename colnames to resolve sample mixup
colnames(tab1)[2:9] <- as.character(s2c.mcf10a$condition)[match(colnames(tab1)[2:9], as.character(s2c.mcf10a$sample))]
tab1 <- tab1[, c(1,3,6,4,5,2,7:18)]
colnames(tab1)[seq(2,9,2)] <- paste(colnames(tab1)[seq(2,9,2)], "rep1", sep = "_") 
colnames(tab1)[seq(3,9,2)] <- paste(colnames(tab1)[seq(3,9,2)], "rep2", sep = "_")

tab2 <- resultsCompressed[[2]]$kallisto_table_genes
colnames(tab2)[2:5] <- as.character(s2c.mcf10Ca1a$condition)[match(colnames(tab2)[2:5], as.character(s2c.mcf10Ca1a$sample))]
colnames(tab2)[seq(2,5,2)] <- paste(colnames(tab2)[seq(2,5,2)], "rep1", sep = "_") 
colnames(tab2)[seq(3,5,2)] <- paste(colnames(tab2)[seq(3,5,2)], "rep2", sep = "_")

tab1 <- as_tibble(merge(tab1, tab2[, c(1:3)], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T))
tab1 <- tab1[,c(1:9, 19:20, 10:18)]
write.csv(tab1, "MCF10A_RNA-Seq_results.csv")

# get list of EMT genes (from qPCR array) ---------------------------------
if(!file.exists("hsap.qPCRGenesTab.rda")){
  qPCRGeneList <- readLines(lDir(pathPrefix, "Data/Tremethick/EMT/ChIP-Seq/MDCK qPCR data/genelist.txt"))
  mart <- human <- useEnsembl(biomart = biomart, host = ensemblHost, dataset = dataset)
  hsap.qPCRGenesTab <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = qPCRGeneList, mart = mart)
  save(hsap.qPCRGenesTab, file = "hsap.qPCRGenesTab.rda")
} else {
  load("hsap.qPCRGenesTab.rda")
}

# make volcano plots for shZ & TGFb vs WT ---------------------------------
source(lDir(devPath, "JCSMR-Tremethick-Lab/Breast/R/volcano_plots_TGFb_shZ.R"))

# heatmap of samples using MCF10A_wt as reference -------------------------
source(lDir(devPath, "JCSMR-Tremethick-Lab/Breast/R/heatmap_MCF10A_shZ_TGFb_MCF10Ca1a"))

# prepare tables of DGEs --------------------------------------------------
sapply(resultsCompressed[["MCF10A"]]$sleuth_results.gene, function(x) {
  length(which(x$qval < 0.1 & abs(x$b) > 1))
})



