require(rtracklayer)
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
require(GenomicRanges)

# external functions ------------------------------------------------------
source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")

# local functions ---------------------------------------------------------
lDir <- function(x, y){
  paste(x, y, sep = "/")
}

# get snakemake run configuration -----------------------------------------
runConfig <- jsonlite::fromJSON("~/Development/JCSMR-Tremethick-Lab/Breast/snakemake/configs/config_RNA-Seq.json")
refVersion <- "hg19"
runConfig$references[[refVersion]]$version

# global variables --------------------------------------------------------
if (refVersion == "hg38"){
  ensemblHost <- "mar2016.archive.ensembl.org"
  } else if (refVersion == "hg19"){
  ensemblHost <- "grch37.ensembl.org"
}

dataset <- "hsapiens_gene_ensembl"
biomart <- "ensembl"
colors <- RColorBrewer::brewer.pal(3, "Set2")

if (amILocal("JCSMR027564ML")){
  pathPrefix <- "~/mount/gduserv"
  mount <- system("mount", intern = T)
  if (length(grep("gduserv", mount)) == 0) {system("sshfs skurscheid@gduserv.anu.edu.au: ~/mount/gduserv/")}
  cpus <- 8
} else {
  pathPrefix <- "~"
  cpus <- 4
  options(width = 137)
}
options(mc.cores = cpus)

setwd(lDir(pathPrefix, 
           "Data/Tremethick/Breast/RNA-Seq_run2/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/R_Analysis/"))
devPath <- "~/Development"


# read in data ------------------------------------------------------------
dataPath <- lDir(pathPrefix, 
                 paste("Data/Tremethick/Breast/RNA-Seq_run2/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/processed_data/",runConfig$references[[refVersion]]$version,"HTSeq/count/", sep = ""))
files <- list.files(path = dataPath, full.names = T)
names(files) <- list.files(path = dataPath, full.names = F)

# preparing annotation data from Ensembl ----------------------------------
ensGenes_file <- paste("ensGenes_", runConfig$references[[refVersion]]$version, ".rda", sep = "")
ensTranscripts_file <- paste("ensTranscripts_", runConfig$references[[refVersion]]$version, ".rda", sep = "")
t2g_file <- paste("t2g_", runConfig$references[[refVersion]]$version, ".rda", sep = "")
myLength_file <- paste("mylength_", runConfig$references[[refVersion]]$version, ".rda", sep = "")
myGC_file <- paste("myGC_", runConfig$references[[refVersion]]$version, ".rda", sep = "")
myBiotypes_file <- paste("myBiotypes_", runConfig$references[[refVersion]]$version, ".rda", sep = "")
myChroms_file <- paste("myChroms_", runConfig$references[[refVersion]]$version, ".rda", sep = "")

if (!file.exists(ensGenes_file)){
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
    "gene_biotype",
    "entrezgene"),
    mart = mart)
  save(ensGenes, file = ensGenes_file)

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
  save(ensTranscripts, file = ensTranscripts_file)
  # create t2g object
  t2g <- ensTranscripts[, c("ensembl_transcript_id", 
                            "ensembl_gene_id", 
                            "external_gene_name", 
                            "version", 
                            "transcript_version")]
  if(refVersion == "hg38"){
    t2g$ensembl_transcript_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep = ".")
  }
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  save(t2g, file = t2g_file)
  
  mylength <- sapply(ensGenes$ensembl_gene_id, function(x){
    y <- ensTranscripts[which(ensTranscripts$ensembl_gene_id == x), ]
    y <- y[which.max(y$transcript_length), ]$transcript_length})
  save(mylength, file = myLength_file)
  mygc <- ensGenes$percentage_gc_content
  names(mygc) <- ensGenes$ensembl_gene_id
  save(mygc, file = myGC_file)
  mybiotypes <- ensGenes$gene_biotype
  names(mybiotypes) <- ensGenes$ensembl_gene_id
  save(mybiotypes, file = myBiotypes_file)
  mychroms <- data.frame(Chr = ensGenes$chromosome_name, GeneStart = ensGenes$start_position, GeneEnd = ensGenes$end_position)
  save(mychroms, file =myChroms_file)
  } else {
  load(ensGenes_file)
  load(ensTranscripts_file)
  load(myLength_file)
  load(myGC_file)
  load(myBiotypes_file)
  load(myChroms_file)
  load(t2g_file)
}

# load kallisto data with tximport and inspect via PCA -------------------------
base_dir <- paste(pathPrefix, 
                  "Data/Tremethick/Breast/RNA-Seq_run2/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/processed_data", 
                  runConfig$references[[refVersion]]$version, 
                  "kallisto",
                  sep = "/")
sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
#sample_id[c(1,2,7,8)] <- unlist(lapply(strsplit(sample_id[c(1,2,7,8)], "_"), function(x) paste(x[1], "wt", x[2], x[3], sep = "_")))
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
# MCF10A_shZ_rep2 is MCF10A_wt_rep2 and vice versa
condition <- as.character(condition[c(1,6,3:5,2,7:10)])
condition <- as.factor(condition)
condition <- relevel(condition, ref = "MCF10A_wt")
s2c <- data.frame(sample = sample_id, condition = condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
################################################################################
s2c.mcf10a <- s2c[grep("MCF10A_|MCF10Ca1a_wt", s2c$condition),]
s2c.mcf10a$condition <- as.factor(as.character(s2c.mcf10a$condition))
s2c.mcf10a$condition <- relevel(s2c.mcf10a$condition, ref = "MCF10A_wt")
s2c.mcf10a$sample <- as.character(s2c.mcf10a$sample)
s2c.mcf10Ca1a <- s2c[grep("MCF10Ca1a", s2c$condition),]
s2c.mcf10Ca1a$condition <- as.factor(as.character(s2c.mcf10Ca1a$condition))
s2c.mcf10Ca1a$condition <- relevel(s2c.mcf10Ca1a$condition, ref = "MCF10Ca1a_wt")
s2c.mcf10Ca1a$sample <- as.character(s2c.mcf10Ca1a$sample)

# make separate sets for different contrasts
# s2c.mcf10a_vs_mcf10a_shZ <- s2c.mcf10a[1:4,]
# s2c.mcf10a_vs_mcf10a_shZ$condition <- as.character(s2c.mcf10a_vs_mcf10a_shZ$condition)
# s2c.mcf10a_vs_mcf10a_shZ$condition <- as.factor(s2c.mcf10a_vs_mcf10a_shZ$condition)
# s2c.mcf10a_vs_mcf10a_shZ$condition <- relevel(s2c.mcf10a_vs_mcf10a_shZ$condition, "MCF10A_wt")
# s2c.mcf10a_vs_mcf10a_TGFb <- s2c.mcf10a[c(1,4,5:6),] 
# s2c.mcf10a_vs_mcf10a_TGFb$condition <- as.character(s2c.mcf10a_vs_mcf10a_TGFb$condition)
# s2c.mcf10a_vs_mcf10a_TGFb$condition <- as.factor(s2c.mcf10a_vs_mcf10a_TGFb$condition)
# s2c.mcf10a_vs_mcf10a_TGFb$condition <- relevel(s2c.mcf10a_vs_mcf10a_TGFb$condition, "MCF10A_wt")

# collate list()
s2c.list <- list(MCF10A = s2c.mcf10a,
                 MCF10Ca1a = s2c.mcf10Ca1a)
                 # MCF10A_vs_shZ = s2c.mcf10a_vs_mcf10a_shZ,
                 # MCF10A_vs_TGFb = s2c.mcf10a_vs_mcf10a_TGFb)

################################################################################
# actual processing using sleuth------------------------------------------------
sleuth_results_output <- paste("sleuthResults_", runConfig$references[[refVersion]]$version, ".rda", sep = "")

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
    target_mapping <- so$target_mapping
    rownames(target_mapping) <- target_mapping$target_id
    sfInit(parallel = T, cpus = cpus)
    sfExport("target_mapping")
    sfExport("kt_wide")
    l1 <- sfLapply(unique(target_mapping$ens_gene), function(x) {
      s <- apply(kt_wide[target_mapping[target_mapping$ens_gene == x, ]$target_id, ], 2, sum)
    })
    sfStop()
    names(l1) <- unique(target_mapping$ens_gene)
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

sleuth_resultsCompressed_file <- paste("sleuthResultsCompressed_", runConfig$references[[refVersion]]$version, ".rda", sep = "")

if(!file.exists(sleuth_resultsCompressed_file)){
  resultsCompressed <- lapply(names(results), function(x){
    results[[x]][grep("sleuth_object", names(results[[x]]), invert = T)]
  })
  names(resultsCompressed) <- names(results)
  
  resultsCompressed <- lapply(names(resultsCompressed), function(x){
    resultsCompressed[[x]][grep("kallisto_pca", names(resultsCompressed[[x]]), invert = T)]
  })
  names(resultsCompressed) <- names(results)
  save(resultsCompressed, file = sleuth_resultsCompressed_file)
} else {
  load(sleuth_resultsCompressed_file)
}

resultsCompressed <- lapply(names(resultsCompressed), function(x) {
  resultsCompressed[[x]]$kallisto_table_wide <- resultsCompressed[[x]]$kallisto_table_wide[, c("target_id", s2c.list[[x]]$sample)]
  return(resultsCompressed[[x]])
})
names(resultsCompressed) <- names(s2c.list)

# prepare table for output ------------------------------------------------
# load("tab1.rda")
tab1 <- resultsCompressed[[1]]$kallisto_table_genes

#########################################
# rename colnames to resolve sample mixup
# and prepare data for export
colnames(tab1)[2:9] <- as.character(s2c.mcf10a$condition)[match(colnames(tab1)[2:9], as.character(s2c.mcf10a$sample))]
tab1 <- tab1[, c(1,3,6,4,5,2,7:18)]
colnames(tab1)[seq(2,9,2)] <- paste(colnames(tab1)[seq(2,9,2)], "rep1", sep = "_") 
colnames(tab1)[seq(3,9,2)] <- paste(colnames(tab1)[seq(3,9,2)], "rep2", sep = "_")

tab2 <- resultsCompressed[[2]]$kallisto_table_genes
colnames(tab2)[2:5] <- as.character(s2c.mcf10Ca1a$condition)[match(colnames(tab2)[2:5], as.character(s2c.mcf10Ca1a$sample))]
colnames(tab2)[seq(2,5,2)] <- paste(colnames(tab2)[seq(2,5,2)], "rep1", sep = "_") 
colnames(tab2)[seq(3,5,2)] <- paste(colnames(tab2)[seq(3,5,2)], "rep2", sep = "_")

tab3 <- as_tibble(merge(tab1, tab2[, c(1:3)], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T))
tab3 <- tab3[,c(1:9, 19:20, 10:18)]
tab_exportFile <- paste("MCF10A_RNA-Seq_run2_results_", runConfig$references[[refVersion]]$version, ".csv", sep = "")
write.csv(tab3, tab_exportFile)
tab3 <- as.data.frame(tab3)

# load("tab3.rda")
tab4 <- data.frame(cbind(tab3[,1],
                   apply(tab3[,c(2,3)], 1, mean),
                   apply(tab3[,c(4,5)], 1, mean),
                   apply(tab3[,c(6,7)], 1, mean),
                   apply(tab3[,c(8,9)], 1, mean),
                   apply(tab3[,c(10,11)], 1, mean),
                   tab3[,c(12:20)]))
colnames(tab4)[1] <- "ensembl_gene_id"
colnames(tab4)[2:6] <- unlist(lapply(strsplit(colnames(tab3)[seq(2,11,2)], "_"), function(x) paste(x[1:2], collapse = "_")))

sapply(colnames(tab4)[2:6], function(x) {
  df <- data.frame(seqnames = tab4$chromosome_name, 
                   starts = as(tab4$start_position - 1, "integer"), 
                   ends = as(tab4$end_position, "integer"), 
                   names = tab4$ensembl_gene_id, 
                   scores = tab4[,x], 
                   strand = c("+", "-")[match(tab3$strand, c(1, -1))])
  write.table(df, 
              file = paste("~/Data/WorkingData/", x, ".bed", sep = ""), 
              quote = F, 
              sep = "\t", 
              row.names = F, 
              col.names = F)})

gr <- GRanges(seqnames = tab4$chromosome_name,
              IRanges(start = tab4$start_position,
                      end = tab4$end_position,
                      names = tab4$ensembl_gene_id),
              strand = c("+", "-")[match(tab3$strand, c(1, -1))],
              tab4[,c(2:6)])
gr <- sort(gr)


# make Venn diagrams of DE genes ------------------------------------------
sapply(resultsCompressed[["MCF10A"]]$sleuth_results.gene, function(x) table(x$qval < 0.1 & abs(x$b) > 1))

DEGenes <- lapply(names(resultsCompressed[["MCF10A"]]$sleuth_results.gene), function(x){
  dat <- resultsCompressed[["MCF10A"]]$sleuth_results.gene[[x]]
  up <- dat[which(dat$qval < 0.1 & dat$b > 1), "target_id"]
  down <- dat[which(dat$qval < 0.1 & dat$b < 1), "target_id"]
  return(list(up_regulated = up,
              down_regulated = down))
})

names(DEGenes) <- names(resultsCompressed[["MCF10A"]]$sleuth_results.gene)

venn.diagram(x = list(TGFb = DEGenes[["conditionMCF10A_TGFb"]][["down_regulated"]], 
                      shZ = DEGenes[["conditionMCF10A_shZ"]][["down_regulated"]]),
             filename = "MCF10A_TGFb_vs_shZ_down.png",
             imagetype = "png", 
             main = "Down-regulated genes in MCF10A\n TGFb & shZ knockdown")

venn.diagram(x = list(TGFb = DEGenes[["conditionMCF10A_TGFb"]][["up_regulated"]], 
                      shZ = DEGenes[["conditionMCF10A_shZ"]][["up_regulated"]]),
             filename = "MCF10A_TGFb_vs_shZ_up.png",
             imagetype = "png", 
             main = "Up-regulated genes in MCF10A\n TGFb & shZ knockdown")

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



