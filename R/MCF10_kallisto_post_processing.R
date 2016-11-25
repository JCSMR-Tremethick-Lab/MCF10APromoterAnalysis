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
annotationVersion <- runConfig$references[[refVersion]]$version
annotationVersion <- annotationVersion[2]
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
  cpus <- 8
  options(width = 137)
}
options(mc.cores = cpus)

setwd(lDir(pathPrefix, 
           "Data/Tremethick/Breast/RNA-Seq_run2/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/R_Analysis/"))
devPath <- "~/Development"
annotationDataPath <- "~/Data/Tremethick/Breast/RNA-Seq_run2/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/R_Analysis/"

# read in data ------------------------------------------------------------
dataPath <- lDir(pathPrefix, 
                 paste("Data/Tremethick/Breast/RNA-Seq_run2/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/processed_data/", 
                       annotationVersion,
                       "/HTSeq/count/", sep = ""))
files <- list.files(path = dataPath, full.names = T)
names(files) <- list.files(path = dataPath, full.names = F)

# preparing annotation data from Ensembl ----------------------------------
ensGenes_file <- paste(annotationDataPath, "ensGenes_", annotationVersion, ".rda", sep = "")
ensTranscripts_file <- paste(annotationDataPath, "ensTranscripts_", annotationVersion, ".rda", sep = "")
t2g_file <- paste(annotationDataPath, "t2g_", annotationVersion, ".rda", sep = "")
myLength_file <- paste(annotationDataPath, "mylength_", annotationVersion, ".rda", sep = "")
myGC_file <- paste(annotationDataPath, "myGC_", annotationVersion, ".rda", sep = "")
myBiotypes_file <- paste(annotationDataPath, "myBiotypes_", annotationVersion, ".rda", sep = "")
myChroms_file <- paste(annotationDataPath, "myChroms_", annotationVersion, ".rda", sep = "")
annotationFileList <- list(ensGenes_file, ensTranscripts_file, t2g_file, myLength_file, myGC_file, myBiotypes_file, myChroms_file) 

if (all(sapply(annotationFileList, file.exists))){
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
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  if(length(grep("ERCC", ensGenes_file)) > 0){
    erccGenes <- import("~/Data/References/Transcriptomes/ERCC/ERCC92.gtf")
    erccGenes <- data.frame(target_id = erccGenes$gene_id, 
                            ens_gene = erccGenes$gene_id, 
                            ext_gene = erccGenes$gene_id)
    erccGenes$version <- 1
    erccGenes$transcript_version <- 1
    t2g <- rbind(t2g, erccGenes)
  }
  if(refVersion == "hg38"){
    t2g$ensembl_transcript_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep = ".")
  }
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
  save(mychroms, file = myChroms_file)
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
                  annotationVersion, 
                  "kallisto",
                  sep = "/")
sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
#sample_id[c(1,2,7,8)] <- unlist(lapply(strsplit(sample_id[c(1,2,7,8)], "_"), function(x) paste(x[1], "wt", x[2], x[3], sep = "_")))
condition <- c(rep("MCF10A_wt", 6), rep("MCF10A_shZ", 3), rep("MCF10A_TGFb", 3))
names(condition) <- names(files)
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
sd1 <- apply(txi$abundance, 1, sd)
summary(sd1)
# treat D6 and D8 samples separately
pcaD6 <- ade4::dudi.pca(t(txi$abundance[sd1 > 3, grep("D6", colnames(txi$abundance))]), scannf = F, nf = 6)
pdf(paste("PCA_MCF10A_wt_vs_TGFb_", annotationVersion, ".pdf", sep = ""))
ade4::s.arrow(pcaD6$li) # should remove MCF10AD6_3
ade4::s.class(pcaD6$li, fac = as.factor(condition[grep("D6", names(condition))]))
dev.off()

pcaD8 <- ade4::dudi.pca(t(txi$abundance[sd1 > 3, grep("D8", colnames(txi$abundance))]), scannf = F, nf = 6)
pdf(paste("PCA_MCF10A_wt_vs_shZ_", annotationVersion, ".pdf", sep = ""))
ade4::s.arrow(pcaD8$li) # should remove MCF10AD8_1
ade4::s.class(pcaD8$li, fac = as.factor(condition[grep("D8", names(condition))]))
dev.off()

# PCA shows that there is a WT sample at each timepoint with outlier behaviour
quarantinedSamples <- c("MCF10AD8_1", "MCF10AD6_3")
# prepare sample to condition table for sleuth processing -----------------
s2c <- data.frame(sample = sample_id, condition = condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c$sample <- as.character(s2c$sample)
# remove the quarantined samples
selected <- ! s2c$sample %in% quarantinedSamples
s2c <- s2c[selected,]
# make separate sets for different contrasts for shZ and TGFb
s2c.mcf10a_vs_mcf10a_shZ <- s2c[grep("D8", s2c$sample),]
s2c.mcf10a_vs_mcf10a_shZ$condition <- droplevels(s2c.mcf10a_vs_mcf10a_shZ$condition)
s2c.mcf10a_vs_mcf10a_shZ$condition <- relevel(s2c.mcf10a_vs_mcf10a_shZ$condition, "MCF10A_wt")
s2c.mcf10a_vs_mcf10a_TGFb <- s2c[grep("D6", s2c$sample),]
s2c.mcf10a_vs_mcf10a_TGFb$condition <- droplevels(s2c.mcf10a_vs_mcf10a_TGFb$condition)
s2c.mcf10a_vs_mcf10a_TGFb$condition <- relevel(s2c.mcf10a_vs_mcf10a_TGFb$condition, "MCF10A_wt")

# collate list()
s2c.list <- list(MCF10A_vs_shZ = s2c.mcf10a_vs_mcf10a_shZ,
                 MCF10A_vs_TGFb = s2c.mcf10a_vs_mcf10a_TGFb)

################################################################################
# actual processing using sleuth------------------------------------------------
analysis_version <- 2
sleuth_results_output <- paste("sleuthResults_", annotationVersion, "_V", analysis_version, ".rda", sep = "")

if(!file.exists(sleuth_results_output)){
  results <- lapply(names(s2c.list), function(x){
    design <- model.matrix(~ condition, data = s2c.list[[x]])
    print(paste("Processing ", x, " transcript-level analysis",sep = ""))
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
    return(list(sleuth_object = so,
                sleuth_results = rt.list,
                sleuth_results_genes = rt.gene.list,
                kallisto_table = kt,
                kallisto_table_wide = kt_wide,
                kallisto_table_genes = kt_genes))
  })
  names(results) <- names(s2c.list)
  save(results, file = sleuth_results_output)
} else {
  load(sleuth_results_output)
}

# re-formatting of list object --------------------------------------------

sleuth_resultsCompressed_file <- paste("sleuthResultsCompressed_", annotationVersion, "_V", analysis_version, ".rda", sep = "")

if(!file.exists(sleuth_resultsCompressed_file)){
  resultsCompressed <- lapply(names(results), function(x){
    results[[x]][grep("sleuth_object", names(results[[x]]), invert = T)]
  })
  names(resultsCompressed) <- names(results)
  
  resultsCompressed <- lapply(names(resultsCompressed), function(x){
    resultsCompressed[[x]][grep("kallisto_pca", names(resultsCompressed[[x]]), invert = T)]
  })
  resultsCompressed <- lapply(names(resultsCompressed), function(x) {
    resultsCompressed[[x]]$kallisto_table_wide <- resultsCompressed[[x]]$kallisto_table_wide[, c("target_id", s2c.list[[x]]$sample)]
    return(resultsCompressed[[x]])
  names(resultsCompressed) <- names(s2c.list)
  save(resultsCompressed, file = sleuth_resultsCompressed_file)
  })
} else {
  load(sleuth_resultsCompressed_file)
}



# diagnostic boxplot of ERCC spike in RNA abundances ----------------------
ERCCs <- resultsCompressed[["MCF10A_vs_shZ"]][["kallisto_table"]][grep("ERCC", resultsCompressed[["MCF10A_vs_shZ"]][["kallisto_table"]]$target_id),]
ERCCs <- rbind(ERCCs, resultsCompressed[["MCF10A_vs_TGFb"]][["kallisto_table"]][grep("ERCC", resultsCompressed[["MCF10A_vs_TGFb"]][["kallisto_table"]]$target_id),])
Transcripts <- resultsCompressed[["MCF10A_vs_shZ"]][["kallisto_table"]][grep("ENS", resultsCompressed[["MCF10A_vs_shZ"]][["kallisto_table"]]$target_id),]
Transcripts <- rbind(Transcripts, resultsCompressed[["MCF10A_vs_TGFb"]][["kallisto_table"]][grep("ENS", resultsCompressed[["MCF10A_vs_TGFb"]][["kallisto_table"]]$target_id),])

pdf("Diagnostic_boxplots.pdf", paper = "a4r")
p1 <- ggplot(data = ERCCs, mapping = aes(x = sample, y = log2(tpm + 1), fill = condition))
p1 + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "ERCC spike ins") # looks like spike ins have massive variability
p2 <- ggplot(data = Transcripts, mapping = aes(x = sample, y = log2(tpm + 1), fill = condition))
p2 + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "All transcripts")
dev.off()

# prepare table for output ------------------------------------------------
# load("tab1.rda")
tab1 <- resultsCompressed[[1]]$kallisto_table_genes


tab3 <- as_tibble(merge(tab1, tab2[, c(1:3)], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T))
tab3 <- tab3[,c(1:9, 19:20, 10:18)]
tab_exportFile <- paste("MCF10A_RNA-Seq_run2_results_", annotationVersion, ".csv", sep = "")
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
  keep <- grep("ENSG", hsap.qPCRGenesTab$ensembl_gene_id)
  hsap.qPCRGenesTab <- hsap.qPCRGenesTab[keep,]
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



