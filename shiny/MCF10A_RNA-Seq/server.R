# This is the server logic of a Shiny web application

library(shiny)
library(DT)
library(tibble)


# load external data files ------------------------------------------------
load("data/sleuthResultsCompressed_GRCh37_hg19_ensembl75_V2.rda", .GlobalEnv)
load("data/ensGenes_GRCh37_hg19_ensembl75_ERCC.rda", .GlobalEnv)


# prepare data sets -----------------------------------------------------
# global DE data used for tables of all genes
globalDEData <- lapply(resultsCompressed, function(y){
  n <- names(y$sleuth_results_genes)
  s <- y$sleuth_results_genes[[n]]$target_id %in% ensGenes$ensembl_gene_id
  dat <- y$sleuth_results_genes[[n]][s,]
  dat <- merge(dat,
               ensGenes[,c("ensembl_gene_id", "external_gene_name", "description")], 
               by.x = "target_id", 
               by.y = "ensembl_gene_id")
  dat <- dat[order(dat$qval), ]
  dat <- tibble::as_tibble(dat)
  return(list(dataTable = dat))
})
names(globalDEData) <- unlist(lapply(resultsCompressed, function(x) names(x$sleuth_results_genes)))


# global gene expression data used for tables of all genes
globalGeneExpData <- lapply(resultsCompressed, function(y){
  s <- y$kallisto_table_genes$ensembl_gene_id %in% ensGenes$ensembl_gene_id
  dat <- y$kallisto_table_genes[s,]
  dat <- merge(dat,
               ensGenes[,c("ensembl_gene_id", "external_gene_name", "description")], 
               by.x = "ensembl_gene_id", 
               by.y = "ensembl_gene_id")
  dat <- tibble::as_tibble(dat)
  return(list(dataTable = dat))
})
names(globalGeneExpData) <- names(resultsCompressed)

#********************End of data preparation*******************************

# Actual Server logic -----------------------------------------------------
shinyServer(function(input, output, session) {
  selCols <- c("target_id", "qval", "b", "external_gene_name", "description")
  output$global_de_tgfb <- DT::renderDataTable(globalDEData[["conditionMCF10A_TGFb"]]$dataTable[,selCols],
                                               server = T)
  output$global_de_shZ <- DT::renderDataTable(globalDEData[["conditionMCF10A_shZ"]]$dataTable[,selCols],
                                              server = T)
  output$global_expression_tgfb <- DT::renderDataTable(globalGeneExpData[["MCF10A_vs_TGFb"]]$dataTable, 
                                                       server = T)
  output$global_expression_shZ <- DT::renderDataTable(globalGeneExpData[["MCF10A_vs_shZ"]]$dataTable, 
                                                       server = T)
})
