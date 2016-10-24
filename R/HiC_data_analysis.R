# HiTC analysis
require(HiTC)

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
setwd(lDir(pathPrefix, "Data/Tremethick/Breast/HiC/R_Analysis/"))
dataPath <- lDir(pathPrefix, "Data/Tremethick/Breast/HiC/HiCPro_output")
devPath <- "~/Development"

files1 <- list.files(lDir(dataPath, "/MCF10A_WT/hic_results/matrix/MCF10A_WT_sample1_rep1/raw/1000000"), full.names = T, recursive = T)

# import HiCPro processed data
sample1 <- importC(files1[3], files1[1], files1[1])
detail(sample1)
