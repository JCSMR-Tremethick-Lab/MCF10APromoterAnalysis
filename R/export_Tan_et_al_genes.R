txdb <- AnnotationDbi::loadDb("~/Data/References/Annotations/Homo_sapiens/GRCh37_hg19_ensembl75/hsapiens_gene_ensembl_GRCh37_TxDB.sqlite")
sigEMTCells <- readr::read_tsv("~/Data/References/Annotations/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_cellLine.txt")

# load data from Ensembl
ensemblHost <- "grch37.ensembl.org"
dataset <- "hsapiens_gene_ensembl"
biomart <- "ensembl"
mart <- biomaRt::useEnsembl(biomart = biomart, dataset = dataset, host = ensemblHost)
attribs <- biomaRt::listAttributes(mart)
ensGenesSigEMTCells <- biomaRt::getBM(c("ensembl_gene_id", "external_gene_name"),
                                      filters = "external_gene_name",
                                      values = sigEMTCells$cellLine_sig,
                                      mart = mart)
ensGenesSigEMTCells <- merge(ensGenesSigEMTCells,
                             sigEMTCells,
                             by.x = "external_gene_name",
                             by.y = "cellLine_sig",
                             all.x = T)

# subset GRanges and write to BED files
gr.genes <- GenomicFeatures::genes(txdb)

deepToolsUtils::WriteGRangesToBED(gr.genes[ensGenesSigEMTCells[ensGenesSigEMTCells$epi_mes == "epi",]$ensembl_gene_id], 
                                  out_file = "~/Data/References/Annotations/Homo_sapiens/GRCh37_hg19_ensembl75/Tan_et_al_epithelial_genes.bed")

deepToolsUtils::WriteGRangesToBED(gr.genes[ensGenesSigEMTCells[ensGenesSigEMTCells$epi_mes == "mes",]$ensembl_gene_id], 
                                  out_file = "~/Data/References/Annotations/Homo_sapiens/GRCh37_hg19_ensembl75/Tan_et_al_mesenchymal_genes.bed")
