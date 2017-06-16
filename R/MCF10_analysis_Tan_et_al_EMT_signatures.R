# using Tan et al. 2014 EMT signatures ------------------------------------
# supplied as gene symbols :(
sigEMTCells <- readr::read_tsv("~/Data/References/Annotations/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_cellLine.txt")


# try to map to Ensembl IDs -----------------------------------------------
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

# additional genes to include ---------------------------------------------
addGenes <- data.frame(external_gene_name = c("H2AFZ", "TGFB1", "TGFB2", "TGFB3"),
                       ensembl_gene_id = c("ENSG00000164032", "ENSG00000105329", "ENSG00000092969", "ENSG00000119699"),
                       epi_mes = c("other", "mes", "mes", "mes"))
# attach to EMT signature
ensGenesSigEMTCells <- rbind(ensGenesSigEMTCells, addGenes)

# in fact, I need to use the differential transcript data in order to match it to UCSC IDs
ensTranscriptsSigEMTCells <-  biomaRt::getBM(c("ensembl_transcript_id", "external_gene_name", "ucsc"), 
                                             filters = "ensembl_gene_id",
                                             values = ensGenesSigEMTCells$ensembl_gene_id,
                                             mart = mart)
ucscTranscriptsSigEMTCells <- ensTranscriptsSigEMTCells[-which(ensTranscriptsSigEMTCells$ucsc == ""), ]
ucscTranscriptsSigEMTCells <- ucscTranscriptsSigEMTCells[!duplicated(ucscTranscriptsSigEMTCells$ucsc),]
ucscTranscriptsSigEMTCells <- merge(ucscTranscriptsSigEMTCells[,c("external_gene_name", "ucsc")], 
                                    ensGenesSigEMTCells[,c("external_gene_name", "epi_mes")],
                                    by.x = "external_gene_name",
                                    by.y = "external_gene_name",
                                    all.x = T)
ucscTranscriptsSigEMTCells <- ucscTranscriptsSigEMTCells[!duplicated(ucscTranscriptsSigEMTCells$ucsc),]

# volcano plot of cell line EMT genes -----------------------------
emtSignatureData <- lapply(names(resultsCompressed), function(x){
  lapply(names(resultsCompressed[[x]]$sleuth_results_gene), function(y){
    s <- resultsCompressed[[x]]$sleuth_results_gene[[y]]$target_id %in% ensGenesSigEMTCells$ensembl_gene_id
    dat <- resultsCompressed[[x]]$sleuth_results_gene[[y]][s,]
    dat <- merge(dat,
                 ensGenes[,c("ensembl_gene_id", "external_gene_name", "description")], 
                 by.x = "target_id", 
                 by.y = "ensembl_gene_id")
    dat <- merge(dat, 
                 ensGenesSigEMTCells[,c("ensembl_gene_id", "epi_mes")],
                 by.x = "target_id",
                 by.y = "ensembl_gene_id",
                 all.x = T)
    dat <- dat[order(dat$qval), ]
    d1 <- duplicated(dat$target_id)
    dat <- dat[!d1, ]
    return(list(dataTable = dat))
  })
})
names(emtSignatureData) <- names(resultsCompressed)
emtSignatureData <- lapply(names(emtSignatureData), function(x){
  print(x)
  names(emtSignatureData[[x]]) <- names(resultsCompressed[[x]]$sleuth_results_gene)
  return(emtSignatureData[[x]])
})
names(emtSignatureData) <- names(resultsCompressed)

x <- "MCF10A_vs_TGFb"
y <- "conditionMCF10A_TGFb"
dat <- emtSignatureData[[x]][[y]]$dataTable
xAxisMax <- max(abs(dat$b)) + 1
plot(dat$b,
     -log10(dat$qval), 
     axes = F, 
     xlab = "", 
     ylab = "", 
     frame = F,
     cex = 0.3,
     xlim = c(-round(xAxisMax, 0), round(xAxisMax,0)),
     pch = 16, main = paste("Volcano plot\nCondition: ", y, sep = ""))
# points(dat[which(-log10(dat$qval) >= 10), "b"], 
#        -log10(dat[which(-log10(dat$qval) >= 10), "qval"]),
#        col = "red", 
#        pch = 16, 
#        cex = 1.1)
points(dat[dat$epi_mes == "epi","b"],
       -log10(dat[dat$epi_mes == "epi","qval"]),
       col = "blue",
       pch = 16,
       cex = 1)
points(dat[dat$epi_mes == "mes","b"],
       -log10(dat[dat$epi_mes == "mes","qval"]),
       col = "green",
       pch = 16,
       cex = 1)
points(dat[dat$external_gene_name %in% c("H2AFZ", 'TGFB1', 'TGFB2', "TGFB3"),"b"],
       -log10(dat[dat$external_gene_name %in% c("H2AFZ", 'TGFB1', 'TGFB2', "TGFB3"),"qval"]),
       col = "red",
       pch = 16,
       cex = 1)

axis(2,
     pos = 0, 
     lwd = 3)
#                             at = c(seq(0,yAxisMax,10)))
axis(1,
     pos = 0, 
     lwd = 3,
     at = c(seq(-round(xAxisMax, 0), round(xAxisMax,0), 2)))
mtext("-log10(q-value)", 
      side = 2)
mtext("beta value",
      side = 1, 
      line = 2)
# only adding labels to genes with adjusted p-value <= 0.01
q_val_cutoff = 0.001
text(dat[which(dat$qval <= q_val_cutoff), "b"], -log10(dat[which(dat$qval <= q_val_cutoff), "qval"]), 
     labels = dat[which(dat$qval <= q_val_cutoff), "external_gene_name" ],
     cex = 0.7,
     pos = 4, offset = 0.3)
