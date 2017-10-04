library(clusterProfiler)
library(org.Hs.eg.db)
gmtfiles <- list.files("~/Data/References/Annotations/MSigDB/msigdb_v6.0_GMTs", pattern = "entrez", full.names = T)
gmtfiles <- gmtfiles[grep("all", gmtfiles)]
gmtfilesNames <- unlist(lapply(strsplit(gmtfiles, "/"), function(x) x[9]))
keytype <- "ENTREZID"
ontologies <- c("MF", "BP", "CC")
padjust <- "fdr"

x <- rt.gene.list$conditionMCF10A_shZ
ensGenes
x <- merge(x, subset(ensGenes, select = c("ensembl_gene_id", "entrezgene") ), by.x = "gene_id", by.y = "ensembl_gene_id")
table(x$qval < 0.1)

geneList <- subset(x, select = c("target_id", "qval", "b", "entrezgene"))
geneList$category <- as.character(NA)
geneList[qval < 0.1 & b > 0 ]$category <- "upregulated"
geneList[qval < 0.1 & b < 0]$category <- "downregulated"
table(geneList$category)
diffGenes <- geneList[qval < 0.1]
formula_res <- compareCluster(entrezgene ~ category, data=diffGenes, fun="enrichGO", org.Hs.eg.db)
dotplot(formula_res)

upGenes <- as.character(diffGenes[b > 0 & qval < 0.1]$entrezgene)
downGenes <- as.character(diffGenes[b < 0 & qval < 0.1]$entrezgene)
universe <- as.character(geneList$entrezgene)
length(universe)
length(downGenes)
length(upGenes)

GOEnrichment <- lapply(ontologies, function(x){
  up <- enrichGO(gene = upGenes,
                 universe = universe,
                 OrgDb = org.Hs.eg.db,
                 keytype = keytype,
                 ont = x,
                 pAdjustMethod = padjust,
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
  down <- enrichGO(gene         = downGenes,
                   universe = universe,
                   OrgDb         = org.Hs.eg.db,
                   keytype       = keytype,
                   ont           = x,
                   pAdjustMethod = padjust,
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
  return(list("up" = up, "down" = down))
})
names(GOEnrichment) <- ontologies

pdf("MCF10ca1ashZ_GO_Enrichment_plots.pdf", paper = "a4r")
lapply(names(GOEnrichment), function(x){
  lapply(names(GOEnrichment[[x]]), function(y){
    print(y)
    dotplot(GOEnrichment[[x]][[y]], title = paste("MCF10ca1a shZ knockdown: GO", x, y, sep = " "))
  })
})
dev.off()

# for the MSigDB analysis we use the whole list of ranked genes

geneList <- geneList[order(geneList$b, decreasing = T),]
geneList$entrezgene <- as.character(geneList$entrezgene)
gL <- geneList$b
names(gL) <- geneList$entrezgene

gmtList <- lapply(gmtfiles, function(y){
  gmt <- read.gmt(y)
  egmt <- GSEA(geneList = gL,
               TERM2GENE = gmt, 
               nPerm = 1000, 
               exponent = 1, 
               seed = T, 
               pAdjustMethod = "fdr")
})
names(gmtList) <- gmtfilesNames

x <- "MCF10ca1a_H2AZ_KD"
sapply(names(gmtList), function(y){
  df <- data.frame(gmtList[[y]])
  if(nrow(df) > 0) {
    fileName <- paste(x, y, "csv", sep = ".")
    write.csv(df, file = fileName)
    lapply(df$ID, function(z){
      if(df[z,]$qvalues < 0.1){
        setwd("GSEA_plots/")
        pdfFile <- paste(x, y, z, "pdf", sep = ".")
        pdf(pdfFile)
        gseaplot(gmtList[[y]], z, title = paste(z))
        dev.off()
        setwd("..")
      }
    })
  }
})



