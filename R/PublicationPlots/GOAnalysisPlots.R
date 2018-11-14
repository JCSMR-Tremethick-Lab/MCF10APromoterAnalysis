# add GO analysis of MCF10A WT clusters -----------------------------------
library("clusterProfiler")
library("org.Hs.eg.db")

ff = function(x){ 
  if (class(x) == "list") 
    lapply(x, ff) 
  else if(class(x) == "enrichResult")
    if (nrow(x) > 0)
      TRUE
  else
    NULL
  else
    NULL
}

dataDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/Figure_2"

keyType <- "SYMBOL"
ontologies <- c("MF", "BP", "CC")
pAdjust <- "fdr"
pCutoff <- 0.10
analysis <- "clusterProfiler"

geneUniverse <- dt2$extGene

GO.MCF10A_WT.clusters <- lapply(levels(dt2$left.MCF10A_WT.category), function(x){
  geneList <- dt2[left.MCF10A_WT.category == x]$extGene
  MF <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "MF",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  BP <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "BP",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  CC <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "CC",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  return(list(MF = MF, CC = CC, BP = BP))
})

names(GO.MCF10A_WT.clusters) <- levels(dt2$left.MCF10A_WT.category)

load(file.path(dataDir, "GOAnalysisObjects.rda"))

pdf(file = file.path(dataDir, "GO.MCF10A_WT.clusters.pdf"), paper = "a4r")
lapply(names(GOs$MCF10A_WT.clusters), function(x) {
  plotTitle <- x
  z <- GOs$MCF10A_WT.clusters[[x]]
  lapply(names(z), function(y){
    if(nrow(z[[y]]) > 0){
      clusterName <- y
      dotplot(z[[y]], title = paste(plotTitle, "-", clusterName, par(cex = 0.5)))
    }
  })
})
dev.off()


# "right.MCF10A_TGFb.category"    -----------------------------------------
GO.MCF10A_TGFb.clusters <- lapply(levels(dt2$right.MCF10A_TGFb.category), function(x){
  geneList <- dt2[right.MCF10A_TGFb.category == x]$extGene
  MF <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "MF",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  BP <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "BP",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  CC <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "CC",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  return(list(MF = MF, CC = CC, BP = BP))
})

names(GO.MCF10A_TGFb.clusters) <- levels(dt2$right.MCF10A_TGFb.category)

pdf(file = file.path(dataDir, "GO.MCF10A_TGFb.clusters.pdf"), paper = "a4r")
lapply(names(GOs$MCF10A_TGFb.clusters), function(x) {
  plotTitle <- x
  z <- GOs$MCF10A_TGFb.clusters[[x]]
  lapply(names(z), function(y){
    if(nrow(z[[y]]) > 0){
      clusterName <- y
      dotplot(z[[y]], title = paste(plotTitle, "-", clusterName, par(cex = 0.5)))
    }
  })
})
dev.off()


# "right.MCF10A_shZ.category" ---------------------------------------------
GO.MCF10A_shZ.clusters <- lapply(levels(dt2$right.MCF10A_shZ.category), function(x){
  geneList <- dt2[right.MCF10A_shZ.category == x]$extGene
  MF <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "MF",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  BP <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "BP",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  CC <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "CC",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  return(list(MF = MF, CC = CC, BP = BP))
})

names(GO.MCF10A_shZ.clusters) <- levels(dt2$right.MCF10A_shZ.category)

pdf(file = file.path(dataDir, "GO.MCF10A_shZ.clusters.pdf"), paper = "a4r")
lapply(names(GOs$MCF10A_shZ.clusters), function(x) {
  plotTitle <- x
  z <- GOs$MCF10A_shZ.clusters[[x]]
  lapply(names(z), function(y){
    if(nrow(z[[y]]) > 0){
      clusterName <- y
      dotplot(z[[y]], title = paste(plotTitle, "-", clusterName, par(cex = 0.5)))
    }
  })
})
dev.off()


# "right.MCF10CA1A.category" ----------------------------------------------
GO.MCF10_CA1A.clusters <- lapply(levels(dt2$right.MCF10CA1A.category), function(x){
  geneList <- dt2[right.MCF10CA1A.category == x]$extGene
  MF <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "MF",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  BP <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "BP",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  CC <- enrichGO(gene = geneList,
                 universe = geneUniverse,
                 OrgDb = org.Hs.eg.db,
                 keyType = keyType,
                 ont = "CC",
                 pAdjustMethod = pAdjust,
                 pvalueCutoff = pCutoff)
  return(list(MF = MF, CC = CC, BP = BP))
})

names(GO.MCF10_CA1A.clusters) <- levels(dt2$right.MCF10CA1A.category)

pdf(file = file.path(dataDir, "GO.MCF10_CA1A.clusters.pdf"), paper = "a4r")
lapply(names(GOs$MCF10_CA1A.clusters), function(x) {
  plotTitle <- x
  z <- GOs$MCF10_CA1A.clusters[[x]]
  lapply(names(z), function(y){
    if(nrow(z[[y]]) > 0){
      clusterName <- y
      dotplot(z[[y]], title = paste(plotTitle, "-", clusterName, par(cex = 0.5)))
    }
  })
})
dev.off()



GOs <- list(GO.MCF10_CA1A.clusters, GO.MCF10A_shZ.clusters, GO.MCF10A_TGFb.clusters, GO.MCF10A_WT.clusters)
names(GOs) <- c("MCF10_CA1A.clusters", "MCF10A_shZ.clusters", "MCF10A_TGFb.clusters", "MCF10A_WT.clusters")
save(GOs, file = file.path(dataDir, "GOAnalysisObjects.rda"))


