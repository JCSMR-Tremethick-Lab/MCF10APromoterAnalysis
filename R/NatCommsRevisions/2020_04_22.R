library(data.table)

Fig1C <- data.table::fread('~/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/NatComms_revisions/Fig1C_sorting.tsv')
cluster_colors <- distinct(Fig1C[,c('group1', 'color')])
colors <- cluster_colors$color
names(colors) <- cluster_colors$group1
colors <- colors[sort(names(colors))]

colors_mod <- colors
colors_mod[c(1,4:7)] <- 'grey'

colors_mod2 <- colors
colors_mod2[c(1:4,6:7)] <- 'grey'


alluvial(FigS2bTab1[,c(1:2)], freq = FigS2bTab1$n,
         col = colors_mod[match(as.integer(tab1$V2), names(colors))])


alluvial(FigS2bTab1[,c(1:2)], freq = FigS2bTab1$n,
         col = colors_mod2[match(as.integer(tab1$V1), names(colors_mod2))],
         alpha = c(0.3, 0.3, 0.3, 0.3, 1, 0.3, 0.3)[match(tab1$V1, c(1:7))])


a1 <- alluvial(tab1[,c(1:2)], freq = tab1$N,
               col = c("red", "green", "grey")[match(tab1$DE, c("up", "down", "non_DE"))],
               alpha = c(1, 1, 0.3)[match(tab1$DE, c("up", "down", "non_DE"))])


# enricher
gois <- FigS2Data[FigS2Data$shZ.group == 5 & !FigS2Data$wt.group == 5]$extGene
universe <- FigS2Data$extGene
m_t2g <- msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name, human_gene_symbol)
em <- enricher(gois, universe = universe, pvalueCutoff = 0.01, TERM2GENE = m_t2g)

wkpEnrichWT <- lapply(c(1:7), function(x){
  g <- FigS2Data[FigS2Data$wt.group == x]$extGene
  g <- bitr(g, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
  c <- enricher(g$ENTREZID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  return(c)
})
wkpEnrichWT

wkpEnrichshZ <- lapply(c(1:7), function(x){
  g <- FigS2Data[FigS2Data$shZ.group == x]$extGene
  g <- bitr(g, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
  c <- enricher(g$ENTREZID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  return(c)
})
wkpEnrichshZ

wkpEnrichTGFb <- lapply(c(1:7), function(x){
  g <- FigS2Data[FigS2Data$tgfb.group == x]$extGene
  g <- bitr(g, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
  c <- enricher(g$ENTREZID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  return(c)
})
wkpEnrichTGFb

wkpEnrichCa1a <- lapply(c(1:7), function(x){
  g <- FigS2Data[FigS2Data$CA1a.group == x]$extGene
  g <- bitr(g, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
  c <- enricher(g$ENTREZID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  return(c)
})
wkpEnrichCa1a

