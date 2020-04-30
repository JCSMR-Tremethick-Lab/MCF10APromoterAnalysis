rT.shH2AZ$target_id

selected <- rT.shH2AZ$target_id
up <- rT.shH2AZ[rT.shH2AZ$b > 0.5]$target_id
down <- rT.shH2AZ[rT.shH2AZ$b < -0.5]$target_id
setkey(FigS2Data, "extGene")

# all promoters
FigS2Data$wt_shz_de <- 'non_DE'
FigS2Data[up]$wt_shz_de <- 'up'
FigS2Data[down]$wt_shz_de <- 'down'
table(FigS2Data$wt_shz_de)

# for all promoters
FigS2bTab <- data.table::as.data.table(table("WT" = FigS2Data$wt.group, "shH2AZ" = FigS2Data$shZ.group, "DE" = FigS2Data$wt_shz_de))
FigS2bTab[,c('WT', 'shH2AZ', 'DE', 'N')] %>% group_by(WT, shH2AZ, DE) -> tab1


a1 <- alluvial(tab1[,c(1:2)], freq = tab1$N,
               col = c("red", "green", "grey")[match(tab1$DE, c("up", "down", "non_DE"))],
               alpha = c(1, 1, 0.3)[match(tab1$DE, c("up", "down", "non_DE"))])

# only expressed genes
setkey(FigS2Data, "extGene")
FigS2bTab <- FigS2Data[selected]
FigS2bTab <- data.table::as.data.table(table("WT" = FigS2bTab$wt.group, "shH2AZ" = FigS2bTab$shZ.group, "DE" = FigS2bTab$wt_shz_de))
FigS2bTab[,c('WT', 'shH2AZ', 'DE', 'N')] %>% group_by(WT, shH2AZ, DE) -> tab1

a1 <- alluvial(tab1[,c(1:2)], freq = tab1$N,
               col = c("red", "green", "grey")[match(tab1$DE, c("up", "down", "non_DE"))],
               alpha = c(1, 1, 0.3)[match(tab1$DE, c("up", "down", "non_DE"))])


tab1 <- data.table::as.data.table(table(FigS2Data$wt.group, FigS2Data$shZ.group))
tab1 %>% group_by(V1, V2) %>% summarise(n = sum(N)) -> tab2
FigS2bTab <- data.table::as.data.table(table("WT" = FigS2Data[rT.shH2AZ$target_id]$wt.group, "shH2AZ" = FigS2Data[rT.shH2AZ$target_id]$shZ.group))
FigS2bTab %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> FigS2bTab1

pdf(file = file.path(dataDir, "AlluvialPlotS2b.pdf"), paper = "a4r")
dev.off()

# only diff-expresse genes
FigS2bTab <- data.table::as.data.table(table("WT" = FigS2Data[selected]$wt.group, "shH2AZ" = FigS2Data[selected]$shZ.group))
FigS2bTab %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> FigS2bTab1
alluvial(FigS2bTab1[,c(1:2)], freq = FigS2bTab1$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)],
         cw = 0.05)






# only down-regulated genes
FigS2bTab <- data.table::as.data.table(table("WT" = FigS2Data[selected]$wt.group, "shH2AZ" = FigS2Data[selected]$shZ.group))
FigS2bTab %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> FigS2bTab1
alluvial(FigS2bTab1[,c(1:2)], freq = FigS2bTab1$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])

# only up-regulated genes
up <- rT.shH2AZ[rT.shH2AZ$b > 0.5]$target_id
FigS2bTab <- data.table::as.data.table(table("WT" = FigS2Data[selected]$wt.group, "shH2AZ" = FigS2Data[selected]$shZ.group))
FigS2bTab %>% group_by(WT, shH2AZ) %>% summarise(n = sum(N)) -> FigS2bTab1
alluvial(FigS2bTab1[,c(1:2)], freq = FigS2bTab1$n,
         col = mcf10awtCategories$color[match(as.integer(tab2$V1), mcf10awtCategories$group)])
