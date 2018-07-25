library(data.table)
library(ggparallel)


# Figure 1 ----------------------------------------------------------------
setwd("/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/sortingTSVS for Tremethick paper /Figure 1/")
sourceFiles <- list.files(".", pattern = ".tsv")
dataList <- lapply(sourceFiles, function(x){
  tab <- data.table::fread(x)
  return(tab)
})
names(dataList) <- unlist(lapply(strsplit(sourceFiles, "\\."), function(x) x[1]))
d1 <- do.call("data.table", dataList)

d1 <- subset(d1, select = c(1, grep("group1", colnames(d1))))
colNames <- colnames(d1)[2:8]
colNames <- unlist(lapply(strsplit(colNames, "_"), function(x) paste(x[1:2], collapse = "_")))
colnames(d1) <- c("gene", colNames)

ggparallel(list("A_Inp", "A_H2AZ", "CA1a_Inp", "CA1a_H2AZ", "shH2AZ_Inp", "TGFb_Inp", "TGFb_H2AZ"), data = d1)

ggparallel(list("A_Inp", "CA1a_Inp", "shH2AZ_Inp", "TGFb_Inp"), data = d1)
ggparallel(list("A_H2AZ", "CA1a_H2AZ", "TGFb_H2AZ"), data = d1)
ggparallel(list("A_H2AZ", "CA1a_H2AZ", "TGFb_H2AZ"), data = d1, method = "hammock", ratio = 0.1)
ggparallel(list("A_Inp", "CA1a_Inp", "shH2AZ_Inp", "TGFb_Inp"), data = d1, method = "hammock", ratio = 0.1)

# attempt at using plot.ly for interactive vis
library(plotly)
p <- m1 %>%
  plot_ly(width = 1920, height = 1080) %>% 
  add_trace(type = 'parcoords',
            line = list(showscale = TRUE,
                        reversescale = TRUE,
                        color = ~order.x,
                        colorscale = 'Jet',
                        cmin = 0,
                        cmax = 20000),
            dimensions = list(
              list(tickvals = c(1:7),
                   label = "group1",
                   values = ~group1.x),
              list(tickvals = c(1:7),
                   label = "group2",
                   values = ~group1.y),
              list(tickvals = c(1:7),
                   label = "group3",
                   values = ~group1)
            )
  )
as_widget(p)

