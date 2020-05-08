library(plotly)

source <- as.numeric(fig_wt_shz$WT) - 1
target <- 7 + as.integer(fig_wt_shz$shH2AZ) - 1
value <- fig_wt_shz$n
col <- unlist(mcf10awtCategories$color)


fig <- plot_ly(
  type = "sankey",
  orientation = "h",
  
  node = list(
    label = c(paste('WT_cluster_', c(1:7), sep = ''), paste('shH2AZ_cluster_', c(1:7), sep = '')),
    color = c(col, col),
    pad = 15,
    thickness = 20,
    line = list(
      color = "black",
      width = 0.5
    )
  ),
  
  link = list(
    source = source,
    target = target,
    value =  value
  )
)

fig <- fig %>% layout(
  title = "Basic Sankey Diagram",
  font = list(
    size = 10
  )
)

fig
