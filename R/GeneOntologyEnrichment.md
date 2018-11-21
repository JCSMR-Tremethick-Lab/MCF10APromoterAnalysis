GeneOntologyEnrichment
========================================================
author: Sebastian Kurscheid
date: 2018-10-16
autosize: true

Summary
========================================================

Purpose of the analysis:

- perform GO enrichment for each identified cluster in each of the MCF10A samples
- visual presentation of the results using ClusterProfiler

Slide With Code
========================================================


```r
dataDir <- "/home/sebastian/Data/Collaborations/FSU/PromoterSeqCap/PublicationFigures/Figure_2"
load(file.path(dataDir, "GOAnalysisObjects.rda"))
lapply(names(GOs$MCF10A_WT.clusters), function(x) {
  plotTitle <- x
  z <- GOs$MCF10A_WT.clusters[[x]]
  lapply(names(z), function(y){
    if(nrow(z[[y]]) > 0){
      clusterName <- y
      clusterProfiler::dotplot(z[[y]], title = paste(plotTitle, "-", clusterName, par(cex = 0.5)))
    }
  })
})
```

```
[[1]]
[[1]][[1]]
NULL

[[1]][[2]]
NULL

[[1]][[3]]
NULL


[[2]]
[[2]][[1]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-1.png)

```

[[2]][[2]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-2.png)

```

[[2]][[3]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-3.png)

```


[[3]]
[[3]][[1]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-4.png)

```

[[3]][[2]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-5.png)

```

[[3]][[3]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-6.png)

```


[[4]]
[[4]][[1]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-7.png)

```

[[4]][[2]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-8.png)

```

[[4]][[3]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-9.png)

```


[[5]]
[[5]][[1]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-10.png)

```

[[5]][[2]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-11.png)

```

[[5]][[3]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-12.png)

```


[[6]]
[[6]][[1]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-13.png)

```

[[6]][[2]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-14.png)

```

[[6]][[3]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-15.png)

```


[[7]]
[[7]][[1]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-16.png)

```

[[7]][[2]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-17.png)

```

[[7]][[3]]
```

![plot of chunk unnamed-chunk-1](GeneOntologyEnrichment-figure/unnamed-chunk-1-18.png)

Slide With Plot
========================================================

![plot of chunk unnamed-chunk-2](GeneOntologyEnrichment-figure/unnamed-chunk-2-1.png)
