scRNA-Seq Tutorial Part 2: Exploration and Downstream Analysis
================

# Dependencies

``` r
library(Seurat)
```

    ## Loading required package: ggplot2

    ## Loading required package: cowplot

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    ## Loading required package: Matrix

``` r
library(princurve)
library(gam)
```

    ## Loading required package: splines

    ## Loading required package: foreach

    ## Loaded gam 1.15

``` r
library(pheatmap)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:Matrix':
    ## 
    ##     expand

``` r
library(viridis)
```

    ## Loading required package: viridisLite

``` r
library(RColorBrewer)
```

# Load processed data

``` r
pbmc <- readRDS("../data/pbmc_processed.rds")
```

# Exploring the structure of the data

Rather than blindly diving into analyzing the data, it’s always a good
idea to explore it a bit. We’ll use a couple dimensionality reduction
techniques to visualize our data and get a rough idea of what’s going
on.

Typically, I would do a couple visualizations before clustering, but now
that we’re confident with our QC/filtering/normalization it’s safe to do
it now (and I know it will be helpful for some conversation in the next
section)

# Clustering

We cluster on PCA space, but need to know how many PCs capture a
reasonable amount of
variation.

``` r
PCElbowPlot(object = pbmc, num.pc=50)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:30, 
    resolution = 0.5, print.output = 0, save.SNN = TRUE, random.seed=2018)
```

Note that the resolution is a user-defined parameter and can be tweaked
if necessary. Higher resolution finds more clusters.

random.seed is just the seed of the random number generator. This allows
us all to get the same result

## PCA

We had already done this, but we’ll explore it a little deeper this
time.

``` r
PCAPlot(object = pbmc)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Which genes drive a cells position along the first few PCs? Let’s look
at the genes with the highest loadings for the first few principal
components. We’ll also make it show up the top results for both positive
and negative
loadings.

``` r
VizPCA(object = pbmc, pcs.use = 1:2, do.balanced=T)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
VizPCA(object = pbmc, pcs.use = 3:4, do.balanced=T)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

If we take the genes with the top loadings and colour a PCA plot by
their expression, the relationship between PC loadings and expression
levels should become clear

``` r
#PC1 genes
FeaturePlot(pbmc, features.plot=c("RPL23A", "RPL3", "MNDA", "CSTA"), reduction.use="pca",
            cols.use=viridis(100))
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
#PC2 genes
FeaturePlot(pbmc, features.plot=c("CD79A", "MS4A1", "IL32", "TRAC"), reduction.use="pca",
            cols.use=viridis(100))
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

One problem with using PCA for complex data, however, is that it can
only capture a limited amount of variation that exists in the data. You
can start to see this in the first PCA plot of this section, where some
clusters on the left of the plot actually overlap quite a bit (which
would mean they’re crappy clusters if PCA was showing you everything).
Part of this is because PCA is a linear dimensionality technique, and
only so much variation can be captured by linear projections of the
data. This is apparent from the PC elbow plot above.

We’ve only been visualizing PC1 and 2, which will always capture the
most variation in the data, but you can see that PCs 3 and 4 also
capture a large amount of variation, and even PCs 5-7 have more than the
plateau that follows, suggesting there may be biological signal in there
as well. We’ve been ignoring these components in our exploration because
we haven’t been able to visualize them, but they likely contain
information that we’re interested in.

We can make more plots, visualizing these components if we want:

``` r
#PC3 genes from above
FeaturePlot(pbmc, features.plot=c("NKG7", "PRF1", "RPL21", "RPS18"), dim.1=3, dim.2=4,
            reduction.use="pca", cols.use=viridis(100))
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
#PC4 genes from above
FeaturePlot(pbmc, features.plot=c("TMSB10", "FCGR3A", "PTCRA", "SDPR"), dim.1=3, dim.2=4,
            reduction.use="pca", cols.use=viridis(100))
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

Wow, that’s interesting. It looks like PC3 is driven by a handful of
cells that are quite different from the rest of the population. When you
first see something like this, I strongly encourage you to not get super
excited and run to your PI about this cool, rare population of cells in
your data. It is possible that this is the case, but I would immediately
take a look at technical artifacts/variables that could drive
this.

``` r
FeaturePlot(pbmc, features.plot=c("nUMI", "nGene", "percent.mito"), dim.1=3, dim.2=4,
            reduction.use="pca", cols.use=viridis(100))
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
PCAPlot(pbmc, dim.1=3, dim.2=4)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

It doesn’t seem to be associated with any technical variable that we’ve
measured, so perhaps it is something interesting. From here, we could
take a look at the genes with the lowest PC4 loadings and see if they
correspond with some specific immune cell type. We could also look to
see if there are associated with something like apoptosis/cell stress.
For now, we won’t get too worked about this and will continue with our
exploration and see if this population pops out in clustering

To deal with the issue of missing out on information in our 2D plots, we
can take advtange of non-linear dimensionality reduction techniques.
Ultimately, any technique that tries to project complex relationships
into 2-3 dimensions will lose information in the process, but these
methods tend to improve quite a bit from PCA, often at the cost of
interpretability of the cells’ positions on the plot.

## tSNE

tSNE is slow and it’s often unreasonably to run it on the entire data
set. To get around this, we first perform PCA, find how many PCs provide
a significant amount of information (based on the elbow plot from
above), and use the matrix of PC positions for each cell across those
components. ie. we run tSNE on a matrix containing ~30 values (PC
coordinates) for each of the 8k cells, rather than the thousands of gene
expression
values.

``` r
pbmc <- RunTSNE(object = pbmc, dims.use = 1:30, do.fast = TRUE, perplexity=30)
```

``` r
TSNEPlot(pbmc)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## Force-directed graphs

We won’t execute this code because the igraph.to.gexf function can take
a long time and it always seems to crash R right after it finishes. This
said, I figured I should talk about graph layouts.

This code construct a shared nearest neighbour graph between all cells
in PCA space (using 20 PCs) and then we convert it to a gexf file that
can be loaded into the program Gephi, which can be used to visualize
networks.

``` r
#library(rgexf)
#pbmc@meta.data$Cluster <- pbmc@ident
#pbmc_knn <- scran::buildSNNGraph(t(pbmc@dr$pca@cell.embeddings[,1:20]), k=10, d=NA)
#gexf <- igraph.to.gexf(pbmc_knn)
#write.gexf(nodes = gexf$nodes, 
           #edges = gexf$edges[,c(2,3)], 
           #nodesAtt = pbmc@meta.data, 
           #output = "../data/pbmc_graph.gexf")
```

After importing the gexf file into Gephi and applying the Forceatlas2
layout algorithm to it, we get our data looks like this: \[PBMC Graph\]
(<https://github.com/dpcook/scrna_seq_workshop_2018/blob/master/figs/pbmc_forceatlas2.png>)

# Annotating clusters

Often after clustering our data we’re interested in finding out which
cell types each cluster corresponds to. Lucky, many of the cell types
comprising PBMCs have known markers that can serve as landmarks for us:

CD4 T Cells: IL7R CD8 T Cells: CD8A B Cells: MS4A1 CD14+ Monocytes:
CD14, LYZ FCGR3A+ Monocytes: FCGR3A, MS4A7 NK Cells: GNLY, NKG7
Dendritic Cells: FCER1A, CST3 Megakaryocytes: PPBP

So lets see if these genes mark clusters pretty
cleanly.

``` r
FeaturePlot(object = pbmc, features.plot = c("IL7R", "CD8A", "MS4A1", "CD14", "LYZ",
                                             "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A",
                                             "CST3", "PPBP"), 
            cols.use = viridis(100), reduction.use = "tsne",
            no.axes=T, pt.size=0.5)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
TSNEPlot(pbmc)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

# Identifying new markers

With scRNA-Seq, we’re often trying to find marker genes that define each
of our populations/clusters.

Two approaches: 1) genes with high fold change vs. all other
populations, or 2) genes with high classification power (ie. if we find
a cell that expresses a marker, there’s a good chance that it’s a given
cell type)

We can do this a cluster at a time with the FindMarkers function, or we
can find markers for all clusters in one run using FindAllMarkers. This
function offers 9 different tests that can be used. The current default
is the Wilcoxon rank sum test. A recent paper benchmarked a large amount
of differential expression tests, which you can check out
[here](LINK%20PAPER). Wilcoxon test ranked quite highly in their list.

The only.pos option, when TRUE, will only return positive markers of
each cluster. Min.pct sets what percentage of cells the gene must be
detected in. Thresh.use is a log2 fold-change requirement to be
considered a marker.

``` r
#This will take 10-15 mins
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.2, 
    thresh.use = 0.25, print.bar=T, random.seed=2018)
```

``` r
#How many markers per cluster
table(pbmc.markers$cluster)
```

    ## 
    ##   0   1   2   3   4   5   6   7   8   9  10  11  12 
    ## 628 103 125 102 232 199 181 206 667 657 445 895 501

Let’s just make a smaller data frame with the top markers (by logFC) of
each
group.

``` r
top_markers <- pbmc.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)
```

If you look at the fold changes for all the markers, you’ll see that we
found some pretty strong markers for most populations. Cluster 1 and 3
don’t seem to have genes that are a *lot* higher, but it’s still pretty
good.

Let’s visualize the top marker expression across the cells

``` r
FeaturePlot(object = pbmc, features.plot = top_markers$gene, 
            cols.use = viridis(100), reduction.use = "tsne",
            no.axes=T, pt.size=0.5)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Heatmaps are also a good way at looking at the stratification of
expression patterns across cells. Let’s make one, showing the top 10
markers of each cluster

``` r
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Or we can just print them out here to look at

``` r
pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
```

    ## # A tibble: 130 x 7
    ## # Groups:   cluster [13]
    ##    p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene         
    ##    <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>        
    ##  1     0      3.01 0.997 0.563         0 0       S100A8       
    ##  2     0      2.94 0.998 0.66          0 0       S100A9       
    ##  3     0      2.57 0.998 0.619         0 0       LYZ          
    ##  4     0      2.42 0.91  0.125         0 0       S100A12      
    ##  5     0      2.25 0.973 0.147         0 0       RP11-1143G9.4
    ##  6     0      2.13 0.985 0.209         0 0       FCN1         
    ##  7     0      2.04 0.979 0.164         0 0       MNDA         
    ##  8     0      2.00 0.998 0.385         0 0       TYROBP       
    ##  9     0      1.92 0.867 0.101         0 0       VCAN         
    ## 10     0      1.90 0.962 0.129         0 0       CSTA         
    ## # ... with 120 more rows

# Exploring dynamics within individual subtypes

## Monocytes

Monocytes are an abundant immune cell type in PBMCs that differentiate
into macrophages. CD68 and LYZ are common markers of monocytes.
CD14+/CD16- (FCGR3A) has been used to classify “classical
monocytes”

``` r
TSNEPlot(pbmc)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
FeaturePlot(object = pbmc, features.plot = c("CD14", "FCGR3A", "LYZ", "CD68"), 
            cols.use = viridis(100), reduction.use = "tsne",
            no.axes=T, pt.size=1)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

CD68 has been considered a pan-monocyte marker, and it seems to light up
in cluster 0 and 9. This corresponds nicely to the other markers as well
(Lyz is also expressed in dendritic cells, explaining its presence in
cluster 10). We also see that the cluster 0 is CD14+/CD16-, which is
consistent with classical monocytes. The CD16+ population (cluster 9) is
interesting though\!

Let’s filter our data down to only include clusters 0 and 9

``` r
monocytes <- SubsetData(pbmc, ident.use=c(0,9), do.center=T, do.scale=T)
```

    ## Scaling data matrix

We can now “zoom in” on this population and redo PCA/tSNE. We’ll have to
recalculate our highly-variable
genes.

``` r
monocytes <- FindVariableGenes(object = monocytes, mean.function = ExpMean, 
                          dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 5, 
                          y.cutoff = 0.5)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
monocytes <- RunPCA(object = monocytes, pc.genes = monocytes@var.genes, 
                    do.print = TRUE, genes.print=10, pcs.print=1:3,
                    pcs.compute=50)
```

    ## [1] "PC1"
    ##  [1] "FCGR3A"   "RHOC"     "CDKN1C"   "HES4"     "HLA-DPA1" "MTSS1"   
    ##  [7] "SIGLEC10" "ABI3"     "RPS19"    "MS4A7"   
    ## [1] ""
    ##  [1] "S100A8"        "S100A12"       "RP11-1143G9.4" "VCAN"         
    ##  [5] "MNDA"          "CD14"          "S100A6"        "FCN1"         
    ##  [9] "MS4A6A"        "SELL"         
    ## [1] ""
    ## [1] ""
    ## [1] "PC2"
    ##  [1] "FTH1"      "S100A4"    "AIF1"      "LST1"      "S100A6"   
    ##  [6] "CST3"      "DUSP1"     "HLA-DRA"   "LINC01272" "NEAT1"    
    ## [1] ""
    ##  [1] "CD2"    "CD3E"   "TRAC"   "CD3D"   "IL32"   "TRAT1"  "IL7R"  
    ##  [8] "CD3G"   "OCIAD2" "IFITM1"
    ## [1] ""
    ## [1] ""
    ## [1] "PC3"
    ##  [1] "FCGR3A"         "CDKN1C"         "S100A4"         "ALOX5AP"       
    ##  [5] "CKB"            "STXBP2"         "S100A12"        "RP11-1008C21.1"
    ##  [9] "MEG3"           "S100A6"        
    ## [1] ""
    ##  [1] "HLA-DRA"  "CD74"     "HLA-DRB1" "LGALS2"   "HLA-DPB1" "HLA-DQA1"
    ##  [7] "HLA-DMA"  "HLA-DPA1" "CST3"     "PID1"    
    ## [1] ""
    ## [1] ""

``` r
PCAPlot(object = monocytes)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

Perhaps as expected, the largest source of variation (PC1) separates the
cells based on their cluster. What’s kind of interesting is that it’s a
continuum and not distinct populations. Does this represent some
interesting biological spectrum? Let’s look at the expression of the
genes with the highest loadings for PC1 (printed out above)

``` r
FeaturePlot(object = monocytes, reduction.use="pca",
            features.plot = c("FCGR3A", "RHOC", "CDKN1C", "S100A8", "S100A12", "CD14"),
            cols.use=viridis(100)) 
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

There is actually quite a clear gradient\! Very cool. Let’s just check
the structure of the other
PCs

``` r
PCAPlot(object = monocytes, dim.1=1, dim.2=3)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
PCAPlot(object = monocytes, dim.1=1, dim.2=4)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->

``` r
PCAPlot(object = monocytes, dim.1=1, dim.2=5)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-26-3.png)<!-- -->

``` r
PCAPlot(object = monocytes, dim.1=1, dim.2=6)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-26-4.png)<!-- -->

PC3 is the only other PC that seems to have a shape that varies in a
structured way (doesn’t mean the others aren’t important, but PC3 seems
like it could be related to the phenotypic continuum along PC1)

### Pseudotime w/ principal curve

``` r
curve <- principal.curve(monocytes@dr$pca@cell.embeddings[,c(1,3)])
plot(curve)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

The lambda vector returned in the “curve” variable contains a value of
each cell’s *distance* along length of the curve. The tag vector is each
cell’s *order* along the curve. We’ll use the lambda value to represent
each cells’ pseudotemporal position

We’ll add the lambda values into the metadata table of our Seurat object

``` r
monocytes@meta.data$Pseudotime <- curve$lambda
```

And now we can plot our PCA, colouring each cell by their pseudotemporal
value

``` r
FeaturePlot(monocytes, reduction.use="pca", cols.use=viridis(100),
            features.plot="Pseudotime", dim.1=1, dim.2=3)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

### Pseudotime: Differential expression

Now that we’ve applied a pseudotime value to all cells, we can look for
genes whose expression changes as a function of pseudotemporal
progression.

We’ll use the LOESS method to model smooth, non-linear expression
patterns using the “gam” package.

``` r
exp <- monocytes@scale.data
pseudo_results <- apply(exp,1,function(g){
  data <- data.frame(Gene=g, Pseudotime=monocytes@meta.data$Pseudotime)
  tmp <- gam(g ~ lo(Pseudotime), data=data) #the lo() tells it to use LOESS
  p <- summary(tmp)[4][[1]][1,5]
  p
})
pseudo_results <- as.data.frame(pseudo_results)
colnames(pseudo_results) <- "gam.pval"
pseudo_results$q.value <- p.adjust(pseudo_results$gam.pval)
pseudo_results$gene <- rownames(pseudo_results)
nrow(filter(pseudo_results, q.value <= 0.05))
```

    ## [1] 1538

And we can visualize these genes in a heatmap, ordering cells by their
Pseudotime value

``` r
sig_genes <- filter(pseudo_results, q.value <= 0.05)$gene
cell_order <- curve$tag
exp <- monocytes@scale.data[sig_genes, cell_order]
hist(exp, breaks=50)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

Note that the range is quite wild. For a more-effective visualization on
the heatmap, we can set upper and lower-bounds. You’ll typically see
this set somewhere between 2-3

``` r
exp.z <- exp #one we can manipulate a little here
exp.z[exp.z>2] <- 2
exp.z[exp.z<(-2)] <- (-2)
```

``` r
bluered <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))#Make a color palette
heatmap <- pheatmap(exp.z,
                    color=bluered,
                    cluster_cols=F, #we manually ordered the cells--don't cluster
                    cutree_row=5,
                    cluster_rows=T, #we want genes to cluster
                    clustering_method="ward.D2",
                    show_rownames=F,
                    show_colnames=F)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

And we can retrieve the genes that fall in each cluster

``` r
gene_clusters <- as.data.frame(cutree(heatmap$tree_row, k=5))
colnames(gene_clusters) <- "Cluster"
gene_clusters$Cluster <- factor(gene_clusters$Cluster)
gene_clusters$Gene <- rownames(gene_clusters)
```

What is the order of the clusters down the heatmap?

``` r
unique(gene_clusters$Cluster[heatmap$tree_row$order])
```

    ## [1] 5 3 2 4 1
    ## Levels: 1 2 3 4 5

How many genes per cluster

``` r
summary(gene_clusters$Cluster)
```

    ##    1    2    3    4    5 
    ## 1126  190  168   37   17

We can save the results if we want to do gene ontology enrichment
analysis outside of R

``` r
write.csv(gene_clusters, file="../output/monocytes_gene_clusters.csv")
```

We can also look at the expression of individual genes throughout
pseudotime

``` r
plot_gene <- function(x){
  gene_data <- data.frame(Pseudotime=monocytes@meta.data$Pseudotime,
                          Exp = exp[x,])
  gene_plot <- ggplot(gene_data, aes(x=Pseudotime, y=Exp)) +
    geom_point(colour="black", alpha=0.25) +
    geom_smooth(colour="firebrick") +
    ylab(paste0(x, " Expression")) +
    theme_classic()
  ggsave(gene_plot, file=paste0("../figs/monocytes_",x,"_pseudotime.png"),
         width=3.5, height=2.5)
  gene_plot
}
```

``` r
head(arrange(pseudo_results, q.value))
```

    ##        gam.pval       q.value    gene
    ## 1  0.000000e+00  0.000000e+00  S100A9
    ## 2  0.000000e+00  0.000000e+00 S100A12
    ## 3  0.000000e+00  0.000000e+00  S100A8
    ## 4  0.000000e+00  0.000000e+00  FCGR3A
    ## 5  0.000000e+00  0.000000e+00     LYZ
    ## 6 1.218369e-310 2.014207e-306    RHOC

``` r
plot_gene("S100A9")
```

    ## `geom_smooth()` using method = 'gam'
    ## `geom_smooth()` using method = 'gam'

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

``` r
plot_gene("S100A12")
```

    ## `geom_smooth()` using method = 'gam'
    ## `geom_smooth()` using method = 'gam'

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-40-2.png)<!-- -->

``` r
plot_gene("S100A8")
```

    ## `geom_smooth()` using method = 'gam'
    ## `geom_smooth()` using method = 'gam'

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-40-3.png)<!-- -->

``` r
plot_gene("FCGR3A")
```

    ## `geom_smooth()` using method = 'gam'
    ## `geom_smooth()` using method = 'gam'

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-40-4.png)<!-- -->

``` r
plot_gene("LYZ")
```

    ## `geom_smooth()` using method = 'gam'
    ## `geom_smooth()` using method = 'gam'

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-40-5.png)<!-- -->

``` r
plot_gene("HLA-DPB1")
```

    ## `geom_smooth()` using method = 'gam'
    ## `geom_smooth()` using method = 'gam'

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-40-6.png)<!-- -->

#### Why tSNE can be misleading

I just wanted to bring up that tSNE seems to be peoples’ default
visualization for large single-cell datasets, but it can be misleading\!
Notice that tSNE does not pull out any sort of continuous structure in
the data
here.

``` r
monocytes <- RunTSNE(object = monocytes, dims.use = 1:10, do.fast = TRUE, perplexity=30)
```

``` r
TSNEPlot(monocytes)
```

![](02_downstream_analysis_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->
