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
library(viridis)
```

    ## Loading required package: viridisLite

# Load processed data

``` r
pbmc <- readRDS("../data/pbmc_processed.rds")
```

# Exploring the structure of the data

Rather than blindly diving into analyzing the data, it’s always a good
idea to explore it a bit. We’ll use a couple dimensionality reduction
techniques to visualize our data and get a rough idea of what’s going
on.

## PCA

We had already done this, but we’ll explore it a little deeper this
time.

## tSNE

## Diffusion Map

# Clustering

# Identifying markers

Note that with complex populations, differential expression is a little
different than typical. If you just run a generic differential
expression test, it is likely that the majority of genes will be
differentially expressed because we would simply be looking for
variation anywhere across the population. A more useful direction in
this type of “cell type atlasing” experiment is to identify genes whose
expression uniquely defines a given population.

Traditional differential expression can be beneficial when making
comparisons across subsets of cell types, or between a given population
across two experimental conditions (eg. Wild-type vs. knockout)

# Pseudotime analysis

I want to see if I can fit this in here. Any trajectories exist in
PBMCs?

# Functional enrichment (GO Term and Pathway analysis)