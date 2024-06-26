
# DGEAR

<!-- badges: start -->
<!-- badges: end -->

The goal of DGEAR is to help the researchers find differentially expressed gene from microarray gene expression data or from the RNA seq count data the easiest way possible.

## Installation

You can install the development version of DGEAR like so:

``` r
# Simple installation:
install.packages("DGEAR")
```

## Library calling


``` r
library(DGEAR)
```
## Data format

DGEAR has it's own example data that can be accessible like so: 
``` r
# Data will be loaded with lazy loading and can be accessible when needed.
data("gene_exp_data")
head(gene_exp_data)
```
## How it works
Here's an example code to run the example dataset.

``` r
library(DGEAR)
data("gene_exp_data")
DGEAR(dataframe = gene_exp_data, con1 = 1, con2 = 10,
  exp1 = 11, exp2 = 20, alpha = 0.05, votting_cutoff = 2)
```
Try to find the DEGs from the example dataset with all the default arguments and understand it's working. 
