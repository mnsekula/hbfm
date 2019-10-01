# hbfm

**hbfm** is an R package for performing gene co-expression analysis in discrete single-cell RNA sequencing data. 

The details for the stochastic EM algorithm and MCMC sampler included in this package are provided in the "A sparse Bayesian factor model for the construction of gene co-expression networks from single-cell RNA sequencing count data" manuscript.

## Installation
This package can be installed from GitHub with the following code:

```{r, warning=FALSE}
# Install from GitHub using devtools
require(devtools)
devtools::install_github("mnsekula/hbfm")
```

Note: To install **hbfm** properly, the package `Rcpp` must be version 1.0.1 or later.

## Getting Started
The package vignette demonstrates how to use the **hbfm** package to construct a gene co-expression matrix from single-cell RNA-sequencing data. This vignette can be viewed online [here](http://htmlpreview.github.io/?https://github.com/mnsekula/hbfm/blob/master/hbfm.html).