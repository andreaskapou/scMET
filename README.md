# scMET
Bayesian modelling of DNA methylation heterogeneity at single-cell resolution

## Background
Here we introduce __scMET__, a Bayesian framework for the analysis of single cell DNA methylation data. This modelling approach combines a hierarchical beta-binomial specification with a generalised linear model framework with the aim of capturing biological overdispersion and overcome data sparsity by sharing information across cells and genomic features.

To disentangle technical from biological variability and overcome data sparsity, scMET couples a hierarchical BB model with a GLM framework (Fig.1a-b). For each cell i and region j, the input for scMET is the number of CpG sites that are observed to be methylated (__Y__) and the total number of sites for which methylation status was recorded (__n__). The BB model uses feature-specific mean parameters __mu__ to quantify overall DNAm across all cells and biological _overdispersion_ parameters __gamma__ as a proxy for cell-to-cell DNAm heterogeneity. These parameters capture the amount of variability that is not explained by binomial sampling noise, which would only account for technical variation.

The GLM framework is incorporated at two levels. Firstly, to introduce feature-specific covariates __x__ (e.g. CpG density) that may explain differences in mean methylation mu across features. Secondly, similar to [Eling2018](https://pubmed.ncbi.nlm.nih.gov/30172840/), we use a non-linear regression framework to capture the mean-overdispersion trend that is typically observed in high throughput sequencing data, such as scBS-seq (Fig.1c). This trend is used to derive _residual overdispersion_ parameters __epsilon__: a measure of cell-to-cell variability that is not confounded by mean methylation. Feature-specific parameters are subsequently used for: (i) feature selection, to identify highly variable features (HVFs) that drive cell-to-cell epigenetic heterogeneity (Fig.1d) and (ii) differential methylation testing, to highlight features that show differences in DNAm mean or variability between specified groups of cells (Fig.1e). 

Overview of the `scMET` model is shown below:

![](inst/figures/scmet-motivation.png)

## Installation
```R
# Install stable version from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scMET")

## Or development version from Github
# install.packages("remotes")
remotes::install_github("andreaskapou/scMET")
```

### Installation issue requiring the V8 library (outdated)
__This is an old issue with Rstan, keeping here for reference__.

The scMET package depends heavily on `Rstan`, whose newer version depends on the V8 library (see this issue: [https://github.com/stan-dev/rstan/issues/831](https://github.com/stan-dev/rstan/issues/831)). For users that don't have the system-level V8 dependency pre-installed, there are two approaches.

1. (Recommended) Download a static libv8 library when installing on Linux:
```R
Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
install.packages("V8")
```
see this blogpost for more details: [https://ropensci.org/blog/2020/11/12/installing-v8/](https://ropensci.org/blog/2020/11/12/installing-v8/).

2. If the above approach does not work, then install the `rstan` package from this branch that removes the dependency and the code that calls it and makes no other changes:
```R
install_github("makoshark/rstan", ref="develop", subdir="rstan/rstan")
```
and then proceed installing scMET as above. If during the installation you still get the error about installing V8, try and set `dependencies = FALSE` when calling the `install_github` function.

## Online vignette
scMET is not yet part of the Bioconductor. Until then, an online vignette can be found in
[https://rpubs.com/cakapourani/scmet-analysis](https://rpubs.com/cakapourani/scmet-analysis).


# Citation:
Kapourani, C. A., Argelaguet, R., Sanguinetti, G., & Vallejos, C. A. (2021). scMET: Bayesian modeling of DNA methylation heterogeneity at single-cell resolution. Genome biology, 22(1), 1-21.

[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02329-8](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02329-8)
