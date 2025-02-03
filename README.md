# Trex
Using TCR sequences for sequence embedding

<!-- badges: start -->
[![R-CMD-check](https://github.com/BorchLab/Trex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BorchLab/Trex/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/BorchLab/Trex/graph/badge.svg)](https://app.codecov.io/gh/BorchLab/Trex?branch=master)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.borch.dev/uploads/screpertoire/articles/trex)
<!-- badges: end -->

## Introduction

<img align="right" src="https://github.com/BorchLab/Trex/blob/main/www/trex_hex.png" width="305" height="352">

Single-cell sequencing is now a integral tool in the field of immunology and oncology that allows researchers to couple RNA quantification and other modalities, 
like immune cell receptor profiling at the level of an individual cell. Towards this end, we developed the [scRepertoire](https://github.com/BorchLab/scRepertoire) 
R package to assist in the interaction of immune receptor and gene expression sequencing. However, utilization of clonal indices for more complex analyses are still lacking, specifically in using clonality in embedding of single-cells. To this end, we developed an R package that uses deep learning to vectorize TCR sequences using order or translating the sequence into amino acid properties.

### TCRex
If you are looking for the (very cool) TCR-epitope prediction algorithm **TCRex**, check out their website [here](https://tcrex.biodatamining.be/).

<img align="center" src="https://github.com/BorchLab/Trex/blob/dev/www/graphicalAbstract.png">

# System requirements 

Trex has been tested on R versions >= 4.0. Please consult the DESCRIPTION file for more details on required R packages - it is specifically designed to work with single-cell objects that have had TCRs added using [scRepertoire](https://github.com/BorchLab/scRepertoire). Trex has been tested on OS X and Windows platforms.

**keras** is necessary to use the autoencoder function (this includes the set up of the tensorflow environment in R):

```r
##Install keras
install.packages("keras")

##Setting up Tensor Flow
library(reticulate)
conda_create("r-reticulate") ##If first time using reticulate
use_condaenv(condaenv = "r-reticulate", required = TRUE)
library(keras)
install_keras()
```

An alternative to this approach above (especially if you want to avoid conda) is to use reticulate to generate a virtualenv, using ```virtualenv_create()``` and subsequently installing the above python packages using ```virtualenv_install()```.

# Model Information 
<img align="center" src="https://github.com/BorchLab/Trex/blob/dev/www/training_info.png">

Each of the models available in Trex follow similar architecture with depth and width of input layers, epochs, batch size, and early stopping calls. The major difference is the size of the input layer, depending on the method chosen with **encoder.input**. 

# Installation

To run Trex, open R and install Trex from github: 

```r
devtools::install_github("BorchLab/Trex")
```
# Usage/Demos

Trex should be able to be run in popular R-based single-cell workflows, including Seurat and Bioconductor/Single-Cell Experiment formats.

## Quick Start 

Check out this [vignette](https://www.borch.dev/uploads/screpertoire/articles/trex) for a quick start tutorial. 

## Autoencoded Matrix

The Trex algorithm allows users to select TCR-based metrics to return autoencoded values to be used in dimensional reduction. If single-cell objects are not filtered for T cells with TCR,  `maTrex()` will still return values, however TREX_1 will be based on the disparity of TCR-containing and TCR-non-containing cells based on the Trex algorithm. 

```r
library(Trex)
my_trex <- maTrex(singleObject)
```

## Seurat or Single-Cell Experiment

You can run Trex within your Seurat or Single-Cell Experiment workflow. **Importantly** `runTrex()` will automatically filter single-cells that do not contain TCR information in the meta data of the single-cell object. 

```r
seuratObj_Tonly <- runTrex(seuratObj, #The single cell object
                       chains = "TRB", #Use of "TRA" or "TRB" 
                       method = "encoder", #Use "encoder" for CNNs or "geometric" geometric-based transformation
                       encoder.model = "VAE" #"VAE" (variational autoencoder) or "AE" (autoencoder)
                       encoder.input = "AF" #Inputs into encoder - "AF", "KF", "both", "OHE
                       reduction.name = "Trex", #Name designation for slot in single-cell object)
                   
seuratObj_Tonly <- runTrex(seuratObj, reduction.name = "Trex")
```

## After Running Trex

From here, you can generate a tSNE/UMAP using the Trex values, similar to the PCA values based on variable gene expression.

```r
seuratObj <- RunTSNE(seuratObj, reduction = "Trex",  reduction.key = "Trex_")
seuratObj <- RunUMAP(seuratObj, reduction = "Trex",  reduction.key = "Trex_")
```

If using Seurat package, the Trex embedding information and gene expression PCA can be used to find the [Weighted Nearest Neighbors](https://pubmed.ncbi.nlm.nih.gov/34062119/). Before applying the WNN approach, best practice would be to remove the TCR-related genes from the list of variable genes and rerunning the PCA analysis. 

### Recalculate PCA without TCR genes with queitTCRgenes() function in Trex.
```r
seuratObj <- quietTCRgenes(seuratObj)
seuratObj <- RunPCA(seuratObj)
```

### Running WNN approach
```r
seuratObj <- FindMultiModalNeighbors(seuratObj, 
                                     reduction.list = list("pca", "Trex"), 
                                     dims.list = list(1:30, 1:20), 
                                     modality.weight.name = "RNA.weight")
                                     
seuratObj <- RunUMAP(seuratObj, 
                     nn.name = "weighted.nn", 
                     reduction.name = "wnn.umap", 
                     reduction.key = "wnnUMAP_")
```
***
## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/BorchLab/Trex/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 
Alternatively, an example with the internal **trex_example** would 
be extremely helpful.

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/BorchLab/Trex/issues).

#### [Pull Requests](https://github.com/BorchLab/Trex/pulls) are welcome for bug fixes, new features, or enhancements.

## Citation
If using the Trex package, please cite our [manuscript](https://pubmed.ncbi.nlm.nih.gov/39164479/). This is also a good place to find more information about the models.


