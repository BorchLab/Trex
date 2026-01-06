# Trex

T Cell Receptor Sequence Embedding for Single-Cell Analysis

<!-- badges: start -->
[![R-CMD-check](https://github.com/BorchLab/Trex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BorchLab/Trex/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/BorchLab/Trex/graph/badge.svg)](https://app.codecov.io/gh/BorchLab/Trex?branch=master)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.borch.dev/uploads/screpertoire/articles/trex)
<!-- badges: end -->

## Introduction

<img align="right" src="https://github.com/BorchLab/Trex/blob/main/www/trex_hex.png" width="305" height="352">

Single-cell sequencing is now an integral tool in the field of immunology and oncology that allows researchers to couple RNA quantification and other modalities, like immune cell receptor profiling at the level of an individual cell. Towards this end, we developed the [scRepertoire](https://github.com/BorchLab/scRepertoire) R package to assist in the interaction of immune receptor and gene expression sequencing.

However, utilization of clonal indices for more complex analyses are still lacking, specifically in using clonality in embedding of single-cells. To this end, we developed Trex - an R package that uses deep learning to vectorize TCR sequences using sequence order or translating the sequence into amino acid properties.

### TCRex

If you are looking for the (very cool) TCR-epitope prediction algorithm **TCRex**, check out their website [here](https://tcrex.biodatamining.be/).

<img align="center" src="https://github.com/BorchLab/Trex/blob/dev/www/graphicalAbstract.png">

## Installation

### Bioconductor (Recommended)

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("Trex")
```

### Development Version

```r
devtools::install_github("BorchLab/Trex")
```

## System Requirements

Trex has been tested on R versions >= 4.4.0. Please consult the DESCRIPTION file for more details on required R packages. The package is specifically designed to work with single-cell objects that have had TCRs added using [scRepertoire](https://github.com/BorchLab/scRepertoire). Trex has been tested on macOS and Linux platforms.

### Python Environment

Trex uses [basilisk](https://bioconductor.org/packages/basilisk/) to automatically manage Python dependencies (TensorFlow, Keras). **No manual Python setup is required** - the environment is created automatically on first use.

### Model Downloads

By default, Trex will automatically pull deep learning models from a [Zenodo repository](https://zenodo.org/records/14636032) and cache them locally. The cache location follows platform conventions using `tools::R_user_dir("Trex", which = "cache")`.

## Model Information

<img align="center" src="https://github.com/BorchLab/Trex/blob/dev/www/training_info.png">

Each of the models available in Trex follow similar architecture with depth and width of input layers, epochs, batch size, and early stopping calls. The major difference is the size of the input layer, depending on the method chosen with **encoder.input**.

**Training data:**
- TRA: 905,611 sequences
- TRB: 1,276,462 sequences

**Architecture:** 256-128-30-128-256 with batch normalization

## Quick Start

### Basic Usage

```r
library(Trex)

# Generate embedding matrix
Trex_values <- maTrex(singleObject,
                      chains = "TRB",
                      method = "encoder",
                      encoder.model = "VAE",
                      encoder.input = "AF")
```

### Integration with Seurat

```r
# Run Trex and add to single-cell object
seuratObj <- runTrex(seuratObj,
                     chains = "TRB",
                     method = "encoder",
                     encoder.model = "VAE",
                     encoder.input = "AF",
                     reduction.name = "Trex")

# Generate UMAP from Trex embedding
seuratObj <- RunUMAP(seuratObj, reduction = "Trex", reduction.key = "Trex_")
```

### Weighted Nearest Neighbors

Combine Trex with RNA for multimodal analysis:

```r
# Remove TCR genes from variable features first
seuratObj <- quietTCRgenes(seuratObj)
seuratObj <- RunPCA(seuratObj)

# WNN approach
seuratObj <- FindMultiModalNeighbors(seuratObj,
                                     reduction.list = list("pca", "Trex"),
                                     dims.list = list(1:30, 1:20),
                                     modality.weight.name = "RNA.weight")

seuratObj <- RunUMAP(seuratObj,
                     nn.name = "weighted.nn",
                     reduction.name = "wnn.umap",
                     reduction.key = "wnnUMAP_")
```

## Encoding Methods

### Deep Learning (encoder)
- **VAE**: Variational autoencoder
- **AE**: Traditional autoencoder

### Input Types
- **AF**: Atchley factors (5 physicochemical properties)
- **KF**: Kidera factors (10 physicochemical properties)
- **OHE**: One-hot encoding

### Geometric Transformation
Alternative BLOSUM62-based encoding that doesn't require Python:

```r
Trex_values <- maTrex(singleObject,
                      chains = "TRB",
                      method = "geometric",
                      theta = pi)
```

## Documentation

Check out the full [vignette](https://www.borch.dev/uploads/screpertoire/articles/trex) for a complete tutorial.

## Bug Reports / Feature Requests

If you run into any issues or bugs please submit a [GitHub issue](https://github.com/BorchLab/Trex/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). Alternatively, an example with the internal **trex_example** would be extremely helpful.

Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/BorchLab/Trex/issues).

[Pull Requests](https://github.com/BorchLab/Trex/pulls) are welcome for bug fixes, new features, or enhancements.

## Citation

If using the Trex package, please cite our [manuscript](https://pubmed.ncbi.nlm.nih.gov/39164479/).
