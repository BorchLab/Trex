#' Example single-cell dataset with TCR information
#'
#' A Seurat v5 object containing 3,200 single T cells from a CD3+ sorted
#' severe COVID-19 patient sample. This dataset is used for examples and
#' vignettes demonstrating Trex functionality.
#'
#' @format A Seurat object with:
#' \describe{
#'   \item{RNA assay}{Gene expression counts}
#'   \item{ADT assay}{Antibody-derived tag counts}
#'   \item{pca}{PCA reduction from RNA}
#'   \item{apca}{PCA reduction from ADT}
#'   \item{CTaa}{CDR3 amino acid sequence (from scRepertoire)}
#'   \item{CTgene}{V and J gene usage (from scRepertoire)}
#' }
#'
#' @source
#' Derived from GSE167118, randomly sampling cells from patient 17.
#'
#' @references
#' Unterman A et al. (2022) Single-cell multi-omics reveals dyssynchrony of
#' the innate and adaptive immune system in progressive COVID-19.
#' Nat Commun 13(1):440
#'
#' @name trex_example
#' @docType data
#' @keywords datasets
NULL
