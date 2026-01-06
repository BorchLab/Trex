#' Remove TCR genes from variable features
#'
#' Most single-cell workflows use highly-expressed and highly-variable
#' genes for initial PCA and subsequent dimensional reduction. This function
#' removes T cell receptor genes from the variable features in a Seurat object
#' or from a character vector of gene names.
#'
#' Removing TCR genes before dimensional reduction prevents the embedding
#' from being dominated by clonotype information when the goal is to
#' capture transcriptional state. This is especially important when
#' combining TCR-based embeddings (from Trex) with RNA-based embeddings
#' using weighted nearest neighbors (WNN) approaches.
#'
#' @param sc Seurat object or character vector of variable gene names.
#' @param assay Character. For Seurat objects, the assay slot from which to
#'   remove TCR genes. Default (\code{NULL}) uses the default assay.
#'
#' @return Seurat object with TCR genes removed from variable features,
#'   or a character vector with TCR genes filtered out.
#'
#' @seealso \code{\link{runTrex}} for TCR embedding
#'
#' @examples
#' trex_example <- quietTCRgenes(trex_example)
#'
#' @export
#' @importFrom SeuratObject DefaultAssay VariableFeatures
#' @importFrom utils getFromNamespace
#'
#' @author Nicky de Vrij, Nikolaj Pagh, Nick Borcherding
quietTCRgenes <- function(sc,
                          assay = NULL) {
    `VariableFeatures<-` <- utils::getFromNamespace(
        "VariableFeatures<-",
        "SeuratObject"
    )
    unwanted_genes <- "^TR[ABDG][VDJ][^D]"

    if (inherits(x = sc, what = "Seurat")) {
        if (is.null(assay)) {
            assay <- DefaultAssay(sc)
        }
        unwanted_genes <- grep(
            pattern = unwanted_genes,
            x = VariableFeatures(sc, assay = assay),
            value = TRUE
        )
        VariableFeatures(sc, assay = assay) <- VariableFeatures(
            sc,
            assay = assay
        )[VariableFeatures(sc, assay = assay) %!in% unwanted_genes]
    } else {
        # Character vector input (e.g., from Bioconductor scran workflows)
        unwanted_genes <- grep(
            pattern = unwanted_genes,
            x = sc,
            value = TRUE
        )
        sc <- sc[sc %!in% unwanted_genes]
    }
    return(sc)
}
