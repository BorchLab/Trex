#' Run Trex embedding on single-cell object
#'
#' Run the Trex algorithm and integrate the resulting embedding into a
#' Seurat or SingleCellExperiment object. This function will automatically
#' filter cells that do not have TCR information in the metadata.
#'
#' @param sc Single-cell object in Seurat or SingleCellExperiment format.
#'   Must contain TCR information added via scRepertoire.
#' @param chain Character. TCR chain to use for encoding:
#'   \itemize{
#'     \item \code{"TRA"}: T cell receptor alpha chain
#'     \item \code{"TRB"}: T cell receptor beta chain
#'   }
#' @param method Character. Encoding method to use:
#'   \itemize{
#'     \item \code{"encoder"}: Deep learning autoencoders for sequence embedding
#'     \item \code{"geometric"}: Geometric transformation based on substitution
#'       matrices
#'   }
#' @param encoder.model Character. Model architecture when
#'   \code{method = "encoder"}:
#'   \itemize{
#'     \item \code{"VAE"}: Variational autoencoder
#'     \item \code{"AE"}: Dense autoencoder
#'   }
#' @param encoder.input Character. Input encoding when
#'   \code{method = "encoder"}:
#'   \itemize{
#'     \item \code{"atchleyFactors"}: Atchley factors (5 physicochemical
#'       properties)
#'     \item \code{"kideraFactors"}: Kidera factors (10 physicochemical
#'       properties)
#'     \item \code{"OHE"}: One-hot encoding (21 dimensions per position)
#'   }
#' @param geometric.method Character or matrix. The substitution matrix to use
#'   for geometric encoding (e.g., \code{"BLOSUM62"}, \code{"BLOSUM45"},
#'   \code{"PAM250"}), or a custom 20x20 numeric matrix.
#'   Default is \code{"BLOSUM62"}.
#' @param geometric.theta Numeric. Rotation angle in radians for geometric
#'   transformation when \code{method = "geometric"}. Default is \code{pi}.
#' @param reduction.name Character. Name for the dimensional reduction slot.
#'   Useful when running Trex with multiple parameter combinations.
#'   Default is \code{"Trex"}.
#'
#' @return Seurat or SingleCellExperiment object with Trex dimensions placed
#'   into the dimensional reduction slot.
#'
#' @seealso \code{\link{Trex_matrix}} for the matrix-only interface;
#'   \code{\link{quietTCRgenes}} for removing TCR genes before WNN analysis
#'
#' @examples
#' trex_example <- runTrex(trex_example,
#'     chain = "TRA",
#'     method = "encoder",
#'     encoder.model = "VAE",
#'     encoder.input = "atchleyFactors"
#' )
#'
#' @export
runTrex <- function(sc,
                    chain = "TRA",
                    method = "encoder",
                    encoder.model = "VAE",
                    encoder.input = "atchleyFactors",
                    geometric.method = "BLOSUM62",
                    geometric.theta = pi,
                    reduction.name = "Trex") {
    checkSingleObject(sc)
    sc <- filter.object(sc)
    reduction <- Trex_matrix(
        input.data = sc,
        chain = chain,
        method = method,
        encoder.model = encoder.model,
        encoder.input = encoder.input,
        geometric.method = geometric.method,
        geometric.theta = geometric.theta,
        verbose = TRUE
    )
    sc <- adding.DR(sc, reduction, reduction.name)
    return(sc)
}
