#' Generate TCR embedding matrix
#'
#' Use this function to run the Trex algorithm and return latent vectors
#' representing TCR CDR3 sequences as a matrix. This is the core function
#' of the Trex algorithm and can be used independently of single-cell objects.
#'
#' @param sc Single-cell object (Seurat or SingleCellExperiment) or the output
#'   of \code{\link[scRepertoire]{combineTCR}} from the scRepertoire package.
#' @param chains Character. TCR chain to use for encoding:
#'   \itemize{
#'     \item \code{"TRA"}: T cell receptor alpha chain
#'     \item \code{"TRB"}: T cell receptor beta chain
#'   }
#' @param method Character. Encoding method to use:
#'   \itemize{
#'     \item \code{"encoder"}: Deep learning autoencoders for sequence embedding
#'     \item \code{"geometric"}: Geometric transformation based on BLOSUM62
#'       substitution matrix
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
#'     \item \code{"AF"}: Atchley factors (5 physicochemical properties)
#'     \item \code{"KF"}: Kidera factors (10 physicochemical properties)
#'     \item \code{"OHE"}: One-hot encoding (21 dimensions per position)
#'   }
#' @param theta Numeric. Rotation angle in radians for geometric transformation
#'   when \code{method = "geometric"}. Default is \code{pi}.
#'
#' @return A data frame of Trex-encoded values with cells as rows and
#'   embedding dimensions as columns (30 dimensions for encoder method,
#'   20 dimensions for geometric method).
#'
#' @seealso \code{\link{runTrex}} for integration with single-cell objects;
#'   \code{\link{CoNGAfy}} for clonotype-based reduction
#'
#' @examples
#' Trex_values <- maTrex(trex_example,
#'     chains = "TRA",
#'     method = "encoder",
#'     encoder.model = "VAE",
#'     encoder.input = "AF"
#' )
#'
#' Trex_values <- maTrex(trex_example,
#'     chains = "TRA",
#'     method = "geometric",
#'     theta = pi
#' )
#'
#' @export
#' @importFrom SeuratObject CreateDimReducObject
maTrex <- function(sc,
                   chains = "TRA",
                   method = "encoder",
                   encoder.model = "VAE",
                   encoder.input = "AF",
                   theta = pi) {
    TCR <- getTCR(sc, chains)
    checkLength(TCR[[1]])

    if (method == "encoder" &&
        encoder.input %in% c("AF", "KF", "both", "all", "OHE")) {
        message("Calculating the encoding values...")
        reduction <- .encoder(TCR, encoder.input, encoder.model)
    } else if (method == "geometric") {
        message("Performing geometric transformation...")
        reduction <- .geometric.encoding(TCR, theta)
    }
    return(reduction)
}

#' Run Trex embedding on single-cell object
#'
#' Run the Trex algorithm and integrate the resulting embedding into a
#' Seurat or SingleCellExperiment object. This function will automatically
#' filter cells that do not have TCR information in the metadata.
#'
#' @inheritParams maTrex
#' @param sc Single-cell object in Seurat or SingleCellExperiment format.
#'   Must contain TCR information added via scRepertoire.
#' @param reduction.name Character. Name for the dimensional reduction slot.
#'   Useful when running Trex with multiple parameter combinations.
#'   Default is \code{"Trex"}.
#'
#' @return Seurat or SingleCellExperiment object with Trex dimensions placed
#'   into the dimensional reduction slot.
#'
#' @seealso \code{\link{maTrex}} for generating embedding matrix only;
#'   \code{\link{quietTCRgenes}} for removing TCR genes before WNN analysis
#'
#' @examples
#' trex_example <- runTrex(trex_example,
#'     chains = "TRA",
#'     method = "encoder",
#'     encoder.model = "VAE",
#'     encoder.input = "AF"
#' )
#'
#' @export
runTrex <- function(sc,
                    chains = "TRA",
                    method = "encoder",
                    encoder.model = "VAE",
                    encoder.input = "AF",
                    theta = pi,
                    reduction.name = "Trex") {
    checkSingleObject(sc)
    sc <- filter.object(sc)
    reduction <- maTrex(
        sc = sc,
        chains = chains,
        method = method,
        encoder.model = encoder.model,
        encoder.input = encoder.input,
        theta = theta
    )
    sc <- adding.DR(sc, reduction, reduction.name)
    return(sc)
}
