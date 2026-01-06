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
#' @seealso \code{\link{Trex_matrix}} for the full-featured matrix interface;
#'   \code{\link{maTrex}} for generating embedding matrix only;
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
                    geometric.method = "BLOSUM62",
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
        geometric.method = geometric.method,
        theta = theta
    )
    sc <- adding.DR(sc, reduction, reduction.name)
    return(sc)
}
