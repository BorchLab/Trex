#' @keywords internal
.onLoad <- function(libname, pkgname) {
    # Set up cache directory for Trex
    cache_dir <- tools::R_user_dir("Trex", which = "cache")

    # Create basilisk subdirectory for Python environment
    basilisk_dir <- file.path(cache_dir, "basilisk")
    dir.create(basilisk_dir, recursive = TRUE, showWarnings = FALSE)

    # Set basilisk options
    options(
        basilisk.external.dir = basilisk_dir,
        dir.expiry.dir = basilisk_dir
    )

    # Register global variables to avoid R CMD check NOTEs
    utils::globalVariables(c(
        "array_reshape",
        "is",
        "reducedDim<-",
        "na.omit",
        "median",
        "slot",
        "get.adjacency",
        "nn",
        "data",
        "Trex.database",
        "colData<-",
        "TR",
        "graph.edgelist",
        "f"
    ))

    invisible()
}
