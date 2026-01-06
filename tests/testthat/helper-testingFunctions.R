#' Load test data from the testdata directory
#'
#' @param dir Subdirectory within testdata
#' @param name Name of the RDS file (without extension)
#' @return The loaded R object
#' @keywords internal
getdata <- function(dir, name) {
    path <- file.path("testdata", dir, paste0(name, ".rds"))
    if (!file.exists(path)) {
        stop("Test data file not found: ", path)
    }
    readRDS(path)
}