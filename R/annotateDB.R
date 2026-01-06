#' Annotate TCR sequences using reference database
#'
#' Use a curated database of TCR-epitope associations to annotate CDR3
#' sequences in single-cell data. The database includes entries from
#' VDJdb, McPAS-TCR, IEDB, and PIRD.
#'
#' @param sc Single-cell object in Seurat or SingleCellExperiment format.
#'   Must contain TCR information added via scRepertoire.
#' @param chains Character. TCR chain to annotate:
#'   \itemize{
#'     \item \code{"TRA"}: T cell receptor alpha chain
#'     \item \code{"TRB"}: T cell receptor beta chain
#'   }
#' @param edit.distance Integer. Maximum Levenshtein edit distance for matching
#'   CDR3 sequences to the reference database. Default is \code{0} (exact
#'   match only). Values of 1-2 are recommended when allowing mismatches.
#'
#' @return Seurat or SingleCellExperiment object with epitope annotation
#'   columns added to the metadata:
#'   \itemize{
#'     \item \code{<chain>_Epitope.target}: Target antigen name
#'     \item \code{<chain>_Epitope.sequence}: Epitope peptide sequence
#'     \item \code{<chain>_Epitope.species}: Source organism
#'     \item \code{<chain>_Tissue}: Tissue of origin
#'     \item \code{<chain>_Cell.type}: Cell type annotation
#'     \item \code{<chain>_Database}: Source database
#'   }
#'
#' @seealso \code{\link{runTrex}} for TCR embedding;
#'   \code{\link{Trex.database}} for database details
#'
#' @examples
#' trex_example <- annotateDB(trex_example,
#'     chains = "TRB"
#' )
#'
#' @export
#' @importFrom stringdist stringdist
annotateDB <- function(sc,
                       chains = "TRB",
                       edit.distance = 0) {
    TCR <- getTCR(sc, chains)
    relevent.db <- as.data.frame(Trex.database[[chains]])
    TCR <- TCR[[1]]

    cdr3.match <- TCR[TCR$cdr3_aa %in% relevent.db[, 1], ]
    cdr3.match <- merge(
        cdr3.match,
        relevent.db,
        by.x = "cdr3_aa",
        by.y = "CDR3",
        all.x = TRUE
    )

    barcodes <- cdr3.match$barcode
    cdr3.match <- cdr3.match[, c(
        "Epitope.target", "Epitope.sequence",
        "Epitope.species", "Tissue", "Cell.type", "Database"
    )]
    colnames(cdr3.match) <- paste0(chains, "_", colnames(cdr3.match))
    rownames(cdr3.match) <- barcodes

    if (edit.distance > 0) {
        annotated.cdr3 <- unique(TCR[TCR$barcode %in% barcodes, 2:5])[, 1]

        additions <- lapply(annotated.cdr3, function(x) {
            pos <- which(
                stringdist(x, TCR$cdr3_aa, method = "lv") <= edit.distance &
                    stringdist(x, TCR$cdr3_aa, method = "lv") > 0
            )
            pos
        })

        for (i in seq_along(additions)) {
            if (length(additions[[i]]) == 0) {
                next()
            }
            new.annotation <- cdr3.match[
                rownames(cdr3.match) %in%
                    TCR[TCR$cdr3_aa %in% annotated.cdr3[i], ]$barcode[1],
            ]
            tmp <- as.data.frame(t(matrix(
                nrow = ncol(new.annotation),
                ncol = length(TCR[additions[[i]], ]$barcode),
                data = unlist(new.annotation)
            )))

            rownames(tmp) <- TCR[additions[[i]], ]$barcode
            colnames(tmp) <- colnames(cdr3.match)
            cdr3.match <- rbind(cdr3.match, tmp)
        }
    }

    sc <- add.meta.data(sc, cdr3.match)
    return(sc)
}
