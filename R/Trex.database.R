#' TCR-epitope annotation database
#'
#' Curated database of TCR CDR3 sequences and their associated epitopes,
#' compiled from multiple public repositories. Used by \code{\link{annotateDB}}
#' to annotate single-cell TCR data.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{TRA}{Data frame of TRA (alpha chain) CDR3 sequences and annotations}
#'   \item{TRB}{Data frame of TRB (beta chain) CDR3 sequences and annotations}
#' }
#'
#' Each data frame contains:
#' \describe{
#'   \item{CDR3}{CDR3 amino acid sequence}
#'   \item{Epitope.target}{Target antigen name}
#'   \item{Epitope.sequence}{Epitope peptide sequence}
#'   \item{Epitope.species}{Source organism}
#'   \item{Tissue}{Tissue of origin}
#'   \item{Cell.type}{Cell type annotation}
#'   \item{Database}{Source database}
#' }
#'
#' @source
#' Data compiled from:
#' \itemize{
#'   \item VDJdb (\url{https://vdjdb.cdr3.net/})
#'   \item McPAS-TCR (\url{http://friedmanlab.weizmann.ac.il/McPAS-TCR/})
#'   \item IEDB (\url{https://www.iedb.org/})
#'   \item PIRD (\url{https://db.cngb.org/pird/home/})
#' }
#'
#' @references
#' VDJdb: Shugay et al. (2018) Nucleic Acids Res 46(D1):D419-D427
#'
#' McPAS-TCR: Tickotsky et al. (2017) Bioinformatics 33(18):2924-2929
#'
#' IEDB: Vita et al. (2019) Nucleic Acids Res 47(D1):D339-D343
#'
#' PIRD: Zhang et al. (2020) Nucleic Acids Res 48(D1):D723-D730
#'
#' @name Trex.database
#' @docType data
#' @keywords datasets
NULL
