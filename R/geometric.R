#' Geometric encoding of TCR sequences
#'
#' Encode CDR3 sequences using BLOSUM62 matrix with geometric rotation
#' transformation. This provides an alternative to deep learning-based
#' encoding that does not require Python dependencies. Uses immApex's
#' geometricEncoder internally.
#'
#' @param TCR List containing TCR data with barcode and cdr3_aa columns
#' @param theta Rotation angle in radians (default pi)
#'
#' @return Data frame with 20-dimensional geometric embeddings
#' @keywords internal
#' @importFrom immApex geometricEncoder
.geometric.encoding <- function(TCR,
                                theta = pi) {
    membership <- TCR[[1]]
    cells <- unique(membership[, "barcode"])
    sequences <- membership$cdr3_aa[match(cells, membership$barcode)]

    # Handle NA sequences
    valid_idx <- !is.na(sequences)
    valid_sequences <- sequences[valid_idx]

    if (length(valid_sequences) > 0) {
        # Use immApex geometricEncoder with BLOSUM62
        encoded <- geometricEncoder(
            input.sequences = valid_sequences,
            method = "BLOSUM62",
            theta = theta,
            verbose = FALSE
        )
        score <- encoded$summary
    } else {
        score <- matrix(0, nrow = 0, ncol = 20)
    }

    # Create full result matrix with zeros for NA sequences
    full_score <- matrix(0, nrow = length(cells), ncol = 20)
    full_score[valid_idx, ] <- score

    rownames(full_score) <- cells
    colnames(full_score) <- paste0("Trex_", seq_len(ncol(full_score)))
    return(full_score)
}
