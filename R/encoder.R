#' Encode TCR CDR3 sequences using deep learning models
#'
#' Internal function that transforms CDR3 amino acid sequences into
#' fixed-dimensional vectors using pre-trained autoencoder models.
#' Uses immApex for sequence encoding before passing to the model.
#'
#' @param TCR List containing TCR data with barcode and cdr3_aa columns
#' @param encoder.input Encoding method: "AF" (Atchley factors),
#'   "KF" (Kidera factors), or "OHE" (one-hot encoding)
#' @param encoder.model Model type: "AE" (autoencoder) or "VAE"
#'   (variational autoencoder)
#'
#' @return Data frame with encoded vectors (30 dimensions for AE/VAE)
#' @keywords internal
#' @importFrom reticulate array_reshape
#' @importFrom basilisk basiliskRun
#' @importFrom immApex onehotEncoder propertyEncoder
.encoder <- function(TCR,
                     encoder.input = encoder.input,
                     encoder.model = encoder.model) {
    chain <- names(TCR)
    aa.model <- quiet(aa.model.loader(chain, encoder.input, encoder.model))
    membership <- TCR[[1]]
    cells <- unique(membership[, "barcode"])

    # Map encoder.input to immApex property.set
    property.set <- switch(encoder.input,
        "AF" = "atchleyFactors",
        "KF" = "kideraFactors",
        NULL
    )

    score <- NULL
    for (n in seq_len(length(cells))) {
        tmp.CDR <- membership[membership$barcode == cells[n], ]$cdr3_aa

        if (encoder.input == "OHE") {
            # Use immApex onehotEncoder
            encoded <- onehotEncoder(
                input.sequences = tmp.CDR,
                max.length = 60,
                verbose = FALSE
            )
            array.reshape.tmp <- array_reshape(encoded$flattened, 1260)
        } else {
            # Use immApex propertyEncoder for AF or KF
            encoded <- propertyEncoder(
                input.sequences = tmp.CDR,
                property.set = property.set,
                max.length = 60,
                verbose = FALSE
            )
            array.reshape.tmp <- array_reshape(
                encoded$flattened,
                ncol(encoded$flattened)
            )
        }

        score.tmp <- auto.embedder(array.reshape.tmp, aa.model, encoder.input)
        score <- rbind(score, score.tmp)
    }

    score <- data.frame(unique(membership[, "barcode"]), score)
    barcodes <- score[, 1]
    score <- score[, -1]
    rownames(score) <- barcodes
    colnames(score) <- paste0("Trex_", seq_len(ncol(score)))
    return(score)
}
