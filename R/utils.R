"%!in%" <- Negate("%in%")

#' @keywords internal
aa.eval <- function(x) {
    x %in% c("AF", "KF", "other")
}

#' Add metadata to single-cell object
#'
#' @param sc Single-cell object (Seurat or SingleCellExperiment)
#' @param meta Metadata to add
#' @param header Column names for metadata
#'
#' @return Updated single-cell object
#' @keywords internal
#' @importFrom rlang %||%
#' @importFrom SingleCellExperiment colData
add.meta.data <- function(sc, meta, header) {
    if (inherits(x = sc, what = "Seurat")) {
        col.name <- names(meta) %||% colnames(meta)
        sc[[col.name]] <- meta
    } else {
        rownames <- rownames(colData(sc))
        colData(sc) <- cbind(
            colData(sc),
            meta[rownames, ]
        )[, union(colnames(colData(sc)), colnames(meta))]
        rownames(colData(sc)) <- rownames
    }
    return(sc)
}

#' Extract TCR information from single-cell object
#'
#' @param sc Single-cell object or scRepertoire output
#' @param chains TCR chain(s) to extract
#'
#' @return List of TCR data frames by chain
#' @keywords internal
#' @importFrom stringr str_split
getTCR <- function(sc, chains) {
    if (inherits(x = sc, what = "Seurat") ||
        inherits(x = sc, what = "SingleCellExperiment")) {
        meta <- grabMeta(sc)
    } else {
        meta <- do.call(rbind, sc)
        rownames(meta) <- meta[, "barcode"]
    }
    tmp <- data.frame(
        barcode = rownames(meta),
        str_split(meta[, "CTaa"], "_", simplify = TRUE),
        str_split(meta[, "CTgene"], "_", simplify = TRUE)
    )
    if (length(chains) == 1 && chains != "both") {
        if (chains %in% c("TRA")) {
            pos <- list(c(2, 4))
        } else if (chains %in% c("TRB")) {
            pos <- list(c(3, 5))
        }
    } else {
        pos <- list(one = c(2, 4), two = c(3, 5))
        chains <- c("TRA", "TRB")
    }
    TCR <- NULL
    for (i in seq_along(pos)) {
        sub <- as.data.frame(tmp[, c(1, pos[[i]])])
        colnames(sub) <- c("barcode", "cdr3_aa", "genes")
        sub$v <- str_split(sub$genes, "[.]", simplify = TRUE)[, 1]
        sub$j <- str_split(sub$genes, "[.]", simplify = TRUE)[, 2]
        sub[sub == ""] <- NA
        TCR[[i]] <- sub
        sub <- NULL
    }
    names(TCR) <- chains
    return(TCR)
}

#' Extract metadata from single-cell object
#'
#' @param sc Single-cell object (Seurat or SingleCellExperiment)
#'
#' @return Data frame of metadata
#' @keywords internal
#' @importFrom SingleCellExperiment colData
grabMeta <- function(sc) {
    if (inherits(x = sc, what = "Seurat")) {
        meta <- data.frame(sc[[]], slot(sc, "active.ident"))
        if ("cluster" %in% colnames(meta)) {
            colnames(meta)[length(meta)] <- "cluster.active.ident"
        } else {
            colnames(meta)[length(meta)] <- "cluster"
        }
    } else if (inherits(x = sc, what = "SingleCellExperiment")) {
        meta <- data.frame(colData(sc))
        rownames(meta) <- sc@colData@rownames
        clu <- which(colnames(meta) == "ident")
        if (length(clu) > 0) {
            if ("cluster" %in% colnames(meta)) {
                colnames(meta)[clu] <- "cluster.active.idents"
            } else {
                colnames(meta)[clu] <- "cluster"
            }
        }
    }
    return(meta)
}

#' Suppress output from expressions
#'
#' @param x Expression to evaluate quietly
#'
#' @return Result of expression (invisibly)
#' @keywords internal
quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}

#' Validate single-cell object type
#'
#' @param sc Object to validate
#'
#' @return NULL (stops with error if invalid)
#' @keywords internal
checkSingleObject <- function(sc) {
    if (!inherits(x = sc, what = "Seurat") &&
        !inherits(x = sc, what = "SummarizedExperiment")) {
        stop(
            "Object indicated is not of class 'Seurat' or ",
            "'SummarizedExperiment', make sure you are using ",
            "the correct data."
        )
    }
}

#' Check CDR3 sequence length
#'
#' @param x Data frame with cdr3_aa column
#'
#' @return NULL (stops with error if sequences too long)
#' @keywords internal
checkLength <- function(x) {
    if (any(na.omit(nchar(x[, "cdr3_aa"])) > 60)) {
        stop(
            "Models have been trained on CDR3 sequences ",
            "less than 60 amino acid residues. Please ",
            "filter the larger sequences before running."
        )
    }
}

#' Load pre-trained autoencoder model from Zenodo
#'
#' Downloads and caches model files from Zenodo repository.
#'
#' @param chain TCR chain ("TRA" or "TRB")
#' @param encoder.input Input encoding type ("AF", "KF", or "OHE")
#' @param encoder.model Model type ("AE" or "VAE")
#'
#' @return Path to the cached model file
#' @keywords internal
aa.model.loader <- function(chain, encoder.input, encoder.model) {
    # Construct model filename
    model_name <- paste0(chain, "_", encoder.input, "_", encoder.model, ".h5")

    # Check metadata for valid model names
    valid_models <- c(
        "TRA_AF_AE.h5", "TRA_AF_VAE.h5",
        "TRA_KF_AE.h5", "TRA_KF_VAE.h5",
        "TRA_OHE_AE.h5", "TRA_OHE_VAE.h5",
        "TRB_AF_AE.h5", "TRB_AF_VAE.h5",
        "TRB_KF_AE.h5", "TRB_KF_VAE.h5",
        "TRB_OHE_AE.h5", "TRB_OHE_VAE.h5"
    )

    if (!model_name %in% valid_models) {
        stop(
            "Invalid model specification. Please check chain, ",
            "encoder.input, and encoder.model parameters."
        )
    }

    # Set up cache directory
    cache_dir <- tools::R_user_dir("Trex", which = "cache")
    models_dir <- file.path(cache_dir, "models")
    dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

    model_path <- file.path(models_dir, model_name)

    # Download if not cached
    if (!file.exists(model_path)) {
        base_url <- "https://zenodo.org/records/14636032/files"
        url <- paste0(base_url, "/", model_name)

        message("Downloading model: ", model_name)
        tryCatch(
            {
                download.file(
                    url = url,
                    destfile = model_path,
                    mode = "wb",
                    quiet = FALSE
                )
            },
            error = function(e) {
                stop(
                    "Failed to download model from Zenodo. ",
                    "Error: ", conditionMessage(e)
                )
            }
        )
    }

    return(normalizePath(model_path))
}

#' Add dimensional reduction to single-cell object
#'
#' @param sc Single-cell object
#' @param reduction Reduction matrix
#' @param reduction.name Name for the reduction
#'
#' @return Updated single-cell object
#' @keywords internal
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim<- reducedDim
adding.DR <- function(sc, reduction, reduction.name) {
    if (inherits(sc, "Seurat")) {
        DR <- suppressWarnings(CreateDimReducObject(
            embeddings = as.matrix(reduction),
            loadings = as.matrix(reduction),
            projected = as.matrix(reduction),
            stdev = rep(0, ncol(reduction)),
            key = reduction.name,
            jackstraw = NULL,
            misc = list()
        ))
        sc[[reduction.name]] <- DR
    } else if (inherits(sc, "SingleCellExperiment")) {
        reducedDim(sc, reduction.name) <- reduction
    }
    return(sc)
}

#' Filter object to cells with TCR data
#'
#' @param sc Single-cell object
#'
#' @return Filtered single-cell object
#' @keywords internal
filter.object <- function(sc) {
    meta <- grabMeta(sc)
    cells.chains <- rownames(meta[!is.na(meta[["CTaa"]]), ])
    if (inherits(x = sc, what = "Seurat")) {
        sc <- subset(sc, cells = cells.chains)
    } else if (inherits(x = sc, what = "SingleCellExperiment")) {
        sc <- sc[, colnames(sc) %in% cells.chains]
    }
    return(sc)
}
