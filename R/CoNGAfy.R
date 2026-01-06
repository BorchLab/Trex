#' Reduce single-cell object to representative cells by clonotype
#'
#' Generate a single-cell object that has a single representative cell
#' for each clonotype. This approach was first introduced in
#' \href{https://pubmed.ncbi.nlm.nih.gov/34426704/}{CoNGA} and was
#' adapted to R. Please read and cite the original work.
#'
#' @param sc Single-cell object in Seurat or SingleCellExperiment format.
#'   Must contain TCR information added via scRepertoire.
#' @param method Character. Method for selecting representative cells:
#'   \itemize{
#'     \item \code{"mean"}: Generate mean expression values across all cells
#'       sharing the same clonotype
#'     \item \code{"dist"}: Select the cell with minimal Euclidean distance
#'       from the clonotype centroid in PCA space
#'   }
#' @param features Character vector. Genes to include in the output.
#'   Default (\code{NULL}) returns all genes.
#' @param assay Character. Name of the assay(s) to return. Default is
#'   \code{"RNA"}. Can be a vector for multiple assays.
#' @param meta.carry Character vector. Metadata variables to carry over from
#'   the original single-cell object. Default is \code{c("CTaa", "CTgene")}.
#'
#' @return Single-cell object with one cell representing each unique clonotype.
#'
#' @seealso \code{\link{runTrex}} for embedding the reduced object;
#'   \code{\link{annotateDB}} for epitope annotation
#'
#' @examples
#' CoNGA.seurat <- CoNGAfy(trex_example,
#'     method = "dist",
#'     features = NULL
#' )
#'
#' @export
#' @importFrom SeuratObject CreateSeuratObject CreateAssayObject
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay assay<- colData<-
CoNGAfy <- function(sc,
                    method = "dist",
                    features = NULL,
                    assay = "RNA",
                    meta.carry = c("CTaa", "CTgene")) {
    sc <- filter.object(sc)
    conga <- NULL

    if (method == "mean") {
        for (x in seq_along(assay)) {
            conga[[x]] <- CoNGA.mean(sc, features, assay[x])
        }
    } else if (method == "dist") {
        for (x in seq_along(assay)) {
            conga[[x]] <- CoNGA.dist(sc, features, assay[x])
        }
    }

    names(conga) <- assay

    if (inherits(x = sc, what = "Seurat")) {
        sc.output <- CreateSeuratObject(
            conga[[1]],
            assay = names(conga)[1],
            project = "Trex"
        )
        if (length(conga) > 1) {
            for (y in 2:length(conga)) {
                sc.output[[names(conga)[y]]] <- CreateAssayObject(conga[[y]])
            }
        }
        CTge <- unique(sc[[]][, c(meta.carry)])
    } else if (inherits(x = sc, what = "SingleCellExperiment")) {
        sc.output <- SingleCellExperiment(assay = conga[[1]])
        if (length(conga) > 1) {
            for (y in 2:length(conga)) {
                assay(sc.output, names(conga)[y]) <- conga[[y]]
            }
        }
        sc.output$CTaa <- rownames(sc.output@colData)
        CTge <- data.frame(unique(sc@colData[, c(meta.carry)]))
    }

    CTge <- meta.handler(CTge, meta.carry)
    clones <- unique(CTge$CTaa)
    rownames(CTge) <- clones
    colnames(CTge) <- meta.carry
    sc.output <- add.meta.data(sc.output, CTge, colnames(CTge))
    return(sc.output)
}

#' Select representative cells by distance
#'
#' For single clones, uses actual RNA values. For multiplets, selects
#' the cell with minimal Euclidean distance in PCA space. For doublets,
#' automatically selects the first cell.
#'
#' @param sc Single-cell object
#' @param features Features to include
#' @param assay Assay name
#'
#' @return Sparse matrix of counts for representative cells
#' @keywords internal
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment reducedDim
CoNGA.dist <- function(sc, features, assay) {
    if (inherits(x = sc, what = "Seurat")) {
        if (assay == "RNA") {
            data.use <- sc[["pca"]]@cell.embeddings
        } else if (assay == "ADT") {
            data.use <- sc[["apca"]]@cell.embeddings
        }
    } else if (inherits(x = sc, what = "SingleCellExperiment")) {
        data.use <- reducedDim(sc, "PCA")
    }

    meta <- grabMeta(sc)
    data <- as.data.frame(meta[, "CTaa"])
    colnames(data) <- "CTaa"
    rownames(data) <- rownames(meta)

    all.clones <- table(data$CTaa)
    unique.clones <- all.clones[which(all.clones > 1)]
    single.clones <- all.clones[all.clones %!in% unique.clones]
    barcodes <- rownames(data)[which(data$CTaa %in% names(single.clones))]

    for (i in seq_along(unique.clones)) {
        loc <- which(data$CTaa == names(unique.clones)[i])
        dist <- as.matrix(dist(data.use[loc, ]))
        cell <- names(which(rowSums(dist) == min(rowSums(dist))))
        if (length(cell) > 1) {
            cell <- cell[1]
        }
        barcodes <- c(barcodes, cell)
    }

    if (inherits(x = sc, what = "Seurat")) {
        assay.use <- sc[[assay]]$counts
    } else if (inherits(x = sc, what = "SingleCellExperiment")) {
        assay.use <- assay(sc)
    }

    features.to.avg <- features %||% rownames(x = assay.use)
    features.assay <- intersect(x = features.to.avg, y = rownames(x = assay.use))
    data.return <- assay.use[
        rownames(assay.use) %in% features.assay,
        colnames(assay.use) %in% barcodes
    ]
    colnames(data.return) <- data$CTaa[match(barcodes, rownames(data))]
    return(data.return)
}

#' Calculate mean expression by clonotype
#'
#' Adapted from the AverageExpression() function in Seurat.
#'
#' @param sc Single-cell object
#' @param features Features to include
#' @param assay Assay name
#'
#' @return Matrix of mean expression values by clonotype
#' @keywords internal
#' @importFrom rlang %||%
#' @importFrom Matrix sparse.model.matrix colSums
#' @importFrom SummarizedExperiment assay
#' @importFrom stats as.formula
CoNGA.mean <- function(sc, features, assay) {
    if (inherits(x = sc, what = "Seurat")) {
        data.use <- sc[[assay]]$counts
        data.use <- expm1(x = data.use)
    } else if (inherits(x = sc, what = "SingleCellExperiment")) {
        data.use <- assay(sc, name = assay)
    }

    features.to.avg <- features %||% rownames(x = data.use)
    features.assay <- intersect(x = features.to.avg, y = rownames(x = data.use))

    meta <- grabMeta(sc)
    data <- as.data.frame(meta[, "CTaa"])
    colnames(data) <- "CTaa"
    rownames(data) <- rownames(meta)

    # Filter cells without CTaa
    data.use <- data.use[, !is.na(data[, "CTaa"])]
    data <- data[which(rowSums(x = is.na(x = data)) == 0), , drop = FALSE]

    for (i in seq_len(ncol(x = data))) {
        data[, i] <- as.factor(x = data[, i])
    }

    num.levels <- vapply(
        X = seq_len(ncol(x = data)),
        FUN = function(i) {
            length(x = levels(x = data[, i]))
        },
        FUN.VALUE = integer(1)
    )

    category.matrix <- sparse.model.matrix(
        object = as.formula(
            object = paste0(
                "~0+",
                paste0(
                    "data[,",
                    seq_len(length(x = "CTaa")),
                    "]",
                    collapse = ":"
                )
            )
        )
    )

    colsums <- Matrix::colSums(x = category.matrix)
    category.matrix <- category.matrix[, colsums > 0]
    colsums <- colsums[colsums > 0]

    category.matrix <- sweep(
        x = category.matrix,
        MARGIN = 2,
        STATS = colsums,
        FUN = "/"
    )

    colnames(x = category.matrix) <- vapply(
        X = colnames(x = category.matrix),
        FUN = function(name) {
            name <- gsub(
                pattern = "data\\[, [1-9]*\\]",
                replacement = "",
                x = name
            )
            return(paste0(
                rev(x = unlist(x = strsplit(x = name, split = ":"))),
                collapse = "_"
            ))
        },
        FUN.VALUE = character(1)
    )

    data.return <- data.use %*% category.matrix
    return(data.return)
}

#' Handle metadata for clonotype consolidation
#'
#' For duplicated clonotypes with different metadata values,
#' concatenates unique values with semicolons.
#'
#' @param meta Metadata data frame
#' @param meta.carry Variables to process
#'
#' @return Processed metadata data frame
#' @keywords internal
#' @importFrom stringr str_sort
meta.handler <- function(meta, meta.carry) {
    unique.clones <- unique(meta[, "CTaa"])
    duplicated.clones <- na.omit(unique(
        meta[, "CTaa"][which(duplicated(meta[, "CTaa"]))]
    ))
    new.meta <- NULL

    for (i in seq_along(duplicated.clones)) {
        meta.tmp <- meta[meta[, "CTaa"] == duplicated.clones[i], ]
        concat.strings <- lapply(meta.carry, function(x) {
            paste0(str_sort(unique(na.omit(meta.tmp[, x]))), collapse = ";")
        })
        new.meta <- rbind(new.meta, unlist(concat.strings))
    }

    old.meta <- meta[meta[, "CTaa"] %!in% duplicated.clones, meta.carry]
    if (length(new.meta) != 0) {
        colnames(new.meta) <- meta.carry
        total.meta <- rbind(old.meta, new.meta)
    } else {
        total.meta <- old.meta
    }
    return(total.meta)
}
