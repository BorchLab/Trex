# test script for annotateDB.R - testcases are NOT comprehensive!

test_that("annotateDB works with Seurat objects", {
    data("trex_example")

    # Test default (TRB chain)
    trex_example <- annotateDB(trex_example)

    expect_equal(
        trex_example@meta.data,
        getdata("annotateDB", "annotateDB_TRB_meta.data")
    )

    # Test TRA chain with edit distance
    trex_example <- annotateDB(trex_example,
        chains = "TRA",
        edit.distance = 1
    )

    expect_equal(
        trex_example@meta.data,
        getdata("annotateDB", "annotateDB_TRA_meta.data")
    )
})

test_that("annotateDB works with SingleCellExperiment objects", {
    skip_if_not_installed("SingleCellExperiment")

    sce <- getdata("runTrex", "SCE.object")

    sce <- annotateDB(sce, chains = "TRA")

    expect_s4_class(sce, "SingleCellExperiment")
    expect_true("annotateDB.TRA" %in% colnames(SummarizedExperiment::colData(sce)))
})