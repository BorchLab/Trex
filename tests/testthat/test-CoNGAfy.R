# test script for CoNGAfy.R - testcases are NOT comprehensive!

test_that("CoNGAfy works with Seurat objects", {
    data("trex_example")

    # Test default method (dist)
    conga_reduction <- CoNGAfy(trex_example)

    expect_equal(
        conga_reduction@meta.data,
        getdata("CoNGAfy", "CoNGAfy_meta.data")
    )

    expect_equal(
        conga_reduction[["RNA"]]$counts,
        getdata("CoNGAfy", "CoNGAfy_counts"),
        tolerance = 1e-2
    )

    # Test mean method
    conga_mean_reduction <- CoNGAfy(trex_example,
        method = "mean"
    )

    expect_equal(
        conga_mean_reduction@meta.data,
        getdata("CoNGAfy", "CoNGAfy_mean_meta.data")
    )

    expect_equal(
        conga_mean_reduction[["RNA"]]$counts,
        getdata("CoNGAfy", "CoNGAfy_mean_counts"),
        tolerance = 1e-2
    )

    # Test with ADT assay
    conga_ADT_reduction <- CoNGAfy(trex_example,
        assay = "ADT"
    )

    expect_equal(
        conga_ADT_reduction[["ADT"]]$counts,
        getdata("CoNGAfy", "CoNGAfy_ADT_counts"),
        tolerance = 1e-2
    )
})

test_that("CoNGAfy works with SingleCellExperiment objects", {
    skip_if_not_installed("SingleCellExperiment")

    sce <- getdata("runTrex", "SCE.object")

    conga_sce <- CoNGAfy(sce)

    expect_s4_class(conga_sce, "SingleCellExperiment")
    expect_true("CTaa" %in% colnames(SummarizedExperiment::colData(conga_sce)))
})