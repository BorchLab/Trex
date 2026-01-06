# test script for quietTCRgenes.R - testcases are NOT comprehensive!

test_that("quietTCRgenes works with character vector", {
    data("trex_example")

    features <- rownames(trex_example[["RNA"]]$counts)

    expect_equal(
        quietTCRgenes(features),
        getdata("quietTCRgenes", "quietTCRgenes_feature.vector")
    )
})

test_that("quietTCRgenes works with Seurat object", {
    data("trex_example")

    features <- rownames(trex_example[["RNA"]]$counts)
    SeuratObject::VariableFeatures(trex_example) <- features

    trex_example <- quietTCRgenes(trex_example)

    expect_equal(
        quietTCRgenes(features),
        SeuratObject::VariableFeatures(trex_example)
    )
})

test_that("quietTCRgenes removes TCR genes", {
    # Test that TCR genes are properly removed
    test_genes <- c("TRAV1-1", "TRBV2-1", "CD4", "CD8A", "TRDV1", "HLA-A")

    result <- quietTCRgenes(test_genes)

    expect_false("TRAV1-1" %in% result)
    expect_false("TRBV2-1" %in% result)
    expect_false("TRDV1" %in% result)
    expect_true("CD4" %in% result)
    expect_true("CD8A" %in% result)
    expect_true("HLA-A" %in% result)
})