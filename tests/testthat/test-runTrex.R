# test script for runTrex.R - testcases are NOT comprehensive!

test_that("runTrex works with encoder method on Seurat objects", {
    data("trex_example")
    set.seed(42)

    # Test TRA chain with VAE and atchleyFactors input
    trex_example <- runTrex(trex_example,
        chain = "TRA",
        method = "encoder",
        encoder.model = "VAE",
        encoder.input = "atchleyFactors",
        reduction.name = "TRA_VAE_AF"
    )

    expect_equal(
        trex_example@reductions$TRA_VAE_AF@cell.embeddings,
        getdata("runTrex", "runTrex_TRA_VAE_AF_reduction"),
        tolerance = 1e-2
    )

    # Test TRB chain with AE and kideraFactors input
    trex_example <- runTrex(trex_example,
        chain = "TRB",
        method = "encoder",
        encoder.model = "AE",
        encoder.input = "kideraFactors",
        reduction.name = "TRB_AE_KF"
    )

    expect_equal(
        trex_example@reductions$TRB_AE_KF@cell.embeddings,
        getdata("runTrex", "runTrex_TRB_AE_kF_reduction"),
        tolerance = 1e-2
    )

    # Test TRB chain with VAE and OHE input
    trex_example <- runTrex(trex_example,
        chain = "TRB",
        method = "encoder",
        encoder.model = "VAE",
        encoder.input = "OHE",
        reduction.name = "TRB_VAE_OHE"
    )

    expect_equal(
        trex_example@reductions$TRB_VAE_OHE@cell.embeddings,
        getdata("runTrex", "runTrex_TRB_VAE_OHE_reduction"),
        tolerance = 1e-2
    )
})

test_that("runTrex works with geometric method on Seurat objects", {
    data("trex_example")
    set.seed(42)

    trex_example <- runTrex(trex_example,
        chain = "TRB",
        method = "geometric",
        reduction.name = "TRB_Geometric"
    )

    expect_equal(
        trex_example@reductions$TRB_Geometric@cell.embeddings,
        getdata("runTrex", "runTrex_TRB_geometric_reduction"),
        tolerance = 1e-2
    )
})

test_that("runTrex works with SingleCellExperiment objects", {
    skip_if_not_installed("SingleCellExperiment")

    sce <- getdata("runTrex", "SCE.object")

    sce <- runTrex(sce,
        chain = "TRA",
        method = "geometric",
        reduction.name = "TRA_Geometric"
    )

    expect_equal(
        SingleCellExperiment::reducedDim(sce, "TRA_Geometric"),
        getdata("runTrex", "runTrex_SCE_TRA_geometric_reduction"),
        tolerance = 1e-2
    )
})

test_that("runTrex adds reduction to Seurat object", {
    data("trex_example")
    set.seed(42)

    trex_example <- runTrex(trex_example,
        chain = "TRA",
        method = "geometric",
        reduction.name = "test_reduction"
    )

    expect_true("test_reduction" %in% names(trex_example@reductions))
    expect_s4_class(trex_example@reductions$test_reduction, "DimReduc")
})

test_that("runTrex accepts custom geometric.method", {
    data("trex_example")
    set.seed(42)

    # Should not error with different substitution matrix
    trex_example <- runTrex(trex_example,
        chain = "TRA",
        method = "geometric",
        geometric.method = "BLOSUM45",
        reduction.name = "test_blosum45"
    )

    expect_true("test_blosum45" %in% names(trex_example@reductions))
})
