# test script for Trex_matrix.R

test_that("Trex_matrix works with encoder method", {
    data("trex_example")

    set.seed(42)

    # Test TRB chain with VAE and atchleyFactors
    result <- Trex_matrix(trex_example,
        chain = "TRB",
        method = "encoder",
        encoder.model = "VAE",
        encoder.input = "atchleyFactors",
        verbose = FALSE
    )

    expect_true(is.data.frame(result))
    expect_true(nrow(result) > 0)
    expect_true(all(grepl("^Trex_", colnames(result))))

    expect_equal(
        result,
        getdata("runTrex", "maTrex_TRB_VAE_AF"),
        tolerance = 1e-2
    )
})

test_that("Trex_matrix works with different encoder inputs", {
    data("trex_example")

    set.seed(42)

    # Test TRA chain with AE and OHE
    result <- Trex_matrix(trex_example,
        chain = "TRA",
        method = "encoder",
        encoder.model = "AE",
        encoder.input = "OHE",
        verbose = FALSE
    )

    expect_equal(
        result,
        getdata("runTrex", "maTrex_TRA_AE_OHE"),
        tolerance = 1e-2
    )
})

test_that("Trex_matrix works with geometric method", {
    data("trex_example")

    set.seed(42)

    # Test geometric with default BLOSUM62
    result <- Trex_matrix(trex_example,
        chain = "TRA",
        method = "geometric",
        verbose = FALSE
    )

    expect_true(is.data.frame(result))
    expect_true(nrow(result) > 0)
    expect_equal(ncol(result), 20)
    expect_true(all(grepl("^Trex_", colnames(result))))

    expect_equal(
        result,
        getdata("runTrex", "maTrex_Geometric"),
        tolerance = 1e-2
    )
})

test_that("Trex_matrix allows custom geometric.method", {
    data("trex_example")

    set.seed(42)

    # Test with different substitution matrix
    result_blosum62 <- Trex_matrix(trex_example,
        chain = "TRA",
        method = "geometric",
        geometric.method = "BLOSUM62",
        verbose = FALSE
    )

    result_blosum45 <- Trex_matrix(trex_example,
        chain = "TRA",
        method = "geometric",
        geometric.method = "BLOSUM45",
        verbose = FALSE
    )

    # Results should be different with different matrices
    expect_false(identical(result_blosum62, result_blosum45))
})

test_that("Trex_matrix allows custom theta", {
    data("trex_example")

    set.seed(42)

    result_pi <- Trex_matrix(trex_example,
        chain = "TRA",
        method = "geometric",
        geometric.theta = pi,
        verbose = FALSE
    )

    result_pi_half <- Trex_matrix(trex_example,
        chain = "TRA",
        method = "geometric",
        geometric.theta = pi / 2,
        verbose = FALSE
    )

    # Results should differ with different theta values
    expect_false(identical(result_pi, result_pi_half))
})

test_that("Trex_matrix validates chain argument", {
    data("trex_example")

    expect_error(
        Trex_matrix(trex_example, chain = "invalid"),
        "'arg' should be one of"
    )
})

test_that("Trex_matrix validates method argument", {
    data("trex_example")

    expect_error(
        Trex_matrix(trex_example, method = "invalid"),
        "'arg' should be one of"
    )
})

test_that("Trex_matrix returns correct dimensions", {
    data("trex_example")

    set.seed(42)

    # Encoder method should return 30 dimensions
    result_encoder <- Trex_matrix(trex_example,
        chain = "TRB",
        method = "encoder",
        verbose = FALSE
    )
    expect_equal(ncol(result_encoder), 30)

    # Geometric method should return 20 dimensions
    result_geometric <- Trex_matrix(trex_example,
        chain = "TRB",
        method = "geometric",
        verbose = FALSE
    )
    expect_equal(ncol(result_geometric), 20)
})
