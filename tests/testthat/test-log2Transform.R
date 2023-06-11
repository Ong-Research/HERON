test_that("log2Transform works", {
    data("heffron2021_wuhan")

    #Test bad arguments
    expect_error(log2Transform(assay(heffron2021_wuhan,"exprs")))
    expect_error(log2Transform(.HERONEpitopeDataSet()))

    #Test result
    test_se <- heffron2021_wuhan
    assay(test_se, "exprs") <- 2^assay(test_se, "exprs")

    expect_no_error(test_se <- log2Transform(test_se))
    expect_equal(assay(test_se, "exprs"), assay(heffron2021_wuhan, "exprs"))
})
