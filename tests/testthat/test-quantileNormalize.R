test_that("quantileNormalize works", {
    data("heffron2021_wuhan")
    expect_error(quantileNormalize(assay(heffron2021_wuhan,"exprs")))
    expect_no_error(quantileNormalize(heffron2021_wuhan))
})
