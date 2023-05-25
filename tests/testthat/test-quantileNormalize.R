test_that("quantileNormalize works", {

    expect_error(
        {
            data("heffron2021_wuhan")
            quantileNormalize(assay(heffron2021_wuhan,"exprs"))
        }
    )
})
