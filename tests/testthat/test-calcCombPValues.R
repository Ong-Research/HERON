test_that("calcCombPvalues", {
    data("heffron2021_wuhan")
    exprs <- assay(heffron2021_wuhan, "exprs")
    expect_error(calcCombPValues(exprs))
})
