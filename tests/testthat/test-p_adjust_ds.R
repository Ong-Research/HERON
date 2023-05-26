test_that("p_adjust_ds works", {
    data("heffron2021_wuhan")
    exprs <- assay(heffron2021_wuhan, "exprs")
    expect_error(p_adjust_ds(exprs))
    expect_error(p_adjust_ds(heffron2021_wuhan))
})
