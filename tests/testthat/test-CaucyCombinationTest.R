test_that("CaucyCombinationTest works", {
    bad_pvals <- c(1,1,1,1.1)
    one_pvals <- c(1,1,1,1)
    p5_pvals <- c(0.5, 0.5, 0.5)
    expect_error(CaucyCombinationTest(bad_pvals))
    expect_warning(CaucyCombinationTest(one_pvals))
    expect_equal(CaucyCombinationTest(c(p5_pvals)), 0.5)
})
