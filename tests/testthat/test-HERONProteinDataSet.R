test_that("HERONProteinDataSet works", {
    pvalues <- matrix(1, nrow=2, ncol=2)

    expect_true(validObject(.HERONProteinDataSet()))
    expect_true(validObject(HERONProteinDataSet(pvalues = pvalues)))
})
