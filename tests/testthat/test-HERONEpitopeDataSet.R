test_that("HERONEpitopeDataSet works", {

    pvalues <- matrix(1, nrow=2, ncol=2)

    expect_true(validObject(.HERONEpitopeDataSet()))
    expect_true(validObject(HERONEpitopeDataSet(pvalues = pvalues)))
})
