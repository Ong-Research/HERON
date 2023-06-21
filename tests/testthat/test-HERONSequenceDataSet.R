test_that("HERONSequenceDataSet", {
    data("heffron2021_wuhan")
    exprs<-assay(heffron2021_wuhan, "exprs")
    expect_true(validObject(.HERONSequenceDataSet()))
    expect_true(validObject(HERONSequenceDataSet(exprs=exprs)))
    expect_error(HERONSequenceDataSet())
})
