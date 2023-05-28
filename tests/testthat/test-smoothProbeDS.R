test_that("smoothProbeDS works", {
    data("heffron2021_wuhan")
    pm <- attr(heffron2021_wuhan, "probe_meta")
    pds <- convertSequenceDSToProbeDS(heffron2021_wuhan, pm)

    expect_error(smoothProbeDS(heffron2021_wuhan))
    expect_no_error(
        {
            smoothProbeDS(pds)
        }
    )
    expect_error(
        {
            smoothProbeDS(pds, w=3)
        }
    )

    expect_no_error(
        {
            smoothProbeDS(pds, w=0)
        }
    )

    #Test smoothing with the prescence of missing values
    sds_na <- heffron2021_wuhan
    exprs <- assay(sds_na,"exprs")

    x <- sample(seq_len(nrow(exprs)), size = 1)
    y <- sample(seq_len(ncol(exprs)), size = 1)

    exprs[x, y] <- NA

    assay(sds_na, "exprs") <- exprs
    pds_na <- convertSequenceDSToProbeDS(sds_na, pm)

    expect_no_error(smoothProbeDS(pds_na))

})
