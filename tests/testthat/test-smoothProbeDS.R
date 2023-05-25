test_that("smoothProbeDS works", {
    expect_error({
        data("heffron2021_wuhan")
        smoothProbeDS(heffron2021_wuhan)
    })
    expect_no_error({
        data("heffron2021_wuhan")
        pm <- attr(heffron2021_wuhan, "probe_meta")
        pds <- convertSequenceDSToProbeDS(heffron2021_wuhan, pm)
        smoothProbeDS(pds)
    })
    expect_error({
        data("heffron2021_wuhan")
        pm <- attr(heffron2021_wuhan, "probe_meta")
        pds <- convertSequenceDSToProbeDS(heffron2021_wuhan, pm)
        smoothProbeDS(pds, w=3)
    })

})
