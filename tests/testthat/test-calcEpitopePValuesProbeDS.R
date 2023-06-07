test_that("calcEpitopePValuesProbeDS works", {
    data("heffron2021_wuhan")
    expect_error(calcEpitopePValuesProbeDS(heffron2021_wuhan))

    pval_seq_res <- calcCombPValues(heffron2021_wuhan)
    pval_pr_res <- convertSequenceDSToProbeDS(pval_seq_res)
    calls_res <- makeProbeCallsPDS(pval_pr_res)
    segments_res <- findEpitopeSegmentsPDS(calls_res, "unique")
    expect_no_error(calcEpitopePValuesProbeDS(calls_res, segments_res))
    mmethod_vec <- c(
        "min_bonf"
        , "min"
        , "max"
        , "fisher"
        , "hmp"
        , "tippets"
        , "wmin2"
        , "wmin3"
        , "wmin4"
        , "wmin5"
        , "wmax1"
        , "wmax2"
        , "cct"
    )
    for (mmethod in mmethod_vec) {
        expect_no_error(
            calcEpitopePValuesProbeDS(calls_res, segments_res, mmethod)
        )
    }

    expect_error(calcEpitopePValuesProbeDS(calls_res, segments_res, "bad"))
})
