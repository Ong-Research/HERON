test_that("calcProteinPValues works", {
    data(heffron2021_wuhan)
    pval_seq_res <- calcCombPValues(heffron2021_wuhan)
    pval_pr_res <- convertSequenceDSToProbeDS(pval_seq_res)
    calls_res <- makeProbeCallsPDS(pval_pr_res)
    segments_res <- findEpitopeSegmentsPDS(calls_res, "unique")
    epval_res <- calcEpitopePValuesProbeDS(calls_res, segments_res)
    expect_no_error(calcProteinPValuesEpitopeDS(epval_res))
})
