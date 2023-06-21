test_that("calcProteinPValues works", {
    data(heffron2021_wuhan)
    pval_seq_res <- calcCombPValues(heffron2021_wuhan)
    pval_pr_res <- convertSequenceDSToProbeDS(pval_seq_res)
    calls_res <- makeProbeCalls(pval_pr_res)
    segments_res <- findEpitopeSegments(calls_res, "unique")
    epval_res <- calcEpitopePValues(calls_res, segments_res)
    expect_no_error(calcProteinPValues(epval_res))
})
