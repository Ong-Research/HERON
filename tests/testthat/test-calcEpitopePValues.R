test_that("calcEpitopePValues works", {
    data("heffron2021_wuhan")
    expect_error(calcEpitopePValues(heffron2021_wuhan))

    pval_seq_res <- calcCombPValues(heffron2021_wuhan)
    pval_pr_res <- convertSequenceDSToProbeDS(pval_seq_res)
    calls_res <- makeProbeCalls(pval_pr_res)
    segments_res <- findEpitopeSegments(calls_res, "unique")
    expect_no_error(calcEpitopePValues(calls_res, segments_res))
    mmethod_vec <- c(
        "min_bonf"
        , "min"
        , "max"
        , "fisher"
        , "hmp"
        , "tippetts"
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
            calcEpitopePValues(calls_res, segments_res, mmethod)
        )
    }

    #Test if bad meta p-value parameter passed in
    expect_error(calcEpitopePValues(calls_res, segments_res, "bad"))


    #Test with NA pvalue
    pvals <- assay(calls_res, "pvalue")

    #select rows from the probe pvalue assay that are involved in an epitope
    segments_probe <- getEpitopeIDsToProbeIDs(segments_res)


    #Do these tests a # of times since we are randomly picking cells
    for (rep_idx in seq_len(3)) {
        ridx <- sample(segments_probe$PROBE_ID, 1)
        cidx <- sample(colnames(pvals), 1)
        pvals_new <- pvals

        #Test NA p-value
        pvals_new[ridx, cidx] <- NA
        assay(calls_res, "pvalue") <- pvals_new
        expect_warning(calcEpitopePValues(calls_res, segments_res, "max"))

        #Test >1 p-value
        pvals_new <- pvals
        pvals_new[ridx, cidx] <- 1.1
        assay(calls_res, "pvalue") <- pvals_new
        expect_warning(calcEpitopePValues(calls_res, segments_res, "max"))

        #Test < 0 p-value
        pvals_new <- pvals
        pvals_new[ridx, cidx] <- -0.1
        assay(calls_res, "pvalue") <- pvals_new
        expect_warning(calcEpitopePValues(calls_res, segments_res, "max"))
    }

})
