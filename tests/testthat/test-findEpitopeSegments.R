test_that("findEpitopeSegments works", {
    data(heffron2021_wuhan)
    probe_meta <- attr(heffron2021_wuhan, "probe_meta")
    seq_ds_qn <- quantileNormalize(heffron2021_wuhan)
    seq_pval_res <- calcCombPValues(seq_ds_qn)
    probe_pval_res <- convertSequenceDSToProbeDS(seq_pval_res, probe_meta)
    probe_calls_res <- makeProbeCallsPDS(probe_pval_res)
    expect_no_error(
        findEpitopeSegmentsPDS(probe_calls_res,"unique")
    )
    expect_no_error(
        findEpitopeSegmentsPDS(
            PDS_obj = probe_calls_res,
            segment_method = "hclust",
            segment_score_type = "binary",
            segment_dist_method = "hamming",
            segment_cutoff = "silhouette"
        )
    )
    expect_no_error(findEpitopeSegmentsPDS(probe_calls_res,"skater","binary","hamming", "silhouette"))
    expect_no_error(findEpitopeSegmentsPDS(probe_calls_res,"hclust","zscore","euclidean", "silhouette"))
    expect_no_error(findEpitopeSegmentsPDS(probe_calls_res,"skater","zscore","euclidean", "silhouette"))

})
