test_that("findEpitopeSegments works", {
    data(heffron2021_wuhan)
    seq_ds_qn <- quantileNormalize(heffron2021_wuhan)
    seq_pval_res <- calcCombPValues(seq_ds_qn)
    probe_pval_res <- convertSequenceDSToProbeDS(seq_pval_res)
    probe_calls_res <- makeProbeCalls(probe_pval_res)
    expect_no_error(
        findEpitopeSegments(probe_calls_res,"unique")
    )
    expect_no_error(
        findEpitopeSegments(
            PDS_obj = probe_calls_res,
            segment_method = "hclust",
            segment_score_type = "binary",
            segment_dist_method = "hamming",
            segment_cutoff = "silhouette"
        )
    )
    expect_no_error(
        findEpitopeSegments(
            PDS_obj = probe_calls_res,
            segment_method = "skater",
            segment_score_type = "binary",
            segment_dist_method = "hamming",
            segment_cutoff = "silhouette"
        )
    )
    expect_no_error(
        findEpitopeSegments(
            PDS_obj = probe_calls_res,
            segment_method = "hclust",
            segment_score_type = "zscore",
            segment_dist_method = "euclidean",
            segment_cutoff = "silhouette"
        )
    )
    expect_no_error(
        findEpitopeSegments(
            PDS_obj = probe_calls_res,
            segment_method = "skater",
            segment_score_type = "zscore",
            segment_dist_method = "euclidean",
            segment_cutoff = "silhouette"
        )
    )

})
