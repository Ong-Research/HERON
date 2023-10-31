test_that("makeCalls works", {
    data(heffron2021_wuhan)
    seq_pval_res <- calcCombPValues(heffron2021_wuhan)
    pr_pval_res <- convertSequenceDSToProbeDS(seq_pval_res)
    pr_calls_res <- try(makeProbeCalls(pr_pval_res))

    expect_s4_class(pr_calls_res, "HERONProbeDataSet")

    epi_segments_uniq_res <- findEpitopeSegments(
        PDS_obj = pr_calls_res,
        segment_method = "unique"
    )
    epi_padj_uniq <- calcEpitopePValues(
        probe_pds = pr_calls_res,
        epitope_ids = epi_segments_uniq_res,
        metap_method = "wilkinsons_max1"
    )
    epi_calls_res <- try(makeEpitopeCalls(epi_padj_uniq))
    expect_s4_class(epi_calls_res, "HERONEpitopeDataSet")

    prot_padj_uniq <- calcProteinPValues(
        epitope_ds = epi_padj_uniq,
        metap_method = "tippetts"
    )
    prot_calls_res <- try(makeProteinCalls(prot_padj_uniq))
    expect_s4_class(prot_calls_res, "HERONProteinDataSet")

    expect_no_error(getKofN(pr_calls_res))
    expect_no_error(getKofN(epi_calls_res))
    expect_no_error(getKofN(prot_calls_res))

    expect_no_error(makeProteinCalls(prot_padj_uniq, one_hit_filter = TRUE))


})
