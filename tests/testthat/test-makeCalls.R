test_that("makeCalls works", {
    data(heffron2021_wuhan)
    probe_meta <- attr(heffron2021_wuhan, "probe_meta")
    seq_pval_res <- calcCombPValues(heffron2021_wuhan)
    pr_pval_res <- convertSequenceDSToProbeDS(seq_pval_res, probe_meta)
    pr_calls_res <- try(makeProbeCallsPDS(pr_pval_res))

    expect_s4_class(pr_calls_res, "HERONProbeDataSet")

    epi_segments_uniq_res <- findEpitopeSegmentsPDS(
        PDS_obj = pr_calls_res,
        segment_method = "unique"
    )
    epi_padj_uniq <- calcEpitopePValuesProbeDS(
        probe_pds = pr_calls_res,
        epitope_ids = epi_segments_uniq_res,
        metap_method = "wilkinsons_max1"
    )
    epi_calls_res <- try(makeEpitopeCallsEDS(epi_padj_uniq))
    expect_s4_class(epi_calls_res, "HERONEpitopeDataSet")

    prot_padj_uniq <- calcProteinPValuesEpitopeDS(
        epitope_ds = epi_padj_uniq,
        metap_method = "tippets"
    )
    prot_calls_res <- try(makeProteinCalls(prot_padj_uniq))
    expect_s4_class(prot_calls_res, "HERONProteinDataSet")

    expect_no_error(getKofN(pr_calls_res))
    expect_no_error(getKofN(epi_calls_res))
    expect_no_error(getKofN(prot_calls_res))
})
