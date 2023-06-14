test_that("calcCombPvalues", {
    data("heffron2021_wuhan")
    exprs <- assay(heffron2021_wuhan, "exprs")
    expect_error(calcCombPValues(exprs))

    ## Test completion of default parameters
    expect_no_error(calcCombPValues(heffron2021_wuhan))
    expect_no_error(calcCombPValues(heffron2021_wuhan, use="t"))
    expect_no_error(calcCombPValues(heffron2021_wuhan, use="z"))

    expect_no_error(calcCombPValues(heffron2021_wuhan, t_sd_shift = 1))
    expect_no_error(calcCombPValues(heffron2021_wuhan, z_sd_shift = 1))


    ## Test bad parameter
    expect_error(calcCombPValues(heffron2021_wuhan, use="k"))

    ## Test using probe dataset
    pr_ds <- convertSequenceDSToProbeDS(heffron2021_wuhan)
    expect_no_error(calcCombPValues(pr_ds))

    ## Test using reduced colData
    colData_red <- colData(pr_ds)
    post_idx <- which(colData_red$visit == "post")
    colData_red <- colData_red[-post_idx[1],]
    expect_no_error(calcCombPValues(pr_ds, colData_in = colData_red))

    ## Test paired t-test method
    colData_paired <- colData(heffron2021_wuhan)

    ### Make some samples paired by setting the sample ids
    pre_idx <- which(colData_paired$visit == "pre")
    colData_post <- colData_paired[colData_paired$visit == "post",]
    new_ids <- colData_post$Sample_ID[seq_len(5)]
    colData_paired$ptid[pre_idx[seq_len(5)]] = new_ids

    paired_ds <- heffron2021_wuhan
    colData(paired_ds) <- colData_paired

    ### calculate p-values
    expect_no_error(calcCombPValues(obj = paired_ds, t_paired = TRUE))






})
