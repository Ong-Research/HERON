test_that("calcMetaPvaluesVec works", {
    bad_pvals <- c(1,1,1,1.1)
    one_pvals <- c(1,1,1,1)
    p5_pvals <- c(0.5, 0.5, 0.5)

    type_vec <- c(
        rep("Bad", length(bad_pvals)),
        rep("One", length(one_pvals)),
        rep("p5", length(p5_pvals))
    )

    pvals <- c(bad_pvals, one_pvals, p5_pvals)

    pvals_df <- data.frame(
        pvalue = pvals,
        Type = type_vec
    )

    by_list <- list(Type = type_vec)

    mpvals <- calcMetaPValuesVec(pvals, by_list = by_list, "max")

    expect_equal(mpvals$NElements,c(4,4,3))
    expect_equal(mpvals$Meta.pvalue, c(1.0, 1.0, 0.5))

    mpvals <- calcMetaPValuesVec(pvals, by_list = by_list, "min")
    expect_equal(mpvals$NElements,c(4,4,3))
    expect_equal(mpvals$Meta.pvalue, c(1.0, 1.0, 0.5))

    mpvals <- calcMetaPValuesVec(pvals, by_list = by_list, "wmin1")
    expect_equal(mpvals$NElements,c(4,4,3))
    expect_equal(mpvals$Meta.pvalue, c(1.0, 1.0, 0.875))



})
