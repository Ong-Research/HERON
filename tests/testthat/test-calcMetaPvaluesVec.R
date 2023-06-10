test_that("calcMetaPvaluesVec works", {
    bad_pvals <- c(1,1,1,1.1) #Test case when there are bad pvalues
    one_pvals <- c(1,1,1,1) #Test case when there are all 1 pvalues
    p5_pvals <- c(0.5, 0.5, 0.5) #Test case when there are all 0.5 pvalues
    #Test case when there is one small pvalue
    onesmall_pvals <- c(1e-17, 0.5, 0.5)
    #Test case when there is a single pvalue
    single_pvalue <- c(0.1)
    #Test case when there is no valid pvalues
    na_pvalue <- c(NA)


    type_vec <- c(
        rep("Bad", length(bad_pvals)),
        rep("Na", length(na_pvalue)),
        rep("One", length(one_pvals)),
        rep("OneSmall", length(onesmall_pvals)),
        rep("P5", length(p5_pvals)),
        rep("Single", length(single_pvalue))
    )

    pvals <- c(
        bad_pvals, na_pvalue, one_pvals,
        onesmall_pvals, p5_pvals, single_pvalue
    )

    pvals_df <- data.frame(
        pvalue = pvals,
        Type = type_vec
    )

    by_list <- list(Type = type_vec)

    nelements <- c(
        4, # Bad
        1, # Na
        4, # One
        3, # OneSmall
        3, # P5
        1  # Single
    )

    #Test taking the max
    mpvals <- calcMetaPValuesVec(pvals, by_list = by_list, "max")
    expect_equal(mpvals$NElements, nelements)
    expect_equal(
        mpvals$Meta.pvalue,
        c(
            1.0, #Bad
            1.0, #Na
            1.0, #One
            0.5, #OneSmall
            0.5, #P5
            0.1 #Single
        ),
        tolerance = 1e-15
    )

    #Test taking the min
    mpvals <- calcMetaPValuesVec(pvals, by_list = by_list, "min")
    expect_equal(mpvals$NElements, nelements)
    expect_equal(
        mpvals$Meta.pvalue,
        c(
            1.0, #Bad
            1.0, #Na
            1.0, #One
            1e-17, #OneSmall
            0.5, #P5
            0.1 #Single
        ),
        tolerance = 1e-15
    )

    #Test wilkinson's min on one the smallest element
    mpvals <- calcMetaPValuesVec(pvals, by_list = by_list, "wmin1")
    expect_equal(mpvals$NElements, nelements)
    expect_equal(
        mpvals$Meta.pvalue,
        c(
            1.0, #Bad
            1.0, #Na
            1.0, #One
            2.9e-17, #OneSmall
            0.875, #P5
            0.1 #Single
        ),
        tolerance = 1e-15
    )

    #Test Wilkinson's min on the 2nd smallest element
    mpvals <- calcMetaPValuesVec(pvals, by_list = by_list, "wmin2")
    expect_equal(mpvals$NElements, nelements)
    expect_equal(
        mpvals$Meta.pvalue,
        c(
            1.0, #Bad
            1.0, #Na
            1.0, #All ones, so 1 should be returned
            0.5, #OneSmall, 2nd minimum, so 0.5
            0.5, #P5
            0.1 #Single
        ),
        tolerance = 1e-15
    )

    #Test CaucyCombinationTest
    expect_no_error(calcMetaPValuesVec(pvals, by_list = by_list, "cct"))
    mpvals <- calcMetaPValuesVec(pvals, by_list = by_list, "cct")
    expect_equal(mpvals$NElements, nelements)
    expect_equal(
        mpvals$Meta.pvalue,
        c(
            1.0, #Bad
            1.0, #Na
            1.0, #All ones, so 1 should be returned
            2.9e-17, #OneSmall,
            0.5, #P5
            0.1 #Single
        ),
        tolerance = 1e-15
    )


})
