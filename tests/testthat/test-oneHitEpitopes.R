test_that("oneHitEpitopes works", {
    hit_mat = data.frame(
        row.names = c("A_1_1","A_2_2","A_3_3","A_4_4"),
        sample1 = c(TRUE, FALSE, FALSE, TRUE),
        sample2 = c(TRUE, TRUE, FALSE, FALSE),
        sample3 = c(TRUE, TRUE, FALSE, FALSE)
    )
    oh = oneHitEpitopes(hit_mat)
  expect_equal(oh, "A_4_4")
})
