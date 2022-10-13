test_that("multiplication works", {
    hit_mat = data.frame(
        row.names = c("A;1","A;2","A;3","A;4"),
        sample1 = c(TRUE, FALSE, FALSE, TRUE),
        sample2 = c(TRUE, TRUE, FALSE, FALSE),
        sample3 = c(TRUE, TRUE, FALSE, FALSE)
    )
  expect_equal(oneHitProbes(hit_mat),  c("A;4"))
})
