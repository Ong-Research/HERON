test_that("findBlocksProbeT works", {
  expect_equal(rownames(findBlocksProbeT(c("A;1","A;3"))), c("A_1_1","A_3_3"))
})
