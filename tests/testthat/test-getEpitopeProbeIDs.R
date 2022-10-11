test_that("getEpitopeProbeIDs works", {
  expect_equal(getEpitopeProbeIDs("A_1_5"), c("A;1","A;2","A;3","A;4","A;5"))
})
