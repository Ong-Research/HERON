test_that("getEpitopeID works", {
  expect_equal(getEpitopeID("A", 1, 2), "A_1_2")
})
