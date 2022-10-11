test_that("oneProbeEpitopes works", {
  expect_equal(oneProbeEpitopes(c("A_1_1", "B_1_1","C_1_2")),
               c(TRUE, TRUE, FALSE))
})
