test_that("getEpitopeIDsToProbeIDs works", {
  expect_equal(getEpitopeIDsToProbeIDs(c("A_1_2","C_8_9")),
    data.frame(
        Epitope_ID = c("A_1_2", "A_1_2", "C_8_9", "C_8_9"),
        PROBE_ID = c("A;1", "A;2", "C;8", "C;9")
    )
  )
})
