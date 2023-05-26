test_that("calcEpitopePValuesProbeDS works", {
    data("heffron2021_wuhan")
    expect_error(calcEpitopePValuesProbeDS(heffron2021_wuhan))


})
