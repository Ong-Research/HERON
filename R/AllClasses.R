#' @rdname HERONSequenceDataSet
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.HERONSequenceDataSet <- setClass(
    "HERONSequenceDataSet",
    contains = "SummarizedExperiment")

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
HERONSequenceDataSet <- function(exprs, ...) {
    se <- SummarizedExperiment(assays = list(exprs = exprs), ...)
    .HERONSequenceDataSet(se)
}


#' @rdname HERONProbeDataSet
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.HERONProbeDataSet <- setClass(
    "HERONProbeDataSet",
    contains = "RangedSummarizedExperiment"
)

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
HERONProbeDataSet <- function(...) {
    rse <- SummarizedExperiment(...)
    .HERONProbeDataSet(rse)
}

#' @rdname HERONEpitopeDataSet
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.HERONEpitopeDataSet <- setClass(
    "HERONEpitopeDataSet",
    contains = "SummarizedExperiment"
)

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
HERONEpitopeDataSet <- function(pvalues, ...) {
    se <- SummarizedExperiment(assays = list(pvalue = pvalues), ...)
    .HERONEpitopeDataSet(se)
}


#' @rdname HERONProteinDataSet
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.HERONProteinDataSet <- setClass(
    "HERONProteinDataSet",
    contains = "SummarizedExperiment"
)

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
HERONProteinDataSet <- function(pvalues, ...) {
    se <- SummarizedExperiment(assays = list(pvalue = pvalues), ...)
    .HERONProteinDataSet(se)
}


