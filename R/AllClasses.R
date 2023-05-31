#' @rdname HERONSequenceDataSet
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.HERONSequenceDataSet <- setClass(
    "HERONSequenceDataSet",
    contains = "SummarizedExperiment")

#' HERONSequenceDataSet object and constructors
#'
#' \code{HERONSequenceDataSet} is a subclass of \code{SummarizedExperiment},
#' used to store the expression values, intermediate calculations, and
#' results of a differential binding code.
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
#' HERONProbeDataSet object and constructors
#'
#' \code{HERONProbeDataSet} is a subclass of \code{RangedSummarizedExperiment}
#' used to TODO
#'
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

#' HERONEpitopeDataSet object and constructors
#'
#' TODO
#'
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

#' HERONProteinDataSet object and constructors
#'
#' TODO
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
HERONProteinDataSet <- function(pvalue, ...) {
    se <- SummarizedExperiment(assays = list(pvalue = pvalue), ...)
    .HERONProteinDataSet(se)
}


