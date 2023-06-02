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
#' results of a differential binding code on the seeuqnce-level.
#'
#' @param exprs binding values with rows as sequences and columns as samples
#' @param ... arguments provided to \code{SummarizedExperiment}, including
#' metadata
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#'
#' exprs <- matrix(1:100,ncol=4)
#' sds <- HERONSequenceDataSet(exprs = exprs)
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
#' used to hold assay information on the probe level
#'
#' @param ... arguments provided to \code{SummarizedExperiment}, including
#' metadata.
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' pds <- HERONProbeDataSet()
#'
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
#' \code{HERONEpitopeDataSet} is a subclass of \code{SummarizedExperiment}
#' used to hold assay information on the epitope-level
#'
#' @param pvalue calculate epitope p-value matrix
#' ... arguments provided to \code{SummarizedExperiment}, including
#' metadata
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
HERONEpitopeDataSet <- function(pvalue, ...) {
    se <- SummarizedExperiment(assays = list(pvalue = pvalue), ...)
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
#' \code{HERONProteinDataSet} is a subclass of \code{SummarizedExperiment}
#' used to hold assay information on the protein-level
#'
#' @param pvalue calculated protein p-value matrix
#' @param ... arguments provided to \code{SummarizedExperiment}, including
#' metadata
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
HERONProteinDataSet <- function(pvalue, ...) {
    se <- SummarizedExperiment(assays = list(pvalue = pvalue), ...)
    .HERONProteinDataSet(se)
}


