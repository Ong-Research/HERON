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
#'
#' metadata can contain a probe DataFrame, that maps sequences
#' (column PROBE_SEQUENCE) to probe identifiers ( column PROBE_ID)
#'
#' @return HERONSequenceDataSet object
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' exprs <- matrix(seq_len(100),ncol=4)
#' colnames(exprs) <- c("C1", "C2", "C3", "C4")
#' sds <- HERONSequenceDataSet(exprs = exprs)
HERONSequenceDataSet <- function(exprs, ...) {
    se <- SummarizedExperiment(assays = list(exprs = exprs), ...)
    if (ncol(colData(se)) == 0) {
        colData(se) <- DataFrame(
            row.names = colnames(exprs),
            SampleName = colnames(exprs),
            ptid = colnames(exprs),
            visit = rep("post", ncol(exprs))
        )
    }
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
#' @return HERONProbeDataSet object
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' pds <- HERONProbeDataSet()
#'
HERONProbeDataSet <- function(...) {
    rse <- SummarizedExperiment(...)
    if (ncol(colData(rse)) == 0) {
        colData(rse) <- DataFrame(
            row.names = colnames(exprs),
            SampleName = colnames(exprs),
            ptid = colnames(exprs),
            visit = rep("post", ncol(exprs))
        )
    }
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
#' @param ... arguments provided to \code{SummarizedExperiment}, including
#' metadata
#'
#' @return HERONEpitopeDataSet object
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' pval <- matrix(runif(100),ncol=4)
#' HERONEpitopeDataSet(pvalue = pval)
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
#'
#' @return HERONProteinDataSet object
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' pval <- matrix(runif(100), ncol=4)
#' HERONProteinDataSet(pvalue = pval)
HERONProteinDataSet <- function(pvalue, ...) {
    se <- SummarizedExperiment(assays = list(pvalue = pvalue), ...)
    .HERONProteinDataSet(se)
}


