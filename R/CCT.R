#' An analytical p-value combination method using the Cauchy distribution
#'
#' The \code{CaucyCombinationTest} function takes in a numeric vector of
#' p-values, a numeric vector of non-negative weights, and return the aggregated
#' p-value using Cauchy method.
#' @param pvals a numeric vector of p-values, where each of the element is
#' between 0 to 1, to be combined.
#' @return the aggregated p-value combining p-values from the
#' vector \code{pvals}.
#' @examples pvalues <- c(2e-02,4e-04,0.2,0.1,0.8)
#' @examples CaucyCombinationTest(pvals=pvalues)
#' @references Liu, Y., & Xie, J. (2020). Cauchy combination test:
#' a powerful test
#' with analytic p-value calculation under arbitrary dependency structures.
#' \emph{Journal of the American Statistical Association 115}(529), 393-402.
#' (\href{https://doi.org/10.1080/01621459.2018.1554485}{pub})
#' https://github.com/xihaoli/STAAR
#' @noRd
CaucyCombinationTest <- function(pvals){
    #Use equal weights for now.
    weights <- rep(1/length(pvals),length(pvals))
    #### check if there are very small non-zero p-values
    is.small <- (pvals < 1e-16)
    if (sum(is.small) == 0){
        cct.stat <- sum(weights*tan((0.5-pvals)*pi))
    }else{
        cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
        cct.stat <- cct.stat +
            sum(weights[!is.small] *
            tan((0.5-pvals[!is.small])*pi))
    }
    #### check if the test statistic is very large.
    if(cct.stat > 1e+15){
        pval <- (1/cct.stat)/pi
    }else{
        pval <- 1-stats::pcauchy(cct.stat)
    }
    return(pval)
}
