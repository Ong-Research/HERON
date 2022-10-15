


#' Normalize a matrix/data.frame using quantile normaliztion
#'
#' @param in_mat matrix or data.frame of numeric values to be normalized
#'
#' @return normalized matrix
#' @export
#'
#' @examples
quantileNormalize<-function(in_mat) {
    #Use limma's version of quantile normalization
    norm_mat = limma::normalizeQuantiles(as.matrix(in_mat));
    norm_mat = as.data.frame(norm_mat);
    colnames(norm_mat) = colnames(in_mat);
    rownames(norm_mat) = rownames(in_mat);
    return(norm_mat);
}
