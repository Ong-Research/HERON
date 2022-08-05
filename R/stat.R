#' Calculate minimum FDRs across samples
#'
#' This code calculates the minimum FDR threshold needed to be a hit in X out of N samples
#' Where samples are in the columns
#' When additional stats is selected, also calculates the mean and median
#'
#' @param fdrs matrix of adjusted p-values
#' @param additional_stats indicator of additional stats
#' @param sort sort the resultant matrix by the last N column (increasing)
#'
#' @return
#' @export
#'
#' @examples
calcMinFDR<-function(fdrs, additional_stats=TRUE, sort=TRUE) {
    fdrs2 = as.data.frame(fdrs,stringsAsFactors=FALSE);
    ncols = ncol(fdrs);
    #cat("Calculating minFDRs\n")

    fdrs2 = t(apply(fdrs2, 1, function(l) { return(l[order(l,decreasing=FALSE)])} ))
    fdrs2 = as.data.frame(fdrs2, stringsAsFactors=FALSE)
    colnames(fdrs2) = paste0("n",1:ncol(fdrs2))
    if (additional_stats) {
        fdrs2$meanFDR = rowMeans(fdrs);
        fdrs2$medFDR = matrixStats::rowMedians(fdrs);
    }
    n_col = paste0("n",ncol(fdrs))
    if (sort) {
      fdrs2 = fdrs2[order(fdrs2[,n_col],decreasing=FALSE),]
    }
    return(fdrs2);

}
