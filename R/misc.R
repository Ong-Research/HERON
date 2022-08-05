
#' Get the amino-acid starting position of the probe within the protein.
#'
#' @param probes vector of probes (i.e. c("A;1", "A;2"))
#'
#' @return
#' @export
#'
#' @examples
getProteinStart<-function(probe.ids) {
    probe.list = strsplit(probe.ids, ";");
    protein.start.list = lapply(
        probe.list,
        function(x) {

            return(x[length(x)]);
        }
    );
    return(as.integer(unlist(protein.start.list)));
}


#' Title
#'
#' @param probes vector of probes (i.e. c("A;1", "A;2"))
#'
#' @return
#' @export
#'
#' @examples
getProteinLabel<-function(probes) {

    probe.list = strsplit(probes,";");
    protein.list = lapply(
        probe.list, function(x){
            n = length(x)-1;
            return(paste(x[1:n],sep=";",collapse=";"))
        }
    )
    protein.labels = unlist(protein.list);
    return(protein.labels);

}


#' Get Protein Tiling
#'
#' Given a set of probes, estimate the tiling of the probes
#' across the protein.
#'
#' @param probes vector of probes (i.e. A;1, A;2)
#' @param return.vector Return result as vector or return as data.frame
#'
#' @return
#' @export
#'
#' @examples
getProteinTiling<-function(probes, return.vector=TRUE) {
    Pos = getProteinStart(probes);
    Protein = getProteinLabel(probes);
    tiling.df = stats::aggregate(
        as.integer(Pos),
        by = list(Protein = Protein),
        function(l) {
            if (length(l) > 1) {
                warning("Multiple tiling options detected")
                l = l[order(l, decreasing = FALSE)]
                return(l[2] - l[1])
            }
            return(-1)
        }
    )
    ans = tiling.df
    colnames(ans)[2] = "Tiling";
    if (return.vector) {
        ans.vec = ans$Tiling;
        names(ans.vec) = ans$Protein;
        ans = ans.vec;
    }
    return(ans);
}



