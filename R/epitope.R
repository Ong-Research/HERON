


#' Obtain Protein Id from Epitope ID
#'
#' Format of EpitopeID is A_B_C, where A is the protein label
#' B is the protein start position of the first probe in the epitope and
#' C is the protein start position of the last prboe in the epitope.
#'
#' @param epitope_ids vector of epitope_ids
#'
#' @return vector of protein labels
#' @export
#'
#' @examples
getEpitopeProtein<-function(epitope_ids) {
    ans = unlist(
        lapply(
            strsplit(epitope_ids,"_"),
            function(l) {
                ans = l[1:(length(l)-2)]
                return(paste(ans,collapse="_"))
            }
        )
    )
    return(ans);
}

#' Obtain first probe's protein start position from Epitope ID
#'
#'
#' @param epitope_ids vector of epitope ids
#'
#' @return vector of integers indicating first probe start positions in the
#' epitope(s)
#'
#' @export
#'
#' @examples
getEpitopeStart<-function(epitope_ids) {
    ans = unlist(
        lapply(
            strsplit(epitope_ids,"_"),
            function(l) {
                return(l[length(l)-1])
            }
        )
    )
    return(as.integer(ans));
}

#' Obtain last probe's protein start position from EpitopeID
#'
#' @param epitope_ids vector of epitope ids
#'
#' @return vector of integers indicating the last probe protein start position
#' @export
#'
#' @examples
getEpitopeStop<-function(epitope_ids) {
    ans = unlist(
        lapply(
            strsplit(epitope_ids,"_"),
            function(l) {
                return(l[length(l)])
            }
        )
    )
    return(as.integer(ans));
}

#' Get the vector of probes from an epitope id
#'
#' @param epitope_id EpitopeID to obtain probes from
#' @param tiling Tiling of the probes across the protein (default 1)
#'
#' @return vector of probe_ids that are contained within the epitope
#' @export
#'
#' @examples
getEpitopeProbeIDs<-function(epitope_id, tiling=1) {
    protein = getEpitopeProtein(epitope_id);
    start = getEpitopeStart(epitope_id)
    stop = getEpitopeStop(epitope_id);

    ans = paste0(protein,";",seq(from=start, to=stop, by=tiling));
    return(ans);
}

#' Obtain the epitope ids from the vectors of protein, first and last probes
#'
#' @param protein vector of proteins
#' @param start vector of first probe protein start positions
#' @param stop vector of last probe protein start positions
#'
#' @return vector of epitope ids
#' @export
#'
#' @examples
getEpitopeID<-function(protein, start, stop) {
    return(paste(protein,start,stop, sep="_"))
}
