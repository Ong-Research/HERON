
#' Get the amino-acid starting position of the probe within the protein.
#'
#' @param probes vector of probes (i.e. c("A;1", "A;2"))
#'
#' @return vector of integers indicating the respective starting locations
#'  of the probes with their associated proteins
#' @export
#'
#' @examples
getProteinStart<-function(probes) {
    probe.list = strsplit(probes, ";");
    protein.start.list = lapply(
        probe.list,
        function(x) {

            return(x[length(x)]);
        }
    );
    return(as.integer(unlist(protein.start.list)));
}


#' Get Protein Label from Probe
#'
#' @param probes vector of probes (i.e. c("A;1", "A;2"))
#'
#' @return vector of strings indicating the protein associated with the
#' respective probes
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
#' @return For each protein, the estimating tiling (spacing) of the probes
#' across the amino acid sequence.
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



#' Cap vector at minimum/maximum values
#'
#' @param val vector of values to cap
#' @param min.value minimum value
#' @param max.value maximum value
#'
#' @return vector of capped values
#' @export
#'
#' @examples
min_max<-function(val, min.value, max.value) {
    val[val < min.value] = min.value;
    val[val > max.value] = max.value;

    return(val);
}


#' Concatenate sequences together based upon their start positions.
#'
#' @param positions
#' @param sequences
#' @param debug
#'
#' @return concatenated sequence (character)
#' @export
#'
#' @examples
catSequences <- function (positions, sequences, debug = FALSE)
{
    if (length(sequences) == 1) {
        return(sequences)
    }
    positions = positions - positions[1] + 1
    seq = rep("", 10000) #TODO - make this more tolerant.
    for (idx in 1:length(positions)) {
        pos = positions[idx]
        aa_vec = strsplit(sequences[idx], "")[[1]]
        startp = pos
        endp = pos + length(aa_vec) - 1
        seq[startp:endp] = aa_vec
    }
    seq = seq[seq != ""]
    seqc = paste0(seq, collapse = "", sep = "")
    return(seqc)
}


