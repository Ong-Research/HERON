
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
            return(paste(x[seq_len(n)],sep=";",collapse=";"))
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
                l = l[order(l, decreasing = FALSE)]
                return(l[2] - l[1]) #Assume that the tiling is consistent across
                #the probes.
            }
            #Not enough probes for the protein to calculate a tiling.
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
#' @param positions start positions of probes in protein
#' @param sequences probe sequences of probes
#'
#' @return concatenated sequence (character)
#' @export
#'
#' @examples
catSequences <- function (positions, sequences) {
    if (length(sequences) == 1) {
        return(sequences)
    }
    positions = positions - positions[1] + 1
    seq = rep("", 10000) #TODO - make this more tolerant.
    for (idx in seq_len(length(positions))) {
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


#' Convert p-value matrix to a z-score matrix
#'
#' @param mat.in matrix of p-values
#' @param one.sided p-values one-sided
#' @param log.p are p-values log transformed?
#' @param inf.zscore infinite z-scores are capped to this value
#'
#' @return matrix of z-scores
#' @export
#'
#' @examples
pvalue_to_zscore <- function(mat.in, one.sided = TRUE, log.p = FALSE, inf.zscore = 16) {
    ans = mat.in
    for (col_idx in seq_len(ncol(mat.in))) {
        ans[, col_idx] = stats::qnorm(mat.in[, col_idx], lower.tail = FALSE,
                               log.p = log.p)
        ans[mat.in[, col_idx] > 1, col_idx] = 0
    }
    if (one.sided) {
        ans[ans < 0] = 0
    }
    ans[is.infinite(as.matrix(ans))] = inf.zscore

    return(ans)

}


#' Calculate hamming distance
#'
#' @param X logical matrix of calls
#'
#' @return matrix of hamming distances
#' @export
#'
#' @examples
hamming <- function(X) {
    #https://johanndejong.wordpress.com/2015/09/23/fast-hamming-distance-in-r/
    D <- (1 - X) %*% t(X)
    D + t(D)
}

#' Calculate normalized hamming distance (jaccard?)
#'
#' @param X logical matrix of calls
#'
#' @return distance object normalized by the number of number of columns
#' @export
#'
#' @examples
hamming_dist<-function(X) {
    ans = hamming(X) / ncol(X)
    return(stats::as.dist(ans))
}


#' Create a matrix for sequence to probes
#'
#' @param probe_meta data.frame with PROBE_ID and PROBE_SEQUENCE columns
#' @param debug print out debugging information
#'
#' @return a sparseMatrix that would convert a sequence matrix to probe
#' matrix upon multiplication
#' @export
#'
#' @examples
getSequenceMatToProbeMat<-function(probe_meta, debug = FALSE) {

    if (debug) {message("Generating seq_to_probe\n")}
    umeta = unique(probe_meta[,c("PROBE_ID","PROBE_SEQUENCE")])


    uniq_seq = unique(umeta$PROBE_SEQUENCE);
    probe_idx = seq_len(nrow(umeta))
    seq_idx = seq_len(length(uniq_seq))
    names(seq_idx) = uniq_seq;

    seq_idx2 = seq_idx[umeta$PROBE_SEQUENCE]

    seq_to_probe = Matrix::sparseMatrix(probe_idx, seq_idx2, x =1);
    rownames(seq_to_probe) = umeta$PROBE_ID;
    colnames(seq_to_probe) = uniq_seq;
    return(seq_to_probe);
}


#' Convert a sequence matrix to a probe matrix
#'
#' @param seq_mat matrix where the rows are peptide sequences and the columns are samples
#' @param probe_meta data.frame with the PROBE_SEQUENCE, PROBE_ID columns
#'
#' @return matrix where the rows are probe ids and the columns are samples
#' @export
#'
#' @examples
convertSequenceMatToProbeMat<-function(seq_mat, probe_meta) {
    seq_to_probe = getSequenceMatToProbeMat(probe_meta);
    probe_mat = as.matrix(seq_to_probe[,rownames(seq_mat)] %*% as.matrix(seq_mat));
    return(probe_mat)
}


getProbeMatToSequenceMat<-function(probe_meta, debug=FALSE) {
    seq_to_probe = getSequenceMatToProbeMat(probe_meta);
    if (debug) {print(utils::head(seq_to_probe))}
    probe_to_seq = Matrix::t(seq_to_probe);
    if (debug) {print(utils::head(probe_to_seq));}
    probe_to_seq = probe_to_seq / Matrix::rowSums(probe_to_seq);

    return(probe_to_seq);
}






