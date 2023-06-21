
#' Get the amino-acid starting position of the probe within the protein.
#'
#' @param probes vector of probes (i.e. c("A;1", "A;2"))
#'
#' @return starting locations of the probes with their associated proteins
#' @export
#'
#' @examples
#' getProteinStart("A;1")
#' getProteinStart("B;2")
#' getProteinStart(c("A;1","B;2"))
getProteinStart<-function(probes) {
    probe.list <- strsplit(probes, ";")
    protein.start.list <- lapply(
        probe.list,
        function(x) {

            return(x[length(x)])
        }
    )
    return(as.integer(unlist(protein.start.list)))
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
#' getProteinLabel("A;1")
#' getProteinLabel("B;2")
#' getProteinLabel(c("A;1","B;2"))
getProteinLabel<-function(probes) {
    probe.list <- strsplit(probes,";")
    protein.list <- lapply(
        probe.list, function(x){
            n <- length(x)-1
            return(paste(x[seq_len(n)],sep=";",collapse=";"))
        }
    )
    protein.labels <- unlist(protein.list)
    return(protein.labels)
}

#' Get Protein Tiling
#'
#' Given a set of probes, estimate the tiling of the probes
#' across the protein.
#' Usually, you will want to calculate this on all the probes available in the
#' dataset.
#'
#' @param probes vector of probes (i.e. A;1, A;2)
#' @param return.vector Return result as vector or return as data.frame
#'
#' @return For each protein, the estimating tiling (spacing) of the probes
#' across the amino acid sequence.
#' @export
#'
#' @examples
#' getProteinTiling(c("A;1","A;2","A;3", "B;2","B;3", "C;1", "C;3"))
getProteinTiling<-function(probes, return.vector=TRUE) {
    Pos <- getProteinStart(probes)
    Protein <- getProteinLabel(probes)
    tiling.df <- stats::aggregate(
        as.integer(Pos),
        by = list(Protein = Protein),
        function(l) {
            if (length(l) > 1) {
                l <- l[order(l, decreasing = FALSE)]
                return(l[2] - l[1]) #Assume that the tiling is consistent across
                #the probes.
            }
            #Not enough probes for the protein to calculate a tiling.
            return(-1)
        }
    )
    ans <- tiling.df
    colnames(ans)[2] <- "Tiling"
    if (return.vector) {
        ans.vec <- ans$Tiling
        names(ans.vec) <- ans$Protein
        ans <- ans.vec
    }
    return(ans)
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
#' min_max(10, 1, 5)
min_max<-function(val, min.value, max.value) {
    val[val < min.value] <- min.value
    val[val > max.value] <- max.value
    return(val)
}

#' Concatenate sequences together based upon their start positions.
#' Assumes the probe sequences have an overlap.
#' @param positions start positions of probes in protein
#' @param sequences probe sequences of probes
#'
#' @return concatenated sequence (character)
#' @export
#'
#' @examples
#' positions <- c(1,2)
#' sequences <- c("MSGSASFEGGVFSPYL", "SGSASFEGGVFSPYLT")
#' catSequences(positions, sequences)
catSequences <- function (positions, sequences) {
    if (length(sequences) == 1) {
        return(sequences)
    }
    positions <- positions - positions[1] + 1
    max_seq_length <- max(positions) + max(nchar(sequences))
    seq <- rep("", max_seq_length)
    for (idx in seq_len(length(positions))) {
        pos <- positions[idx]
        aa_vec <- strsplit(sequences[idx], "")[[1]]
        startp <- pos
        endp <- pos + length(aa_vec) - 1
        seq[seq.int(startp, endp)] <- aa_vec
    }
    seq <- seq[seq != ""]
    seqc <- paste0(seq, collapse = "", sep = "")
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
#' mat <- matrix(runif(100), nrow=10)
#' rownames(mat) <- paste0("A;",seq_len(nrow(mat)))
#' pvalue_to_zscore(mat)
pvalue_to_zscore <- function(
    mat.in,
    one.sided = TRUE,
    log.p = FALSE,
    inf.zscore = 16
) {
    ans <- mat.in
    for (col_idx in seq_len(ncol(mat.in))) {
        ans[, col_idx] <- stats::qnorm(mat.in[, col_idx], lower.tail = FALSE,
            log.p = log.p)
        #pvalue of 1 => z-score of -inf.zscore
        ans[mat.in[, col_idx] >= 1, col_idx] <- -inf.zscore
    }

    ans[is.infinite(as.matrix(ans)) & ans > 0] <- inf.zscore
    ans[is.infinite(as.matrix(ans)) & ans < 0] <- -inf.zscore

    if (one.sided) {
        ans[ans < 0] <- 0
    }

    return(ans)
}

#' Calculate hamming distance
#'
#' @param X logical matrix of calls
#'
#' @return matrix of hamming distances
#'
#' @examples
#' mat = matrix(runif(100) >= 0.5, nrow=10)
#' rownames(mat) = paste0("A;",seq_len(nrow(mat)))
#' hamming(mat)
#' @noRd
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
#'
#' @examples
#' mat = matrix(runif(100) >= 0.5, nrow=10)
#' rownames(mat) = paste0("A;",seq_len(nrow(mat)))
#' hamming_dist(mat)
#' @noRd
hamming_dist<-function(X) {
    ans <- hamming(X) / ncol(X)
    return(stats::as.dist(ans))
}

#' Convert a sequence matrix to a probe matrix
#'
#' @param seq_mat matrix with rows as sequences and the columns as samples
#' @param probe_meta data.frame with the PROBE_SEQUENCE, PROBE_ID columns
#'
#' @return matrix where the rows are probe ids and the columns are samples
#'
#' @noRd
#'
#' @examples
#' data(heffron2021_wuhan)
#' seq_mat <- assay(heffron2021_wuhan, "exprs")
#' pr_meta <- metadata(heffron2021_wuhan)$probe_meta
#' probe_mat = convertSequenceMatToProbeMat(seq_mat, pr_meta)
convertSequenceMatToProbeMat<-function(seq_mat, probe_meta) {

    umeta <- unique(probe_meta[,c("PROBE_SEQUENCE", "PROBE_ID")])
    ans <- merge(umeta, seq_mat, by.x="PROBE_SEQUENCE", by.y = 0)
    rownames(ans) <- ans$PROBE_ID
    ans <- ans[,c(-1,-2), drop=FALSE] #Make sure we maintain the data.frame
    return(ans)
}

#' Convert HERONSequenceDataSet to HERONProbeDataSet
#'
#' @param seq_ds a HERONSequenceDataSet object
#' @param probe_meta optional data.frame with the PROBE_SEQUENCE, PROBE_ID
#' columns
#'
#' the probe meta data frame can be provided within the metadata()$probe_meta
#' or as a argument to the function.  The argument supersedes the metadata
#' list.
#'
#' @return HERONProbeDataSet
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' probe_ds <- convertSequenceDSToProbeDS(heffron2021_wuhan)
#' probe_meta <- metadata(heffron2021_wuhan)$probe_meta
#' probe_ds <- convertSequenceDSToProbeDS(heffron2021_wuhan, probe_meta)
convertSequenceDSToProbeDS<-function(seq_ds, probe_meta) {
    stopifnot(is(seq_ds, "HERONSequenceDataSet"))
    if (missing(probe_meta)) {
        if ("probe_meta" %in% names(metadata(seq_ds))) {
            probe_meta <- metadata(seq_ds)$probe_meta
        } else {
            stop("Missing probe meta data frame. Either put in ",
            "metadata()$probe_meta or pass in the probe_meta argument")
        }
    }
    umeta <- unique(probe_meta[,c("PROBE_SEQUENCE", "PROBE_ID")])
    passay_list <- lapply(
        assays(seq_ds),
        function(a) {
            a <- convertSequenceMatToProbeMat(a, umeta)
            return(a[umeta$PROBE_ID,])
        }
    )
    probe_ds <- HERONProbeDataSet(
        assays = passay_list,
        colData = colData(seq_ds),
        rowRanges = getGRanges(umeta)
    )
    return(probe_ds)
}

#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges GRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges IRanges
getGRanges<-function(umeta) {
    ans <- GRanges(
        seqnames = getProteinLabel(umeta$PROBE_ID),
        ranges = IRanges(
            start = getProteinStart(umeta$PROBE_ID),
            width = nchar(umeta$PROBE_SEQUENCE),
            names = umeta$PROBE_ID,
            PROBE_ID = umeta$PROBE_ID,
            PROBE_SEQUENCE = umeta$PROBE_SEQUENCE
        )
    )
    return(ans)
}


toNumericMatrix<-function(in_obj) {
    in_obj <- as.matrix(in_obj)

    for (col_idx in  seq_len(ncol(in_obj))) {
        in_obj[,col_idx] <- as.numeric(in_obj[,col_idx])
    }
    return(in_obj)
}



