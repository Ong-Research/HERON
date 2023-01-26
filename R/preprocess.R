


#' Normalize a matrix/data.frame using quantile normalization
#'
#' @param in_mat matrix or data.frame of numeric values to be normalized
#'
#' @return normalized data.frame
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' seq_mat_qn = quantileNormalize(heffron2021_wuhan)
#'
quantileNormalize<-function(in_mat) {
    #Use limma's version of quantile normalization
    norm_mat <- limma::normalizeQuantiles(as.matrix(in_mat));
    norm_mat <- as.data.frame(norm_mat);
    colnames(norm_mat) <- colnames(in_mat);
    rownames(norm_mat) <- rownames(in_mat);
    return(norm_mat);
}

#' Smooth probes across protein tiling
#'
#' @param probe_mat probe matrix
#' @param w smoothing width, probes w/2 before and after are used for smoothing
#' @param eps error tolerance for detecting 0 or 1 in smoothing code
#'
#' @return probe matrix with
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' probe_meta <- attr(heffron2021_wuhan, "probe_meta")
#' probe_mat <- convertSequenceMatToProbeMat(heffron2021_wuhan, probe_meta)
#' smoothed_mat <- smoothProbeMat(probe_mat)
smoothProbeMat<-function(probe_mat, w=2, eps = 1e-6) {
    if (w %% 2 != 0) {stop("w needs to be an even number!")}
    smoothing_mat <- createProbeSmoothMat(probes = rownames(probe_mat), w=w);
    ans <- smoothProbeMatInternal(probe_mat, smoothing_mat, eps);
    return(ans);
}

smoothProbeMatInternal<-function(probe_mat, smoothing_mat, eps=1e-6) {
    probe_mat2 <- as.matrix(probe_mat);
    NAs <- is.na(probe_mat2);

    smoothing_mat2 <- smoothing_mat[rownames(probe_mat2), rownames(probe_mat2)];
    if (sum(NAs) == 0) { ## We can use the matrix as-is if no NAs exist.
        probe_smooth_mat <- as.matrix(
            smoothing_mat2 %*% probe_mat2
        );
    } else {
        probe_smooth_mat <- probe_mat2;
        for (col_idx in seq_len(ncol(probe_mat2))) {
            values <- probe_mat2[,col_idx];
            svalues <- NULL;
            vnas <- is.na(values);
            smoothing_mat3 <- smoothing_mat2;
            ##Set the probes with the missing values to zero, then re-normalize
            ##the weights in the matrix. Then multiply the new smoothing matrix
            ## by the vector of with missing values set to 0. Then set missing
            ## value back into result. Copy results back to result matrix.
            if (sum(vnas) > 0) {
                smoothing_mat3[,vnas] <- 0;
                smoothing_mat3[vnas,] <- 0;
                rm <- Matrix::rowSums(smoothing_mat3);
                rind <- which(rm >= eps & rm <= 1-eps) #Only modify if needed
                smoothing_mat3[rind,] <- smoothing_mat3[rind,] / rm[rind]
                values[vnas] <- 0;
                svalues <- smoothing_mat3 %*% values;
                svalues[vnas && rm < eps] <- NA; #These values *should* be NA.
            } else {
                svalues <- smoothing_mat3 %*% values;
            }
            probe_smooth_mat[,col_idx] <- svalues[,1]
        }
    }
    return(probe_smooth_mat);
}

createProbeSmoothMatProtein<-function(protein_meta, lr) {
    sub_meta <- protein_meta[order(protein_meta$POSITION, decreasing=FALSE),];
    sub_probe_mat <- Matrix::Diagonal(nrow(sub_meta), x=0);
    rownames(sub_probe_mat) <- sub_meta$PROBE_ID;
    colnames(sub_probe_mat) <- sub_meta$PROBE_ID;
    for (idx1 in seq_len(nrow(sub_meta))) {
        probe <- sub_meta$PROBE_ID[idx1];
        pos <- sub_meta$POSITION[idx1];
        pos_left <- pos - lr;
        pos_right <- pos + lr;
        idx2 <- sub_meta$POSITION >= pos_left & sub_meta$POSITION <= pos_right
        probes <- sub_meta$PROBE_ID[idx2];
        weights <- rep(1.0,length(probes))
        sub_probe_mat[probe,probes] <- weights / sum(weights);
    }
    return(sub_probe_mat);
}


createProbeSmoothMat<-function(probes, w=2) {
    umeta <- data.frame(
        PROBE_ID = probes,
        SEQ_ID = getProteinLabel(probes),
        POSITION = getProteinStart(probes),
        stringsAsFactors = FALSE
    )
    if (w <=0) {
        probe_probe_mat <- Matrix::Diagonal(nrow(umeta), x=1);
        rownames(probe_probe_mat) <- umeta$PROBE_ID;
        colnames(probe_probe_mat) <- umeta$PROBE_ID;
    } else {
        lr <- w/2;
        split_list <- split(umeta, umeta$SEQ_ID);
        probe_probe_mat_list <- lapply(
            split_list,
            createProbeSmoothMatProtein,
            lr = lr
        )
        probe_probe_mat <- Matrix::bdiag(probe_probe_mat_list);
        cols <- unlist(lapply(probe_probe_mat_list, colnames))
        colnames(probe_probe_mat) <- cols;
        rownames(probe_probe_mat) <- cols;
    }
    return(probe_probe_mat);
}
