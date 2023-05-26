#' Calculate minimum FDRs across samples
#'
#' This code calculates the minimum FDR threshold needed
#' to call K out of N samples, where samples are in the columns
#' When additional stats is selected, also calculates the mean and median
#'
#' @param fdrs matrix of adjusted p-values
#' @param additional_stats indicator to calculate mean and median FDR
#' @param sort sort the resultant matrix by the last N column (increasing)
#'
#' @return matrix of FDR for each K level in n1..N column names
calcMinFDR<-function(fdrs, additional_stats=TRUE, sort=TRUE) {
    fdrs2 <- as.data.frame(fdrs,stringsAsFactors=FALSE)
    ncols <- ncol(fdrs)
    #cat("Calculating minFDRs\n")

    fdrs2 <- t(apply(fdrs2, 1, function(l) {
        return(l[order(l,decreasing=FALSE)])
        }
        )
    )
    fdrs2 <- as.data.frame(fdrs2, stringsAsFactors=FALSE)
    colnames(fdrs2) <- paste0("n",seq_len(ncol(fdrs2)))
    if (additional_stats) {
        fdrs2$meanFDR <- rowMeans(fdrs)
        fdrs2$medFDR <- matrixStats::rowMedians(fdrs)
    }
    n_col <- paste0("n",ncol(fdrs))
    if (sort) {
        fdrs2 <- fdrs2[order(fdrs2[,n_col],decreasing=FALSE),]
    }
    return(fdrs2)

}


#' Calculate Probe-level p-values
#'
#' calculates p-values on the matrix (can be sequence as well), using either a
#' z-score global, t-stat differential, or a combination of both.
#' input a matrix that has samples on the columns and sequence probes on the
#' rows.
#'
#' @param probe_mat numeric matrix, rows as seqs and columns as samples
#' @param pData experiment design data.frame
#' @param t.sd_shift sd shift multiplier for diff t-test
#' @param t.abs_shift abs shift to use for the diff t-test
#' @param t.paired do paired t-test
#' @param z.sdshift sd shift multiplier for global z-test
#' @param use which p-value method(s) to use
#' @param p.adjust.method method for adjusting p-values
#'
#' @return matrix of calculated p-values
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' probe_meta <- attr(heffron2021_wuhan, "probe_meta")
#' pData <- attr(heffron2021_wuhan, "pData")
#' probe_mat = convertSequenceMatToProbeMat(heffron2021_wuhan, probe_meta)
#' pval_res <- calcProbePValuesProbeMat(probe_mat, pData)
calcProbePValuesProbeMat<-function(
        probe_mat,
        pData,
        t.sd_shift = NA,
        t.abs_shift = NA,
        t.paired = FALSE,
        z.sdshift=0,
        use = "both",
        p.adjust.method = "BH") {

    if (use == "both") { use.t <- TRUE; use.z <- TRUE; use.c <- TRUE }
    else if (use == "t") { use.t <- TRUE; use.z <- FALSE; use.c <- FALSE }
    else if (use == "z") { use.t <- FALSE; use.z <- TRUE; use.c <- FALSE }
    else { stop("Unknown use paramater:" , use) }
    if (missing(pData) || is.null(pData)) {
        c_mat <- probe_mat
    } else {
        c_mat <- probe_mat[,pData$TAG]
    }
    if (use.t) {
        if (t.paired) {
            pvaluet_df <- calcProbePValuesTPaired(
                probe_mat = c_mat, pData = pData,
                sd_shift = t.sd_shift, abs_shift = t.abs_shift)
        } else {
            pvaluet_df <- calcProbePValuesTUnpaired(
                c_mat, pData, sd_shift = t.sd_shift, abs_shift = t.abs_shift)
        }
        praw <- pvaluet_df
    }
    if (use.z) {
        pvaluez_df <- calcProbePValuesZ(c_mat, pData, sd_shift = z.sdshift)
        praw <- pvaluez_df
    }
    if (use.c) {
        praw <- combinePValueMatrix(pvaluet_df, pvaluez_df)
    }
    praw[is.na(praw)] <- 1.0 # Conservatively set NAs to p-value = 1
    padj_df <- p_adjust_mat(praw, method = p.adjust.method)

    ans <- padj_df
    attr(ans, "pvalue") <- praw
    attr(ans, "c_mat") <- c_mat
    attr(ans, "pData") <- pData
    if (use.t) { attr(ans, "t") <- pvaluet_df }
    if (use.z) { attr(ans, "z") <- pvaluez_df }
    return(ans)
}

combinePValueMatrix<-function(pmat1, pmat2) {
    use_cols <- intersect(colnames(pmat1), colnames(pmat2))
    if (length(use_cols) != ncol(pmat1)) {
        warning("Combining p-values, some columns are mismatched")
    }
    ans <- pmat1[,use_cols]
    for (col in use_cols) {
        ans[,col] <- pmax(ans[,col], pmat2[,col])
        ##(Wilkinson's max p-value)
        ans[,col] <- stats::pbeta(ans[,col], 2, 1)
    }
    return(ans)
}

#' Calculate probe p-values using a sequence matrix.
#'
#' @param seq_mat matrix of values where the rows are sequence identifiers
#' and the columns are samples
#' @param probe_meta data.frame with the columns PROBE_ID and PROBE_SEQUENCE
#' @param pData design matrix
#' @param t.sd_shift standard deviation shift for differential test
#' @param t.abs_shift absolute shift for differential test
#' @param t.paired run paired analysis
#' @param z.sdshift standard deviation shift for global test
#' @param use use global-test ("z"), differential-test ("t"), or both ("both")
#' @param p.adjust.method p-value adjustment method to use (default BH)
#'
#' @return matrix of adjusted p-values with additional attributes.
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' probe_meta <- attr(heffron2021_wuhan, "probe_meta")
#' pData <- attr(heffron2021_wuhan, "pData")
#' pval_res <- calcProbePValuesSeqMat(heffron2021_wuhan, probe_meta, pData)
calcProbePValuesSeqMat<-function(
        seq_mat,
        probe_meta,
        pData,
        t.sd_shift = NA,
        t.abs_shift = NA,
        t.paired = FALSE,
        z.sdshift = 0,
        use = "both",
        p.adjust.method = "BH"
) {

    seq_results <- calcProbePValuesProbeMat(
        probe_mat = seq_mat,
        pData = pData,
        t.sd_shift = t.sd_shift,
        t.abs_shift = t.abs_shift,
        t.paired = t.paired,
        z.sdshift = z.sdshift,
        use = use,
        p.adjust.method = p.adjust.method
    )
    ans <- convertSequenceMatToProbeMat(
        seq_results,
        probe_meta
    )
    attr(ans, "pvalue") <- convertSequenceMatToProbeMat(
        attr(seq_results, "pvalue"),
        probe_meta
    )
    if ("t" %in% names(attributes(seq_results))) {
        attr(ans, "t") <- convertSequenceMatToProbeMat(
            attr(seq_results, "t"),
            probe_meta
        )
    }
    if ("z" %in% names(attributes(seq_results))) {
        attr(ans, "z") <- convertSequenceMatToProbeMat(
            attr(seq_results, "z"),
            probe_meta
        )
    }
    attr(ans, "seq_results") <- seq_results
    return(ans)

}

#' Calculate p-values using the "exprs" assay from the sequence or probe dataset
#'
#' @param obj HERONSequenceDataSet or HERONProbeDataSet
#' @param t.sd_shift
#' @param t.abs_shift
#' @param t.paired
#' @param z.sdshift
#' @param use
#' @param p.adjust.method
#'
#' @return HERONSequenceDataSet/HERONProbeDataSet with the pvalue assay added
#' @export
#'
#' @examples
calcCombPValues<-function(
    obj,
    t.sd_shift = NA,
    t.abs_shift = NA,
    t.paired = FALSE,
    z.sdshift = 0,
    use = "both",
    p.adjust.method = "BH"
) {
    stopifnot(is(obj, "HERONSequenceDataSet") || is(obj, "HERONProbeDataSet"))
    pval <- calcProbePValuesProbeMat(
        probe_mat = assays(obj)$expr,
        pData = colData(obj),
        t.sd_shift = t.sd_shift,
        t.abs_shift = t.abs_shift,
        t.paired = t.paired,
        z.sdshift = z.sdshift,
        use = use,
        p.adjust.method = "none"
    )

    res <- addPValues(obj, pval)
    res <- p_adjust_ds(res)
    return(res)
}

addPValues<-function(obj, pval) {
    res <- obj
    if (all(dim(assay(obj, "exprs")) == dim(pval)) &&
        all(colnames(assay(obj,"exprs")) == colnames(pval)) &&
        all(rownames(assay(obj,"exprs")) == rownames(pval))
    ) {
        assay(res, "pvalue") <- pval
    } else {
        message("Reduced rows/columns")
        exprs_old <- assay(obj, "exprs")
        colData_old <- colData(obj)

        exprs_new <- exprs_old[rownames(pval), colnames(pval)]
        colData_new <- colData_old[(colnames(pval)),]
        if (is(obj, "HERONSequenceDataSet")) {
            res <- HERONSequenceDataSet(
                exprs = exprs_new,
                colData = colData_new
            )
        } else if (is(obj, "HERONProbeDataSet")) {
            rowRanges_new <- rowRanges(obj)
            res <- HERONProbeDataSet(
                assays = list(exprs = exprs_new),
                colData = colData_new,
                rowRanges = rowRanges_new
            )
        }
        assay(res, "pvalue") <- pval
    }
    return(res)
}

#' Calculate Global p-values Using Normal (z) Distribution
#'
#' @param probe_mat numeric matrix or data.frame of values
#' @param pData design data.frame
#' @param sd_shift standard deviation shift to use when calculating p-values
#' (default 0)
#' @param all calculate p-values for all samples, otherwise report the p-values
#' for just the post samples
#'
#' @return matrix of "p-values"
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' probe_meta <- attr(heffron2021_wuhan, "probe_meta")
#' pData <- attr(heffron2021_wuhan, "pData")
#' pval_res <- calcProbePValuesZ(heffron2021_wuhan, pData)
calcProbePValuesZ<-function(
        probe_mat,
        pData,
        sd_shift = 0,
        all = FALSE
) {

    if (all || missing(pData) || is.null(pData)) {
        message("No pData or all asked for, calculating on all columns")
        all_cols <- colnames(probe_mat)
        post_cols <- all_cols
        post_names <- all_cols
    } else {
        pre_cols <- pData$TAG[tolower(pData$visit)=="pre"]
        post_cols <- pData$TAG[tolower(pData$visit) == "post"]
        post_names <- pData$ptid[tolower(pData$visit) == "post"]
        all_cols <- post_cols
    }

    vals <- unlist(probe_mat[,all_cols])
    g_mean <- mean(vals, na.rm=TRUE)
    g_sd <- stats::sd(vals, na.rm=TRUE)

    zvalues <- matrix(NA, nrow = nrow(probe_mat), ncol=ncol(probe_mat))
    rownames(zvalues) <- rownames(probe_mat)
    colnames(zvalues) <- colnames(probe_mat)
    zvalues[,post_cols] <- as.matrix(probe_mat[,post_cols] - g_mean) / g_sd

    pvalues <- apply(zvalues - sd_shift, 2, stats::pnorm, lower.tail=FALSE)

    pars <- c("mean" = g_mean, "sd" = g_sd)
    attr(pvalues,"pars") <- pars
    attr(pvalues,"zscore") <- zvalues

    return(pvalues)
}


#' Title
#'
#' @param pData
#'
#' @return
#' @export
#'
#' @examples
getPairedMapping<-function(pData) {
    pre_df <- pData[tolower(pData$visit) =="pre",]
    post_df <- pData[tolower(pData$visit) == "post",]

    mapping <- data.frame(
        ptid = pre_df$ptid,
        pre = pre_df$TAG,
        post = rep(NA, nrow(pre_df)),
        stringsAsFactors=FALSE
    )
    rownames(mapping) <- mapping$ptid
    for (idx in seq_len(nrow(post_df))) {
        post_ptid <- post_df$ptid[idx]
        if (post_ptid %in% rownames(mapping)) {
            mapping[post_ptid,"post"] <- post_df$TAG[idx]
        }
    }

    mapping <- stats::na.omit(mapping)
    return(mapping)
}

getPTP<-function(x, stderr, sd_shift, sx, abs_shift, dfree) {
    no_shift <- is.na(sd_shift) & is.na(abs_shift)
    if (no_shift) {
        tstat <- (x)/stderr
    } else if (!is.na(sd_shift)) {
        tstat <- (x-sd_shift*sx)/stderr
    } else {
        tstat <- (x-abs_shift)/stderr
    }
    ans <- stats::pt(tstat, dfree, lower.tail=FALSE)
    return(ans)
}


#' Calculate Probe p-values using a differential paired t-test
#'
#' @param probe_mat numeric matrix or data.frame of values
#' @param pData design data.frame
#' @param sd_shift standard deviation shift to use when calculating p-values.
#' Either sd_shift or abs_shift should be set
#' @param abs_shift absolute shift to use when calculating p-values.
#' @param debug print debugging information
#'
#' @return matrix of p-values on the post columns defined in the pData matrix.
#' Attributes of the matrix are:
#'
#' pars - data.frame parameters used in the paired t-test for each row
#' (e.g. df, sd)
#'
#' mapping - data.frame of mapping used for pre-post column calculation
#' diff_mat - data.frame
#' containing the post-pre differences for each sample (column) and probe (row)
#'
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' probe_meta <- attr(heffron2021_wuhan, "probe_meta")
#' pData <- attr(heffron2021_wuhan, "pData")
#' pre_idx = which(pData$visit == "pre")
#' ## Make some samples paired
#' pData_post = pData[pData$visit == "post",]
#' new_ids = pData_post$Sample_ID[seq_len(5)]
#' pData$ptid[pre_idx[seq_len(5)]] = new_ids
#' pval_res <- calcProbePValuesTPaired(heffron2021_wuhan, pData)
calcProbePValuesTPaired <- function(
        probe_mat,
        pData,
        sd_shift = NA,
        abs_shift = NA,
        debug = FALSE
) {
    if (!is.na(sd_shift) && !is.na(abs_shift)) {
        stop("Either sd or abs can be set. Not both.")
    }
    mapping <- getPairedMapping(pData)
    ans <- matrix(NA, nrow = nrow(probe_mat), ncol=nrow(mapping))
    diff_mat <- probe_mat[,mapping$post] - probe_mat[,mapping$pre]
    colnames(diff_mat) <- mapping$ptid
    rownames(diff_mat) <- rownames(probe_mat)
    pars <- matrix(data = NA, nrow=nrow(probe_mat), ncol = 5)
    colnames(pars) <- c("diff_mean", "diff_sd", "diff_var",
        "diff_stderr", "dfree")
    rownames(pars) <- rownames(probe_mat)
    for (r_idx in seq_len(nrow(probe_mat))) {
        x <- c(t(probe_mat[r_idx,mapping$post] - probe_mat[r_idx,mapping$pre]))
        nx <- length(x)
        mx <- mean(x, na.rm=TRUE)
        sx <- stats::sd(x, na.rm=TRUE)
        vx <- stats::var(x, na.rm=TRUE)
        stderr <- sqrt(vx/nx)
        dfree <- sum(!is.na(x))
        pars[r_idx, "diff_mean"] <- mx
        pars[r_idx, "diff_var"] <- vx
        pars[r_idx, "diff_sd"] <- sx
        pars[r_idx, "diff_stderr"] <- stderr
        pars[r_idx, "dfree"] <- dfree
        for (c_idx in seq_len((nrow(mapping)))) {
            tp <- getPTP(x[c_idx], stderr, sd_shift, sx, abs_shift, dfree)
            ans[r_idx, c_idx] <- tp
        }
    }
    ans <- as.data.frame(ans, stringsAsFactors=FALSE)
    colnames(ans) <- mapping$ptid
    rownames(ans) <- rownames(probe_mat)
    attr(ans, "pars") <- pars
    attr(ans, "mapping") <- mapping
    attr(ans, "diff_mat") <- diff_mat
    return(ans)
}

getPostTVal <- function(
    post_mat, pre_means,
    pre_stderr, pre_sds,
    sd_shift, abs_shift) {
    no_shift <- is.na(sd_shift) && is.na(abs_shift)
    if (no_shift) {
        post_tvalues <- (post_mat - pre_means)/pre_stderr
    } else if (!is.na(sd_shift)) {
        post_tvalues <- (post_mat - pre_means-sd_shift*pre_sds)/pre_stderr
    } else {
        post_tvalues <- (post_mat - pre_means-abs_shift)/pre_stderr
    }
    return(post_tvalues)
}

#' Calculate Probe p-values using a differential unpaired t-test
#'
#' @param probe_mat numeric matrix or data.frame of values
#' @param pData design data.frame
#' @param sd_shift standard deviation shift to use when calculating p-values
#' Either sd_shift or abs_shift should be set
#' @param abs_shift absolute shift to use when calculating p-values
#'
#' @return matrix of p-values on the post columns defined in the pData matrix
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' probe_meta <- attr(heffron2021_wuhan, "probe_meta")
#' pData <- attr(heffron2021_wuhan, "pData")
#' pval_res <- calcProbePValuesTUnpaired(heffron2021_wuhan, pData)
calcProbePValuesTUnpaired<-function(
        probe_mat,
        pData,
        sd_shift=NA,
        abs_shift=NA,
        keep_all_cols = FALSE
) {
    if (!is.na(sd_shift) && !is.na(abs_shift)) {
        stop("Either sd or abs can be set, not both.")
    }
    pre_cols <- pData$TAG[tolower(pData$visit) =="pre"]
    post_cols <- pData$TAG[tolower(pData$visit) == "post"]
    post_names <- pData$ptid[tolower(pData$visit) == "post"]

    ans <- matrix(NA, nrow = nrow(probe_mat),ncol=ncol(probe_mat))
    rownames(ans) <- rownames(probe_mat)
    colnames(ans) <- colnames(probe_mat)
    pre_means <- rowMeans(probe_mat[,pre_cols])
    pre_sds <- matrixStats::rowSds(as.matrix(probe_mat[,pre_cols]))
    post_means <- rowMeans(probe_mat[,post_cols])
    post_sds <- matrixStats::rowSds(as.matrix(probe_mat[,post_cols]))
    pre_var <- matrixStats::rowVars(as.matrix(probe_mat[,pre_cols]))
    n <- length(pre_cols)
    pre_stderr <- sqrt(pre_var / n)
    dfree <- n-1

    pars <- data.frame(
        pre_mean = pre_means, pre_sds = pre_sds,
        post_mean = post_means, post_sds = post_sds,
        pre_stderr = pre_stderr, diff_mean = post_means - pre_means,
        dfree = rep(dfree, nrow(probe_mat)),
        stringsAsFactors=FALSE
    )
    rownames(pars) <- rownames(probe_mat)
    post_mat <- probe_mat[,post_cols]
    post_tv <- getPostTVal(post_mat, pre_means, pre_stderr,
        pre_sds, sd_shift, abs_shift)
    colnames(post_tv) <- post_names
    pars <- cbind(pars, post_tv)
    for (post_col in post_cols) {
        ans[,post_col] <-
            stats::pt(q=post_tv[,post_col], df=dfree, lower.tail=FALSE)
    }
    ans <- as.data.frame(ans,stringsAsFactors=FALSE)
    if (!keep_all_cols) {ans<-ans[,post_cols]}

    attr(ans, "pars") <- pars
    return(ans)
}


#' Calculate protein-level p-values
#'
#' @param epitope_pvalues_mat matrix of epitope-level p-values
#' @param metap_method meta p-value method to use
#' @param p_adjust_method p.adjust method to use
#'
#' @return matrix of protein-level p-values
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' probe_meta <- attr(heffron2021_wuhan, "probe_meta")
#' pData <- attr(heffron2021_wuhan, "pData")
#' pval_res <- calcProbePValuesSeqMat(heffron2021_wuhan, probe_meta, pData)
#' calls_res <- makeProbeCalls(pval_res)
#' segments_res <- findEpitopeSegments(probe_calls = calls_res)
#' epval_res <- calcEpitopePValuesMat(pval_res, segments_res)
#' ppval_res <- calcProteinPValuesMat(epval_res)
calcProteinPValuesMat<-function(
        epitope_pvalues_mat,
        metap_method = "wmin1",
        p_adjust_method = "BH"

) {
    if ("pvalue" %in% names(attributes(epitope_pvalues_mat))) {
        pvalues <- attr(epitope_pvalues_mat, "pvalue")
    } else {
        pvalues <- epitope_pvalues_mat
    }

    by_list <- list(Protein = getEpitopeProtein(rownames(pvalues)))

    protein_pvalues <- calcMetaPValuesMat(
        pvalues_mat = pvalues,
        by_list = by_list,
        method = metap_method
    )
    ans <- protein_pvalues
    ans <- p_adjust_mat(protein_pvalues, p_adjust_method)
    attr(ans, "pvalue") <- protein_pvalues
    return(ans)
}

#' Title
#'
#' @param epitope_ds
#' @param metap_method
#' @param p_adjust_method
#'
#' @return
#' @export
#'
#' @examples
calcProteinPValuesEpitopeDS<-function(
        epitope_ds,
        metap_method = "wmin1",
        p_adjust_method = "BH"

) {
    stopifnot(is(epitope_ds, "HERONEpitopeDataSet"))

    pvalues <- assay(epitope_ds, "pvalue")

    by_list <- list(Protein = getEpitopeProtein(rownames(pvalues)))

    protein_pvalues <- calcMetaPValuesMat(
        pvalues_mat = pvalues,
        by_list = by_list,
        method = metap_method
    )
    res <- HERONProteinDataSet(pvalue = protein_pvalues)
    colData(res) <- colData(epitope_ds)
    res <- p_adjust_ds(res, p_adjust_method)
    return(res)
}


#' Calculate epitope-level p-values
#'
#' @param probe_pvalues matrix of probe p-values
#' @param epitope_ids vector of epitope ids
#' @param metap_method meta p-value method to use
#' @param p_adjust_method what p.adjust method to use.
#'
#' @return matrix of epitope-level p-values
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' probe_meta <- attr(heffron2021_wuhan, "probe_meta")
#' pData <- attr(heffron2021_wuhan, "pData")
#' pval_res <- calcProbePValuesSeqMat(heffron2021_wuhan, probe_meta, pData)
#' calls_res <- makeProbeCalls(pval_res)
#' segments_res <- findEpitopeSegments(probe_calls = calls_res)
#' epval_res <- calcEpitopePValuesMat(attr(pval_res, "pvalue"), segments_res)
calcEpitopePValuesMat<-function(
    probe_pvalues,
    epitope_ids,
    metap_method = "wmax1",
    p_adjust_method = "BH"
) {
    if ("pvalue" %in% names(attributes(probe_pvalues))) {
        pvalues <- attr(probe_pvalues, "pvalue")
    } else {
        pvalues <- probe_pvalues
    }

    epi_probe_df <- getEpitopeIDsToProbeIDs(epitope_ids)
    #Remove any probes that are not in the probe_pvalues matrix
    epi_probe_df <- epi_probe_df[epi_probe_df$PROBE_ID %in%
                                    rownames(probe_pvalues),]

    #Order the p-values to match the epi_probe_df
    probe_pvalues_mat <- pvalues[epi_probe_df$PROBE_ID,]

    #by List will be the epitope_id associated with the probe.
    by_list <- list(EpitopeID = epi_probe_df$Epitope_ID)

    #Do the meta p-value calculation.
    epitope_pvalues <- calcMetaPValuesMat(
        pvalues_mat = probe_pvalues_mat,
        by_list = by_list,
        method = metap_method
    )

    ans <- epitope_pvalues
    ans <- p_adjust_mat(epitope_pvalues, p_adjust_method)
    attr(ans, "pvalue") <- epitope_pvalues
    return(ans)
}

#' Calculate epitope-level p-values
#'
#' @param probe_pds
#' @param epitope_ids vector of epitope ids
#' @param metap_method meta p-value method to use
#' @param p_adjust_method what p.adjust method to use.
#'
#' @return HERONEpitopeDataSet
#' @export
#'
#' @examples
calcEpitopePValuesProbeDS<-function(
        probe_pds,
        epitope_ids,
        metap_method = "wmax1",
        p_adjust_method = "BH"
) {
    stopifnot(is(probe_pds, "HERONProbeDataSet"))
    stopifnot("pvalue" %in% assayNames(probe_pds))
    pvalues_mat <- calcEpitopePValuesMat(
        probe_pvalues = assay(probe_pds, "pvalue"),
        epitope_ids = epitope_ids,
        metap_method = metap_method,
        p_adjust_method = p_adjust_method
    )

    eds <- HERONEpitopeDataSet(pvalue = pvalues_mat)
    colData(eds) <- colData(probe_pds)
    assays(eds)$padj <- p_adjust_mat(pvalues_mat, method = p_adjust_method)
    return(eds)
}

#' Adjust a matrix of p-values column-by-column
#'
#' @param pvalues_mat matrix of p-values
#' @param method what adjustment algorithm to use (see p.adjust)
#'
#' @return matrix of column-by-column adjusted p-values
#' @export
#'
#' @examples
#' mat = matrix(runif(25) >= 0.5, nrow=5)
#' rownames(mat) = paste0("A;",seq_len(nrow(mat)))
#' p_adjust_mat(mat)
p_adjust_mat<-function(pvalues_mat, method = "BH") {
    ans <- pvalues_mat
    for (col_idx in seq_len(ncol(ans))) {
        ans[,col_idx] <- stats::p.adjust(ans[,col_idx], method=method)
    }
    return(ans)
}

#' Adjust an assay of p-values column-by-column
#'
#' @param obj SummarizedExperiment with a "pvalue" assay
#' @param method what adjustment algorithm to use (see p.adjust)
#'
#' @return
#' @export
#'
#' @examples
p_adjust_ds <- function(obj, method = "BH") {
    stopifnot(is(obj, "SummarizedExperiment"))
    stopifnot("pvalue" %in% assayNames(obj))
    assay(obj, "padj") <- p_adjust_mat(assay(obj, "pvalue"), method = method)
    return(obj)
}

