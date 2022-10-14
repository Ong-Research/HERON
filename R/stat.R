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
#' @export
#'
#' @examples
calcMinFDR<-function(fdrs, additional_stats=TRUE, sort=TRUE) {
    fdrs2 = as.data.frame(fdrs,stringsAsFactors=FALSE);
    ncols = ncol(fdrs);
    #cat("Calculating minFDRs\n")

    fdrs2 = t(apply(fdrs2, 1, function(l) {
        return(l[order(l,decreasing=FALSE)])
        }
        )
    )
    fdrs2 = as.data.frame(fdrs2, stringsAsFactors=FALSE)
    colnames(fdrs2) = paste0("n",seq_len(ncol(fdrs2)))
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


#' Calculate Probe-level p-values
#'
#' calculates p-values on the matrix (can be sequence as well), using either a
#' z-score global, t-stat differential, or a combination of both.
#' input a matrix that has samples on the columns and sequence probes on the
#' rows.  pData is a design matrix that describes the layout of the
#' experiment, similar to the pepStat matrix.
#'
#' @param probe_mat matrix of numeric values that where the rows are features
#' and the columns are samples. assumed to be log-transformed data
#' @param pData experiment design data.frame
#' @param t.sd_shift shift to use when calculating the differential t-test
#' p-values using multiples of the calculated standard deviation of the values.
#' Either t.sd_shift or t.abs_shift should be set
#' @param t.abs_shift shift to use when calculating the differential t-test
#' p-values using an absolute shift. Either t.sd_shift or t.abs_shift should be
#' set
#' @param z.sdshift shift to use when calculating the global z-test p-values
#' using multiples of the calculated standard deviation across all values in the
#' matrix.
#' @param use what p-value method to use, t - differential t-test, z - global
#' z-test, both - both the differential t-test and global z-test.
#' @param p.adjust.method method for adjusting p-values (multiple testing)
#'
#' @return matrix of p-values calculating on the matrix values and supplied
#' parameters
#' @export
#'
#' @examples
calcProbePValuesProbeMat<-function(
    probe_mat,
    pData,
    t.sd_shift = NA,
    t.abs_shift = NA,
    z.sdshift=0,
    use = "both",
    p.adjust.method = "BH"
) {

    if (use == "both") { use.t = TRUE; use.z = TRUE; use.c = TRUE; }
    else if (use == "t") { use.t = TRUE; use.z = FALSE; use.c = FALSE; }
    else if (use == "z") { use.t = FALSE; use.z = TRUE; use.c = FALSE; }
    else { stop("Unknown use paramater:" , use); }
    c_mat = probe_mat[,pData$TAG]

    post.cols.orig = pData$TAG[tolower(pData$visit) == "post"]
    post.cols = pData$ptid[tolower(pData$visit) == "post"]

    if (use.t) {
        pvaluet_df = calcProbePValuesTUnpaired(
            c_mat, pData, sd_shift = t.sd_shift, abs_shift = t.abs_shift);
        praw = pvaluet_df;
    }
    if (use.z) {
        pvaluez_df = calcProbePValuesZ(c_mat, pData, sd_shift = z.sdshift)
        praw = pvaluez_df;
    }
    if (use.c) {
        use_cols = intersect(colnames(pvaluet_df),colnames(pvaluez_df))
        praw = pvaluet_df[,use_cols];
        for (col in use_cols) {
            praw[,col] = pmax(praw[,col], pvaluez_df[,col]);
            praw[,col] = stats::pbeta(praw[,col], 2, 1);
            #(Wilkinson's max p-value).
        }
    }
    praw[is.na(praw)] = 1.0 # Conservatively set NAs to p-value 1
    padj_df = p_adjust_mat(praw, method = p.adjust.method);

    ans = padj_df;
    attr(ans, "pvalue") = praw;
    attr(ans, "c_mat") = c_mat;
    attr(ans, "pData") = pData;
    if (use.t) { attr(ans, "t") = pvaluet_df; }
    if (use.z) { attr(ans, "z") = pvaluez_df; }
    return(ans);
}


#' Calculate probe p-vlaues using a sequence matrix.
#'
#' @param seq_mat matrix of values where the rows are sequence identifiers
#' and the columns are samples
#' @param probe_meta data.frame with the columns PROBE_ID and PROBE_SEQUENCE
#' @param pData design matrix
#' @param t.sd_shift standard deviation shift for differential test
#' @param t.abs_shift absolute shift for differential test
#' @param z.sdshift standard deviation shift for global test
#' @param use use global-test ("z"), differential-test ("t"), or both ("both")
#' @param p.adjust.method p-value adjustment method to use (default BH)
#'
#' @return matrix of adjusted p-values with additional attributes.
#' @export
#'
#' @examples
calcProbePValuesSeqMat<-function(
        seq_mat,
        probe_meta,
        pData,
        t.sd_shift = NA,
        t.abs_shift = NA,
        z.sdshift = 0,
        use = "both",
        p.adjust.method = "BH"
) {

    seq_results = calcProbePValuesProbeMat(
        probe_mat = seq_mat,
        pData = pData,
        t.sd_shift = t.sd_shift,
        t.abs_shift = t.abs_shift,
        z.sdshift = z.sdshift,
        use = use,
        p.adjust.method = p.adjust.method
    );
    ans = convertSequenceMatToProbeMat(
        seq_results,
        probe_meta
    );
    attr(ans, "pvalue") = convertSequenceMatToProbeMat(
        attr(seq_results, "pvalue"),
        probe_meta
    )
    if ("t" %in% names(attributes(seq_results))) {
        attr(ans, "t") = convertSequenceMatToProbeMat(
            attr(seq_results, "t"),
            probe_meta
        )
    }
    if ("z" %in% names(attributes(seq_results))) {
        attr(ans, "z") = convertSequenceMatToProbeMat(
            attr(seq_results, "z"),
            probe_meta
        )
    }
    attr(ans, "seq_results") = seq_results;
    return(ans);

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
calcProbePValuesZ<-function(
        probe_mat,
        pData,
        sd_shift = 0,
        all = FALSE
) {

    if (all || missing(pData) || is.null(pData)) {
        message("No pData or all asked for, calculating on all columns");
        all_cols = colnames(probe_mat);
        post_cols = all_cols;
        post_names = all_cols;
    } else {
        pre_cols = pData$TAG[tolower(pData$visit)=="pre"]
        post_cols = pData$TAG[tolower(pData$visit) == "post"];
        post_names = pData$ptid[tolower(pData$visit) == "post"];
        all_cols = c(pre_cols, post_cols);
    }
    ans = matrix(NA, nrow = nrow(probe_mat), ncol=length(post_cols));
    vals = unlist(probe_mat[,all_cols])
    global_mean = mean(vals, na.rm=TRUE);
    global_sd = stats::sd(vals, na.rm=TRUE);
    pars = c("mean" = global_mean, "sd" = global_sd);

    post_mat = probe_mat[,post_cols];
    post_zval = (post_mat - global_mean) / global_sd;
    post_zval2 = post_zval - sd_shift;
    post_pval = apply((post_zval2), 2, stats::pnorm, lower.tail=FALSE);

    colnames(post_pval) = post_names;
    colnames(post_zval) = post_names;

    attr(post_pval,"pars") = pars;
    attr(post_pval,"zscore") = post_zval;

    return(post_pval);

}

getPairedMapping<-function(pData) {
    pre_df = pData[tolower(pData$visit) =="pre",]
    post_df = pData[tolower(pData$visit) == "post",];

    mapping = data.frame(
        ptid = pre_df$ptid,
        pre = pre_df$TAG,
        post = rep(NA, nrow(pre_df)),
        stringsAsFactors=FALSE
    );
    rownames(mapping) = mapping$ptid;
    for (idx in seq_len(nrow(post_df))) {
        post_ptid = post_df$ptid[idx];
        if (post_ptid %in% rownames(mapping)) {
            mapping[post_ptid,"post"] = post_df$TAG[idx];
        }
    }

    mapping = stats::na.omit(mapping);
    return(mapping);
}

getPTP<-function(x, stderr, sd_shift, sx, abs_shift, dfree) {
    no_shift = is.na(sd_shift) & is.na(abs_shift);
    if (no_shift) {
        tstat = (x)/stderr;
    } else if (!is.na(sd_shift)) {
        tstat = (x-sd_shift*sx)/stderr;
    } else {
        tstat = (x-abs_shift)/stderr;
    }
    ans = stats::pt(tstat, dfree, lower.tail=FALSE);
    return(ans);
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
calcProbePValuesTPaired <- function(
        probe_mat,
        pData,
        sd_shift = NA,
        abs_shift = NA,
        debug = FALSE
) {
    if (!is.na(sd_shift) && !is.na(abs_shift)) {
        stop("Either sd or abs can be set. Not both.");
    }
    mapping = getPairedMapping(pData);
    ans = matrix(NA, nrow = nrow(probe_mat), ncol=nrow(mapping))
    diff_mat = probe_mat[,mapping$post] - probe_mat[,mapping$pre];
    colnames(diff_mat) = mapping$ptid;
    rownames(diff_mat) = rownames(probe_mat);
    pars = matrix(data = NA, nrow=nrow(probe_mat), ncol = 5)
    colnames(pars) = c("diff_mean", "diff_sd", "diff_var",
        "diff_stderr", "dfree")
    rownames(pars) = rownames(probe_mat);
    for (r_idx in seq_len(nrow(probe_mat))) {
        x = t(probe_mat[r_idx,mapping$post] - probe_mat[r_idx,mapping$pre]);
        nx = length(x)
        mx = mean(x, na.rm=TRUE);
        sx = stats::sd(x, na.rm=TRUE)
        vx = stats::var(x, na.rm=TRUE);
        stderr = sqrt(vx/nx)
        dfree = sum(!is.na(x))
        pars$diff_mean[r_idx] = mx;
        pars$diff_var[r_idx] = vx;
        pars$diff_sd[r_idx] = sx;
        pars$diff_stderr[r_idx] = stderr;
        pars$dfree[r_idx] = dfree
        for (c_idx in seq_len((nrow(mapping)))) {
            tp = getPTP(x[c_idx], stderr, sd_shift, sx, abs_shift, dfree)
            ans[r_idx, c_idx] = tp
        }
    }
    ans = as.data.frame(ans, stringsAsFactors=FALSE);
    colnames(ans) = mapping$ptid;
    rownames(ans) = rownames(probe_mat);
    attr(ans, "pars") = pars;
    attr(ans, "mapping") = mapping;
    attr(ans, "diff_mat") = diff_mat;
    return(ans);
}

getPostTVal <- function(
    post_mat, pre_means,
    pre_stderr, pre_sds,
    sd_shift, abs_shift) {
    no_shift = is.na(sd_shift) && is.na(abs_shift);
    if (no_shift) {
        post_tvalues = (post_mat - pre_means)/pre_stderr;
    } else if (!is.na(sd_shift)) {
        post_tvalues = (post_mat - pre_means-sd_shift*pre_sds)/pre_stderr;
    } else {
        post_tvalues = (post_mat - pre_means-abs_shift)/pre_stderr;
    }
    return(post_tvalues);
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
calcProbePValuesTUnpaired<-function(
        probe_mat,
        pData,
        sd_shift=NA,
        abs_shift=NA
) {
    if (!is.na(sd_shift) && !is.na(abs_shift)) {
        stop("Either sd or abs can be set, not both.");
    }
    pre_cols = pData$TAG[tolower(pData$visit) =="pre"]
    post_cols = pData$TAG[tolower(pData$visit) == "post"];
    post_names = pData$ptid[tolower(pData$visit) == "post"];

    ans = matrix(NA, nrow = nrow(probe_mat),ncol=length(post_cols))
    pre_means = rowMeans(probe_mat[,pre_cols]);
    pre_sds = matrixStats::rowSds(as.matrix(probe_mat[,pre_cols]));
    post_means = rowMeans(probe_mat[,post_cols]);
    post_sds = matrixStats::rowSds(as.matrix(probe_mat[,post_cols]));
    pre_var = matrixStats::rowVars(as.matrix(probe_mat[,pre_cols]))
    n = length(pre_cols);
    pre_stderr = sqrt(pre_var / n);
    dfree = n-1;

    pars = data.frame(
        pre_mean = pre_means, pre_sds = pre_sds,
        post_mean = post_means, post_sds = post_sds,
        pre_stderr = pre_stderr, diff_mean = post_means - pre_means,
        dfree = rep(dfree, nrow(probe_mat)),
        stringsAsFactors=FALSE
    );
    rownames(pars) = rownames(probe_mat)
    post_mat = probe_mat[,post_cols];
    post_tv = getPostTVal(post_mat, pre_means, pre_stderr,
        pre_sds, sd_shift, abs_shift)
    colnames(post_tv) = post_names;
    pars = cbind(pars, post_tv)
    for (c_idx in seq_len(ncol(post_tv))) {
        ans[,c_idx] = stats::pt(q=post_tv[,c_idx], df=dfree, lower.tail=FALSE);
    }
    ans = as.data.frame(ans,stringsAsFactors=FALSE)
    rownames(ans) = rownames(probe_mat);
    colnames(ans) = post_names;
    attr(ans, "pars") = pars;
    return(ans);
}


#' Calculate protein-level p-values
#'
#' @param epitope_pvalues_mat matrix of epitope-level p-values
#' @param method meta p-value method to use
#'
#' @return matrix of protein-level p-values
#' @export
#'
#' @examples
calcProteinPValuesMat<-function(
        epitope_pvalues_mat,
        method = "maxFDR"

) {

    by_list = list(Protein = getEpitopeProtein(rownames(epitope_pvalues_mat)))

    protein_pvalues = calcMetaPValuesMat(
        pvalues_mat = epitope_pvalues_mat,
        by_list = by_list,
        method = method
    )
    return(protein_pvalues);
}



#' Calculate epitope-level p-values
#'
#' @param probe_pvalues matrix of probe p-values
#' @param epitope_ids vector of epitope ids
#' @param method meta p-value method to use
#'
#' @return matrix of epitope-level p-values
#' @export
#'
#' @examples
calcEpitopePValuesMat<-function(
    probe_pvalues,
    epitope_ids,
    method = "maxFDR"
) {

    epi_probe_df = getEpitopeIDsToProbeIDs(epitope_ids)
    #Remove any probes that are not in the probe_pvalues matrix
    epi_probe_df = epi_probe_df[epi_probe_df$PROBE_ID %in%
                                    rownames(probe_pvalues),]

    #Order the p-values to match the epi_probe_df
    probe_pvalues_mat = probe_pvalues[epi_probe_df$PROBE_ID,]

    #by List will be the epitope_id associated with the probe.
    by_list = list(EpitopeID = epi_probe_df$Epitope_ID)

    #Do the meta p-value calculation.
    epitope_pvalues = calcMetaPValuesMat(
        pvalues_mat = probe_pvalues_mat,
        by_list = by_list,
        method = method
    )
    return(epitope_pvalues);
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
    ans = pvalues_mat;

    for (col_idx in seq_len(ncol(ans))) {
        ans[,col_idx] = stats::p.adjust(ans[,col_idx], method=method);
    }
    return(ans);
}


