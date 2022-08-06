#' Calculate minimum FDRs across samples
#'
#' This code calculates the minimum FDR threshold needed to be a hit in K out of N samples
#' Where samples are in the columns
#' When additional stats is selected, also calculates the mean and median
#'
#' @param fdrs matrix of adjusted p-values
#' @param additional_stats indicator of additional stats
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


#' Calculate Probe-level p-values
#'
#' calculates p-values on the matrix (can be sequence as well), using either a z-score global, t-stat differential, or a combination of both.
#' input a matrix that has samples on the columns and sequence probes on the rows.  pData is a design matrix that describes the layout of the
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
#' @param make.plots make plots of the results
#' @param combine - what combining meta p-value method to use when combining the
#' t- and z- tests (Wilkinson's max (max) only supported for now)
#' @param debug - output debugging information
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
        make.plots = TRUE,
        combine="max",
        debug=FALSE
) {

    current_mat = probe_mat[,pData$TAG]
    use.ez = FALSE;
    if (use == "both") {
        use.t = TRUE;
        use.z = TRUE;
        use.c = TRUE;
    } else if (use == "t") {
        use.t = TRUE;
        use.z = FALSE;
        use.c = FALSE;
    } else if (use == "z") {
        use.t = FALSE;
        use.z = TRUE;
        use.c = FALSE;
    } else
    {
        stop("Unknown use paramater:" , use);
    }

    if (debug) {
        message("use:", use);
        message("combine:",combine);
        message("t.sd_shift:", t.sd_shift);
        message("t.abs_shift:", t.abs_shift);
        message("z.sd_shift:", z.sdshift);
    }
    post.cols.orig = pData$TAG[tolower(pData$visit) == "post"]
    post.cols = pData$ptid[tolower(pData$visit) == "post"]

    if (use.t) {
        message("differential t-test")
        pvaluet_df = calcProbePValuesTUnpaired(
            current_mat,
            pData,
            sd_shift = t.sd_shift,
            abs_shift = t.abs_shift
        ); #Differential t-test
        if (debug) {print(colnames(pvaluet_df));}
        if (make.plots) {plot(current_mat[,post.cols.orig[1]],pvaluet_df[,post.cols[1]], xlab="Signal", ylab="diff p-value", pch='.', log="y")}
        praw = pvaluet_df;
    }
    if (use.z) {
        if (debug) {message("global z-test");}
        pvaluez_df = calcProbePValuesZ(current_mat, pData, sd_shift = z.sdshift); #Global z-test
        parsz = attr(pvaluez_df, "pars")
        if (debug) {print(parsz);}
        if (debug) {message("global mean:",parsz$global_mean, " global sd:",parsz$global_sd);}
        if (make.plots) {plot(current_mat[,post.cols.orig[1]],pvaluez_df[,post.cols[1]], xlab="Signal", ylab="global p-value", pch='.', log="y")}
        praw = pvaluez_df;
    }

    if (use.c) {
        use.cols = intersect(colnames(pvaluet_df),colnames(pvaluez_df))
        message("combining");
        if (combine == "max") {
            praw = pvaluet_df[,use.cols];
            for (col in use.cols) {
                praw[,col] = pmax(praw[,col], pvaluez_df[,col]);
                praw[,col] = stats::pbeta(praw[,col], 2, 1);  # Probability of getting a result smaller than the max. (Wilkinsons max p-value)
            }
        } else {
            stop("Unknown combine operation")
        }
        if (make.plots){plot(pvaluet_df[,post.cols[1]],pvaluez_df[,post.cols[1]], xlab="t p-value", ylab="z p-value", pch='.')}
        if (make.plots){plot(current_mat[,post.cols.orig[1]],praw[,post.cols[1]], xlab="Signal", ylab = "Combined p-value", pch='.', log="y")}
    }
    praw[is.na(praw)] = 1.0 # Conservatively set NAs to p-value 1
    padj_df = praw;
    message("adjusting using BH");
    #print(head(padj_df));

    for (col_idx in 1:ncol(padj_df)) {
        padj_df[,col_idx] = stats::p.adjust(padj_df[,col_idx], method="BH");
    }
    if (make.plots) {plot(current_mat[,post.cols.orig[1]],padj_df[,post.cols[1]], xlab="Signal", ylab="adjusted p-value", pch='.', log="y")}

    ans = padj_df;
    attr(ans, "pvalue") = praw;

    if (debug && use.t) { attr(ans, "t") = pvaluet_df; }
    if (debug && use.z) { attr(ans, "z") = pvaluez_df; }

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

    if (all || missing(pData) || is.null(pData) || all) {
        message("No pData or all asked for, calculating Z-score on all columns");
        ans = matrix(NA, nrow = nrow(probe_mat), ncol=ncol(probe_mat));
        global_mean = mean(unlist(probe_mat), na.rm=TRUE); #Some values may have been removed (NA).
        global_sd = stats::sd(unlist(probe_mat), na.rm=TRUE); #Some values may have been removed (NA).
        zvalues = (probe_mat - global_mean) / global_sd;
        pvalues = apply((zvalues - sd_shift), 2, stats::pnorm, lower.tail=FALSE);
        pars = c(global_mean, global_sd);
        attr(pvalues, "pars") = pars;
        attr(pvalues, "zscore") = zvalues;
        return(pvalues);
    }

    pre_cols = pData$TAG[tolower(pData$visit)=="pre"]
    post_cols = pData$TAG[tolower(pData$visit) == "post"];
    post_names = pData$ptid[tolower(pData$visit) == "post"];

    ans = matrix(NA, nrow = nrow(probe_mat), ncol=length(post_cols));

    global_mean = mean(unlist(probe_mat[,c(pre_cols, post_cols)]), na.rm=TRUE);
    #global_var = var(unlist(probe_mat[,c(pre_cols, post_cols)], na.rm=TRUE);
    global_sd = stats::sd(unlist(probe_mat[,c(pre_cols, post_cols)]), na.rm=TRUE);


    pars = c(global_mean, global_sd);
    names(pars) = c("global_mean","global_sd");

    post_mat = probe_mat[,post_cols];

    post_zvalues = (post_mat - global_mean) / global_sd;

    post_pvalues = apply((post_zvalues - sd_shift), 2, stats::pnorm, lower.tail=FALSE);

    colnames(post_pvalues) = post_names;
    colnames(post_zvalues) = post_names;


    attr(post_pvalues,"pars") = pars;
    attr(post_pvalues,"zscore") = post_zvalues;

    return(post_pvalues);

}


#' calculate Probe p-values using a Differential unpaired t-test
#'
#' @param probe_mat numeric matrix or data.frame of values
#' @param pData design data.frame
#' @param sd_shift standard deviation shift to use when calculating p-values
#' Either sd_shift or abs_shift should be set
#' @param abs_shift absolute shift to use when calculating p-values
#' @param debug output debugging information
#'
#' @return matrix of p-values on the post columns defined in the pData matrix
#' @export
#'
#' @examples
calcProbePValuesTUnpaired<-function(
        probe_mat,
        pData,
        sd_shift=NA,
        abs_shift=NA,
        debug = FALSE
) {

    if (!is.na(sd_shift) && !is.na(abs_shift)) {
        stop("Either sd or abs can be set");
    }

    no_shift = is.na(sd_shift) && is.na(abs_shift);


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
        pre_mean = pre_means,
        pre_sds = pre_sds,
        post_mean = post_means,
        post_sds = post_sds,
        pre_stderr = pre_stderr,
        diff_mean = post_means - pre_means,
        dfree = rep(dfree, nrow(probe_mat)),
        p.value = rep(NA, nrow(probe_mat)),
        stringsAsFactors=FALSE
    );
    rownames(pars) = rownames(probe_mat)
    if (debug) {
        for (idx in 1:nrow(probe_mat)) {
            shift = 0;
            if (!is.na(sd_shift)) {
                shift = pre_sds[idx] * sd_shift
            } else if (!is.na(abs_shift)) {
                shift = abs_shift
            }

            pars$p.value[idx] = stats::t.test(
                x = t(probe_mat[idx, post_cols]),
                y = t(probe_mat[idx, pre_cols]),
                alternative = "greater",
                var.equal = TRUE,
                mu = shift)$p.value;
        }
    }

    post_mat = probe_mat[,post_cols];
    if (no_shift) {
        post_tvalues = (post_mat - pre_means)/pre_stderr;
    } else if (!is.na(sd_shift)) {
        post_tvalues = (post_mat - pre_means-sd_shift*pre_sds)/pre_stderr;
    } else {
        post_tvalues = (post_mat - pre_means-abs_shift)/pre_stderr;
    }


    colnames(post_tvalues) = post_names;
    pars = cbind(pars, post_tvalues)

    for (col_idx in 1:ncol(post_tvalues)) {
        ans[,col_idx] = stats::pt(q=post_tvalues[,col_idx], df=dfree, lower.tail=FALSE);
    }

    ans = as.data.frame(ans,stringsAsFactors=FALSE)
    rownames(ans) = rownames(probe_mat);
    colnames(ans) = post_names;

    attr(ans, "pars") = pars;

    return(ans);

}



