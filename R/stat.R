#' Calculate minimum FDRs across samples
#'
#' This code calculates the minimum FDR threshold needed
#' to call K out of N samples, where samples are in the columns
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

    attr(ans, "freq") = rowSums(ans) / ncol(ans);

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

    if (all || missing(pData) || is.null(pData) || all) {
        message("No pData or all asked for, calculating on all columns");
        global_mean = mean(unlist(probe_mat), na.rm=TRUE);
        global_sd = stats::sd(unlist(probe_mat), na.rm=TRUE);
        zvalues = (probe_mat - global_mean) / global_sd;
        pvalues = apply((zvalues - sd_shift), 2, stats::pnorm, lower.tail=FALSE)
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
    global_sd = stats::sd(
        unlist(probe_mat[,c(pre_cols, post_cols)]), na.rm=TRUE
    );


    pars = c(global_mean, global_sd);
    names(pars) = c("global_mean","global_sd");

    post_mat = probe_mat[,post_cols];

    post_zvalues = (post_mat - global_mean) / global_sd;

    post_pvalues = apply(
        (post_zvalues - sd_shift), 2, stats::pnorm, lower.tail=FALSE
    );

    colnames(post_pvalues) = post_names;
    colnames(post_zvalues) = post_names;


    attr(post_pvalues,"pars") = pars;
    attr(post_pvalues,"zscore") = post_zvalues;

    return(post_pvalues);

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
    pre_df = pData[tolower(pData$visit) =="pre",]
    post_df = pData[tolower(pData$visit) == "post",];
    post_names = pData$ptid[tolower(pData$visit) == "post"];
    if (!is.na(sd_shift) && !is.na(abs_shift)) {
        stop("Either sd or abs can be set. Not both.");
    }
    no_shift = is.na(sd_shift) & is.na(abs_shift);
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
    #print(mapping)
    nx = nrow(mapping)
    dfree = nx-1;

    ans = matrix(NA, nrow = nrow(probe_mat), ncol=nrow(mapping))
    colnames(ans) = mapping$ptid;
    rownames(ans) = rownames(probe_mat);

    diff_mat = probe_mat[,mapping$post] - probe_mat[,mapping$pre];
    colnames(diff_mat) = mapping$ptid;
    rownames(diff_mat) = rownames(probe_mat);

    pars = data.frame(
        diff_mean = rep(NA, nrow(probe_mat)),
        diff_sd = rep(NA, nrow(probe_mat)),
        diff_var = rep(NA, nrow(probe_mat)),
        diff_stderr =rep(NA, nrow(probe_mat)),
        dfree = rep(dfree, nrow(probe_mat)),
        pvalue = rep(NA, nrow(probe_mat)),
        stringsAsFactors=FALSE
    );
    rownames(pars) = rownames(probe_mat);

    for (row_idx in seq_len(nrow(probe_mat))) {
        x = t(probe_mat[row_idx,mapping$post] - probe_mat[row_idx,mapping$pre]);
        #print(x)
        nx = length(x)
        mx = mean(x, na.rm=TRUE);
        sx = stats::sd(x, na.rm=TRUE)
        vx = stats::var(x, na.rm=TRUE);
        dfree = nx - 1
        stderr = sqrt(vx/nx)
        current_df = data.frame(
            Pre = c(t(probe_mat[row_idx, mapping$pre])),
            Post = c(t(probe_mat[row_idx, mapping$post]))
        )
        rownames(current_df) = post_names
        if (!is.na(abs_shift)) {
            current_df$Post = current_df$Post - abs_shift
        } else if (!is.na(sd_shift)) {
            pre_sd = stats::sd(current_df$Pre, na.rm = TRUE);
            current_df$Post = current_df$Post - sx * sd_shift;
        }
        current_df = stats::na.omit(current_df); #Keep only complete pairs

        t.test.res = stats::t.test(
            current_df$Post, current_df$Pre,
            paired=TRUE, alternative="greater"
        );

        pars$diff_mean[row_idx] = mx;
        pars$diff_var[row_idx] = vx;
        pars$diff_sd[row_idx] = sx;
        pars$diff_stderr[row_idx] = stderr;
        pars$pvalue[row_idx] = t.test.res$p.value;
        pars$dfree[row_idx] = nrow(current_df) - 1;

        for (col_idx in seq_len((nrow(mapping)))) {
            if (no_shift) {
                tstat = (x[col_idx])/stderr;
            } else if (!is.na(sd_shift)) {
                tstat = (x[col_idx]-sd_shift*sx)/stderr;
            } else {
                tstat = (x[col_idx]-abs_shift)/stderr;
            }
            ans[row_idx, col_idx] =
                stats::pt(tstat, df = pars$dfree[row_idx], lower.tail=FALSE);
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


#' Calculate Probe p-values using a differential unpaired t-test
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
        stop("Either sd or abs can be set, not both.");
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
        for (idx in seq_len(nrow(probe_mat))) {
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

    for (col_idx in seq_len(ncol(post_tvalues))) {
        ans[,col_idx] = stats::pt(
            q=post_tvalues[,col_idx], df=dfree, lower.tail=FALSE
        );
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


#' Calculate meta p-values on a matrix (column-by-column)
#'
#' @param pvalues_mat matrix of p-values where each row is a feature and
#' each column is a sample
#' @param by_list list of grouping elements (see aggregate)
#' @param method what meta p-value method to use.
#'
#' @return matrix of meta p-values
#' @export
#'
#' @examples
calcMetaPValuesMat<-function(
        pvalues_mat,
        by_list,
        method="min"
) {

    if (sum(is.na(pvalues_mat))) {
        message("NAs detected in pvalues... Correcting to 1")
        pvalues_mat[is.na(pvalues_mat)] = 1;
    }

    if (sum(pvalues_mat>1) > 0) {
        warning("Some p-values are > 1... Correcting to 1.")
        pvalues_mat[pvalues_mat > 1] = 1;
    }
    if (sum(pvalues_mat<0) > 0) {
        warning("Some pvalues are < 0... Correcting to 0.")
        pvalues_mat[pvalues_mat < 0] = 0;
    }


    meta_pvalues = NULL;
    for (col_idx in seq_len(ncol(pvalues_mat))) {
        current = calcMetaPValuesVec(
            pvalues = pvalues_mat[,col_idx],
            method = method,
            by_list = by_list,
            do.sort = FALSE
        );
        meta_pvalues = cbind(meta_pvalues, current$Meta.pvalue);
        meta_row = current[,1];
    }

    NAs = is.na(meta_pvalues); #All NAs have p-value = 1.
    meta_pvalues[NAs] = 1.0;


    rownames(meta_pvalues) = meta_row;
    colnames(meta_pvalues) = colnames(pvalues_mat);

    return(meta_pvalues);
}



#' Calculate meta p-values given a vector of lower-level p-values
#'
#' @param pvalues vector of p-values
#' @param by_list list of groupings for meta-pvalue calculation
#' @param method meta p-value method to use
#' @param do.sort sort in increasing order?
#'
#' @return vector of meta p-values
#' @export
#'
#' @examples
calcMetaPValuesVec<-function(
        pvalues,
        by_list,
        method="min_bonf",
        do.sort = FALSE) {

    if (method == "minFDR") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                l = l[!is.na(l)]
                if (length(l) == 0) {return(1.0);}
                return(min(l))
            }
        );
    }
    else if (method == "maxFDR") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                l = l[!is.na(l)]
                if (length(l) == 0) {return(1.0);}
                return(max(l));
            }
        )
    }
    else if (method == "min") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                l = l[!is.na(l)];
                if (length(l) == 0) {return(1.0);}
                return(min(l));
            }
        );
    } else if (method == "min_bonf") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                l = l[!is.na(l)]
                if (length(l) == 0) { return(1.0);}
                return(stats::p.adjust(min(l,na.rm=TRUE),"bonf",length(l)))
            }
        );
    } else if (method == "fischer" || method == "sumlog") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                l = l[!is.na(l)]
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved for function call
                l = min_max(l, .Machine$double.xmin, 1);
                return(metap::sumlog(l)$p)
            }
        );
    } else if (method == "stouffer" || method == "sumz") {
        ans = stats::aggregate(
            pvalues,
            by = by_list,
            function(l) {
                l = l[!is.na(l)]
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved for function call
                l = min_max(l, .Machine$double.xmin, 1);
                return(metap::sumz(l)$p)
            }
        );
    } else if (method == "meanz") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                l = l[!is.na(l)]
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                l = min_max(l, .Machine$double.xmin, 1);
                return(metap::meanz(l)$p)
            }
        );
    } else if (method == "lancaster") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                l = l[!is.na(l)]
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved for function call
                l = min_max(l, .Machine$double.xmin, 1);
                return(metap::invchisq(l,length(l))$p);
            }
        );

    } else if (method == "invt") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                l = l[!is.na(l)]
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved for function call
                l = min_max(l, .Machine$double.xmin, 1);
                return(metap::invt(l,length(l))$p)
            }
        );

    } else if (method == "logitp") {

        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved
                l = min_max(l, .Machine$double.xmin, 1);
                return(metap::logitp(l)$p)
            }
        );

    } else if (method == "meanp") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved
                l = min_max(l, .Machine$double.xmin, 1);
                return(metap::meanp(l)$p)
            }
        );
    } else if (method == "edgington" || method == "sump") {
        ans = stats::aggregate(
            pvalues,
            by_list,
            function(l) {
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved
                l = min_max(l, .Machine$double.xmin, 1);
                return(metap::sump(l)$p)
            }
        );

    } else if (method == "hmp" || method == "harmonicmeanp") {
        ans = stats::aggregate(
            pvalues, by=by_list,
            function(l) {
                if (length(l) == 0) {
                    return(1);
                }
                if (length(l) == 1) {
                    return(l[1]);
                }
                return(harmonicmeanp::p.hmp(p = l, L=length(l)))
                        }
        );
    } else if (method == "wilkinsons_min1" || method == "tippets") {
        #Wilkinson's with r=1 is tippet's method
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved
                l = min_max(l, .Machine$double.xmin, 1);
                return(metap::minimump(l)$p)
            }
        );


    } else if (method == "wilkinsons_min2" || method == "min2") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved
                l = min_max(l, .Machine$double.xmin, 1);
                return(metap::wilkinsonp(l,r=2)$p)
            }
        );

    } else if (method == "wilkinsons_min3") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved
                l = min_max(l, .Machine$double.xmin, 1);
                if (length(l) < 3) {return(metap::wilkinsonp(l,r=2)$p)}
                return(metap::wilkinsonp(l,r=3)$p)
            }
        );
    } else if (method == "wilkinsons_min4") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved
                l = min_max(l, .Machine$double.xmin, 1);
                if (length(l) < 4) {return(metap::wilkinsonp(l,r=length(l))$p)}
                return(metap::wilkinsonp(l,r=4)$p)
            }
        );
    } else if (method == "wilkinsons_min5") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved
                l = min_max(l, .Machine$double.xmin, 1);

                if (length(l) < 5) {return(metap::wilkinsonp(l,r=length(l))$p)}
                            return(metap::wilkinsonp(l,r=5)$p)
                        }
        );
    } else if (method == "wilkinsons_max1" || method=="max") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            #Make sure p-values are well-behaved
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::maximump(l)$p)
                        }
        );

    } else if (method == "wilkinsons_max2" || method=="max2") {
        ans = stats::aggregate(pvalues, by=by_list,
                        function(l) {

                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(1.0);} #Conservative
                            r = length(l) - 1;
                            #Make sure p-values are well-behaved
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::wilkinsonp(l, r=r)$p)
                        }
        );
    } else if (method == "cct") {
        ans = stats::aggregate(
            pvalues,
            by=by_list,
            function(l) {
                if (length(l) == 0) {return(1);}
                if (length(l) == 1) {return(l[1]);}
                #Prevent a return of 1 when there is one p-value that is 1
                l = min_max(l, .Machine$double.xmin, 1 - 1e-10)
                return(CCT(pvals = l))
            }
        );
    } else {
        stop("Unknown method ",method);
    }


    ansn = stats::aggregate(
        pvalues,
        by = by_list,
        function(l){
            l = l[!is.na(l)]
            if (length(l) == 0){return(1);} # We return a p-value of 1
            return(length(l))
            }
        );


    colnames(ansn)[2] = "NElements";

    ans = cbind(ansn, ans[,2]);

    colnames(ans)[ncol(ans)] = "Meta.pvalue";
    ans$Meta.padj = stats::p.adjust(ans$Meta.pvalue);
    if (do.sort) {
        ans = ans[order(ans$Meta.pvalue, decreasing=FALSE),]
    }
    return(ans);
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


