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
#' @param p.adjust.method method for adjusting p-values (multiple testing)
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
        p.adjust.method = "BH",
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
    message("adjusting using ", p.adjust.method);
    #print(head(padj_df));

    for (col_idx in 1:ncol(padj_df)) {
        padj_df[,col_idx] = stats::p.adjust(padj_df[,col_idx], method=p.adjust.method);
    }
    if (make.plots) {plot(current_mat[,post.cols.orig[1]],padj_df[,post.cols[1]], xlab="Signal", ylab="adjusted p-value", pch='.', log="y")}

    ans = padj_df;
    attr(ans, "pvalue") = praw;

    if (debug && use.t) { attr(ans, "t") = pvaluet_df; }
    if (debug && use.z) { attr(ans, "z") = pvaluez_df; }

    return(ans);
}


#' Title
#'
#' @param seq_mat
#' @param probe_meta
#' @param pData
#' @param t.sd_shift
#' @param t.abs_shift
#' @param z.sdshift
#' @param use
#' @param make.plots
#' @param combine
#'
#' @return
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
        make.plots = TRUE,
        combine = "max" #Maybe eliminate this parameter, or allow user to change the meta p-value method?
) {

    seq_results = calcProbePValuesProbeMat(
        probe_mat = seq_mat,
        pData = pData,
        t.sd_shift = t.sd_shift,
        t.abs_shift = t.abs_shift,
        z.sdshift = z.sdshift,
        use = use,
        make.plots = make.plots,
        combine = combine
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
        message("No pData or all asked for, calculating Z-score on all columns");
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

#' Calculate Epitope-level p-values
#'
#' @param probe_pvalues matrix of probe-level p-values
#' @param epitope_ids vector of epitope ids to calculate p-values from
#' (needs to be a unique set)
#' @param method meta p-value method to use (see calcProteinValues)
#' @param do.sort sort the p-value after calculation?
#' @param data_matrix matrix of values to calculate the co-variances for some
#' meta p-value methods (ebrown and kosts)
#' @param probes speed up parameter
#'
#' @return matrix of epitope p-values where each row is an epitope id
#' @export
#'
#' @examples
calcEpitopePValues<-function(
    probe_pvalues,
    epitope_ids,
    method = "min",
    do.sort=FALSE,
    data_matrix = NULL,
    probes = rownames(probe_pvalues)
    ) {

    if (length(epitope_ids) == 0) {return(data.frame());}

    if (method == "maxFDR" || method == "minFDR") {

    }



    #Call each epitope a "protein" and call the calcProteinPValues code.
    meta = NULL;
    meta_list = list();
    for (epitope_idx in 1:length(epitope_ids)) {
        eprobes = getEpitopeProbeIDs(epitope_ids[epitope_idx]);
        meta_list[[epitope_idx]] =
            data.frame(
                Epitope_ID = rep(epitope_ids[epitope_idx], length(eprobes)),
                PROBE_ID = eprobes,
                stringsAsFactors=FALSE
            );
    }
    meta.dt = data.table::rbindlist(meta_list);
    meta = as.data.frame(meta.dt, stringsAsFactors=FALSE);


    meta = meta[meta$PROBE_ID %in% probes,]; #Need to do this in case of different tilings.

    eprobe_pvalues = probe_pvalues[meta$PROBE_ID,]

    #Make sure pvalues are valid.

    if (sum(is.na(eprobe_pvalues))) {
        message("NAs detected in pvalues... Correcting to 1")
        eprobe_pvalues[is.na(eprobe_pvalues)] = 1;
    }


    if (sum(eprobe_pvalues>1) > 0) {
        message("Some p-values are > 1... Correcting to 1.")
        eprobe_pvalues[eprobe_pvalues > 1] = 1;
    }
    if (sum(eprobe_pvalues<0) > 0) {
        message("Some pvalues are < 0... Correcting to 0.")
        eprobe_pvalues[eprobe_pvalues < 0] = 0;
    }

    edata_matrix = NULL;
    if (!is.null(data_matrix)) {
        edata_matrix = data_matrix[meta$PROBE_ID,];
    }
    ans =
        calcProteinPValuesMat(
            pvalues_mat = eprobe_pvalues,
            method = method,
            data_matrix = edata_matrix,
            probes = meta$PROBE_ID,
            proteins = meta$Epitope_ID #Use the epitope_id as the "protein"
        )

    attr(ans, "emeta") = meta;


    return(ans);


}


#' Calculate protein-level p-values from a matrix of p-values
#'
#' @param pvalues_mat matrix of probe-level p-values
#' @param method meta p-value method to use
#' @param data_matrix matrix of values for Kost's and EBrown's method
#' @param probes probes
#' @param proteins proteins
#'
#' @return matrix of protein-level p-values
#' @export
#'
#' @examples
calcProteinPValuesMat<-function(
        pvalues_mat,
        method="min",
        data_matrix=NULL,
        probes=rownames(pvalues_mat),
        proteins=getProteinLabel(probes)
) {

    if (sum(is.na(pvalues_mat))) {
        message("NAs detected in pvalues... Correcting to 1")
        pvalues_mat[is.na(pvalues_mat)] = 1;
    }

    if (sum(pvalues_mat>1) > 0) {
        message("Some p-values are > 1... Correcting to 1.")
        pvalues_mat[pvalues_mat > 1] = 1;
    }
    if (sum(pvalues_mat<0) > 0) {
        message("Some pvalues are < 0... Correcting to 0.")
        pvalues_mat[pvalues_mat < 0] = 0;
    }


    protein_pvalues = NULL;
    for (col_idx in 1:ncol(pvalues_mat)) {
        current = calcProteinPValuesVec(
            proteins = proteins,
            pvalues = pvalues_mat[,col_idx],
            method = method,
            do.sort = FALSE,
            data_matrix = data_matrix,
            probes=probes
        );
        protein_pvalues = cbind(protein_pvalues, current$Protein.pvalue);
        proteins_row = current$Protein;
    }

    NAs = is.na(protein_pvalues);
    protein_pvalues[NAs] = 1.0;


    rownames(protein_pvalues) = proteins_row;
    colnames(protein_pvalues) = colnames(pvalues_mat);

    return(protein_pvalues);
}


#' Calculate protein p-values given a vector of probe-level p-values
#'
#' @param pvalues named vector of probes/eptiopes
#' @param method meta p-value method to use
#' @param do.sort sort in increasing order?
#' @param data_matrix matrix of values for Kost's and EBROWN's method
#' @param probes probe ids on vector
#' @param proteins protein labels on vector
#'
#' @return vector of protein-level p-values
#' @export
#'
#' @examples
calcProteinPValuesVec<-function(pvalues, method="min_bonf", do.sort = FALSE,
                                data_matrix = NULL, probes = names(pvalues),
                                proteins = getProteinLabel(probes)) {

    if (method == "minFDR") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins), function(l) {return(min(l,na.rm=TRUE))});
    }
    else if (method == "maxFDR") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins), function(l) {return(max(l,na.rm=TRUE))});
    }
    else if (method == "min") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins), function(l) {return(min(l,na.rm=TRUE))});
    } else if (method == "min_bonf") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            return(stats::p.adjust(min(l,na.rm=TRUE),"bonf",length(l)))
                        }
        );
    } else if (method == "fischer" || method == "sumlog") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l[l <= 0] = .Machine$double.xmin;
                            l[l > 1] = 1; # Don't think this ever will happen..

                            return(metap::sumlog(l)$p)
                        }
        );

    } else if (method == "stouffer" || method == "sumz") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l[l <= 0] = .Machine$double.xmin;
                            l[l > 1] = 1; # Don't think this ever will happen, but just in case...
                            return(metap::sumz(l)$p)
                        }
        );

    } else if (method == "meanz") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representive postive number
                            l[l <= 0] = .Machine$double.xmin;
                            l[l > 1] = 1; # Don't think this ever will happen, but just in case...
                            return(metap::meanz(l)$p)
                        }
        );

    } else if (method == "lancaster") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::invchisq(l,length(l))$p);
                        }
        );

    } else if (method == "invt") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::invt(l,length(l))$p)
                        }
        );

    } else if (method == "logitp") {

        ans = stats::aggregate(
            pvalues,
            by=list(Protein = proteins),
            function(l) {
                if (length(l) == 0) {return(1.0);}
                if (length(l) == 1) {return(l[1]);}
                #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                l = min_max(l, .Machine$double.xmin, 1);
                return(metap::logitp(l)$p)
            }
        );

    } else if (method == "meanp") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}

                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::meanp(l)$p)
                        }
        );


    } else if (method == "edgington" || method == "sump") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::sump(l)$p)
                        }
        );

    } else if (method == "hmp" || method == "harmonicmeanp") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
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
    } else if (method == "wilkinsons_min1" || method == "tippets") { #Wilkinson's with r=1 is tippet's method
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::minimump(l)$p)
                        }
        );


    } else if (method == "wilkinsons_min2" || method == "min2") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::wilkinsonp(l,r=2)$p)
                        }
        );

    } else if (method == "wilkinsons_min3") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            if (length(l) == 2) {return(metap::wilkinsonp(l,r=2)$p)}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::wilkinsonp(l,r=3)$p)
                        }
        );

    } else if (method == "wilkinsons_min4") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            if (length(l) <= 4) {return(metap::wilkinsonp(l,r=length(l))$p)}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::wilkinsonp(l,r=4)$p)
                        }
        );
    } else if (method == "wilkinsons_min5") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            if (length(l) <= 5) {return(metap::wilkinsonp(l,r=length(l))$p)}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::wilkinsonp(l,r=5)$p)
                        }
        );
    } else if (method == "wilkinsons_max1" || method=="max") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {
                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(l[1]);}
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::maximump(l)$p)
                        }
        );

    } else if (method == "wilkinsons_max2" || method=="max2") {
        ans = stats::aggregate(pvalues, by=list(Protein = proteins),
                        function(l) {

                            if (length(l) == 0) {return(1.0);}
                            if (length(l) == 1) {return(1.0);} #Conservative, assume if there was another p-value, it would be 1
                            r = length(l) - 1;
                            #Make sure p-values are well-behaved, p-values of value 0 => smallest representative positive number
                            l = min_max(l, .Machine$double.xmin, 1);
                            return(metap::wilkinsonp(l, r=r)$p)
                        }
        );

    } else if (method == "kosts") {
        #cat("Sort\n");
        data_matrix = as.data.frame(data_matrix[probes,]);
        nc = 1:ncol(data_matrix);
        data_matrix$pvalue = pvalues;
        #cat("Split\n");
        data_list = split(data_matrix, proteins);

        ans_list = lapply(data_list, function(l) {
            if (nrow(l) == 0) {
                return(1);
            }
            if (nrow(l) == 1) {
                return(l$pvalue[1]);
            }
            kost.res = kostsMethodFast(l[,nc], l$pvalue);
            return(kost.res);
        }
        )

        ans = data.frame(
            Protein = names(data_list),
            pvalue = unlist(ans_list),
            stringsAsFactors=FALSE
        );
    } else if (method == "ebrown") {
        #cat("Sort\n");
        data_matrix = as.data.frame(data_matrix[probes,]);
        nc = 1:ncol(data_matrix);
        data_matrix$pvalue = pvalues;
        #cat("Split\n");
        data_list = split(data_matrix, proteins);

        ans_list = lapply(data_list, function(l) {
            if (nrow(l) == 0) {
                return(1);
            }
            if (nrow(l) == 1) {
                return(l$pvalue[1]);
            }
            brown.res = EmpiricalBrownsMethod::empiricalBrownsMethod(l[,nc], l$pvalue);
            return(brown.res);
        }
        )

        ans = data.frame(
            Protein = names(data_list),
            pvalue = unlist(ans_list),
            stringsAsFactors=FALSE
        );
    } else if (method == "cct") {
        ans = stats::aggregate(
            pvalues,
            by=list(
                Protein = proteins
            ),
            function(l) {
                if (length(l) == 0) {
                    return(1);
                }
                if (length(l) == 1) {
                    return(l[1]);
                }
                #Prevent a return of 1 when there is one p-value that is 1?  Read the reference!
                #Make sure p-values are well-behaved, p-values of value 0 => smallest representive postive number
                l[l <= 0] = .Machine$double.xmin;
                max.p = 1-1e-10;
                l[l > max.p] = max.p


                return(CCT(pvals = l))
            }
        );
#    } else if (method == "BJ") {
#        ans = callGBJ(pvalues, probes, proteins, data_matrix, method="BJ");
#    } else if (method == "GBJ") {
#        ans = callGBJ(pvalues, probes, proteins, data_matrix, method="GBJ");
#    } else if (method == "HC") {
#        ans = callGBJ(pvalues, probes, proteins, data_matrix, method="HC");
#    } else if (method == "minP") {
#        ans = callGBJ(pvalues, probes, proteins, data_matrix, method="minP");
#    } else if (method == "BJI") {
#        ans = callGBJ(pvalues, probes, proteins, data_matrix, method="BJ", indep=TRUE);
#    } else if (method == "GBJI") {
#        ans = callGBJ(pvalues, probes, proteins, data_matrix, method="GBJ", indep=TRUE);
#    } else if (method == "HCI") {
#        ans = callGBJ(pvalues, probes, proteins, data_matrix, method="HC", indep=TRUE);
#    } else if (method == "minPI") {
#        ans = callGBJ(pvalues, probes, proteins, data_matrix, method="minP", indep=TRUE);
    } else {
        stop("Unknown method ",method);
    }


    ansn = stats::aggregate(pvalues, by=list(Protein = proteins), function(l){return(length(l))});


    colnames(ansn)[2] = "NProbes";

    ans = cbind(ansn, ans[,2]);

    colnames(ans)[ncol(ans)] = "Protein.pvalue";
    ans$Protein.padj = stats::p.adjust(ans$Protein.pvalue);
    if (do.sort) {
        ans = ans[order(ans$Protein.pvalue, decreasing=FALSE),]
    }
    return(ans);



}



