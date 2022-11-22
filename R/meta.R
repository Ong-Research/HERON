
#' Calculate meta p-values on a matrix (column-by-column)
#'
#' @param pvalues_mat matrix of p-values where each row is a feature and
#' each column is a sample
#' @param by_list list of grouping elements (see aggregate)
#' @param method what meta p-value method to use.
#'
#' @return matrix of meta p-values
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


minBonfMeta<-function(pvals) {
    pvals = pvals[!is.na(pvals)]
    if (length(pvals) == 0) {return(1.0)}
    if (length(pvals) == 1) {return(pvals)}
    ans = min(pvals,na.rm=TRUE)
    ans = stats::p.adjust(ans,"bonf",length(pvals))
    return(ans);
}

minMeta<-function(pvals) {
    pvals = pvals[!is.na(pvals)]
    if (length(pvals) == 0) {return(1.0)}
    if (length(pvals) == 1) {return(pvals)}
    ans = min(pvals);
    return(ans);
}

maxMeta<-function(pvals) {
    pvals = pvals[!is.na(pvals)]
    if (length(pvals) == 0) {return(1.0)}
    if (length(pvals) == 1) {return(pvals)}
    ans = max(pvals)
    return(ans);
}

fischerMeta<-function(pvals) {
    pvals = pvals[!is.na(pvals)]
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved for function call
    pvals = min_max(pvals, .Machine$double.xmin, 1);
    return(metap::sumlog(pvals)$p)
}

stoufferMeta<-function(pvals) {
    pvals = pvals[!is.na(pvals)]
    if (length(l) == 0) {return(1.0);}
    if (length(l) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved for function call
    l = min_max(l, .Machine$double.xmin, 1);
    return(metap::sumz(l)$p)
}

meanzMeta<-function(pvals) {
    pvals = pvals[!is.na(pvals)]
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    pvals = min_max(pvals, .Machine$double.xmin, 1);
    return(metap::meanz(pvals)$p)
}

lancasterMeta<-function(pvals) {
    pvals = pvals[!is.na(pvals)]
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved for function call
    l = min_max(pvals, .Machine$double.xmin, 1);
    return(metap::invchisq(pvals, length(pvals))$p);
}

invtMeta<-function(pvals) {
    pvals = pvals[!is.na(pvals)]
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved for function call
    pvals = min_max(pvals, .Machine$double.xmin, 1);
    return(metap::invt(pvals, length(pvals))$p)
}

logitpMeta<-function(pvals) {
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved
    l = min_max(pvals, .Machine$double.xmin, 1);
    return(metap::logitp(pvals)$p)
}

meanpMeta<-function(pvals) {
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved
    pvals = min_max(pvals, .Machine$double.xmin, 1);
    return(metap::meanp(pvals)$p)
}

sumpMeta<-function(pvals) {
    if (length(l) == 0) {return(1.0);}
    if (length(l) == 1) {return(l[1]);}
    #Make sure p-values are well-behaved
    l = min_max(l, .Machine$double.xmin, 1);
    return(metap::sump(l)$p)
}

hmpMeta<-function(pvals) {
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    return(harmonicmeanp::p.hmp(p = pvals, L=length(pvals)))
}

wilkMin1Meta<-function(pvals) {
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved
    pvals = min_max(pvals, .Machine$double.xmin, 1);
    return(metap::minimump(pvals)$p)
}

wilkMin2Meta<-function(pvals) {
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved
    pvals = min_max(pvals, .Machine$double.xmin, 1);
    return(metap::wilkinsonp(pvals,r=2)$p)
}

wilkMin3Meta<-function(pvals) {
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved
    pvals = min_max(pvals, .Machine$double.xmin, 1);
    if (length(pvals) < 3) {return(metap::wilkinsonp(pvals, r=length(pvals))$p)}
    return(metap::wilkinsonp(pvals, r=3)$p)
}

wilkMin4Meta<-function(pvals) {
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved
    l = min_max(l, .Machine$double.xmin, 1);
    if (length(l) < 4) {return(metap::wilkinsonp(pvals, r=length(pvals))$p)}
    return(metap::wilkinsonp(l,r=4)$p)
}

wilkMin5Meta<-function(pvals) {
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved
    l = min_max(pvals, .Machine$double.xmin, 1);

    if (length(pvals) < 5) {return(metap::wilkinsonp(pvals, r=length(pvals))$p)}
    return(metap::wilkinsonp(pvals, r=5)$p)
}

wilkMax1Meta<-function(pvals) {
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Make sure p-values are well-behaved
    pvals = min_max(pvals, .Machine$double.xmin, 1);
    return(metap::maximump(pvals)$p)
}

wilkMax2Meta<-function(pvals) {
    if (length(pvals) == 0) {return(1.0);}
    if (length(pvals) == 1) {return(1.0);} #Conservative
    r = length(pvals) - 1;
    #Make sure p-values are well-behaved
    pvals = min_max(pvals, .Machine$double.xmin, 1);
    return(metap::wilkinsonp(pvals, r=r)$p)
}

cctMeta<-function(pvals) {
    if (length(pvals) == 0) {return(1);}
    if (length(pvals) == 1) {return(pvals[1]);}
    #Prevent a return of 1 when there is one p-value that is 1
    pvals = min_max(pvals, .Machine$double.xmin, 1 - 1e-16)
    return(CCT(pvals));
}

getMetaPFxn<-function(method="min_bonf") {
    ans = switch(
        method,
        "min_bonf" = minBonfMeta,
        "min" = minMeta,
        "minFDR" = minMeta,
        "max" = maxMeta,
        "maxFDR" = maxMeta,
        "fischer" = fischerMeta,
        "sumlog" = fischerMeta,
        "stouffer" = stoufferMeta,
        "sumz" = stoufferMeta,
        "meanz" = meanzMeta,
        "lancaster" = lancasterMeta,
        "invt" = invtMeta,
        "logitp" = logitpMeta,
        "meanp" = meanpMeta,
        "edgington" = sumpMeta,
        "sump" = sumpMeta,
        "hmp" = hmpMeta,
        "harmonicmeanp" = hmpMeta,
        "wilkinsons_min1" = wilkMin1Meta,
        "wmin1" = wilkMin1Meta,
        "tippets" = wilkMin1Meta,
        "wilkinsons_min2" = wilkMin2Meta,
        "wmin2" = wilkMin2Meta,
        "wilkinsons_min3" = wilkMin3Meta,
        "wilkinsons_min4" = wilkMin4Meta,
        "wilkinsons_min5" = wilkMin5Meta,
        "wilkinsons_max1" = wilkMax1Meta,
        "wmax1" = wilkMax1Meta,
        "wilkinsons_max2" = wilkMax2Meta,
        "wmax2" = wilkMax2Meta,
        "cct" = cctMeta
    )
    return(ans);
}


#' Calculate meta p-values given a vector of lower-level p-values
#'
#' @param pvalues vector of p-values
#' @param by_list list of groupings for meta-pvalue calculation
#' @param method meta p-value method to use
#' @param do.sort sort in increasing order?
#'
#' @return vector of meta p-values
calcMetaPValuesVec<-function(
        pvalues,
        by_list,
        method="min_bonf",
        do.sort = FALSE) {

    meta_fxn = getMetaPFxn(method);
    if (!is.null(meta_fxn)) {
        ans = stats::aggregate(
            pvalues,
            by = by_list,
            FUN = meta_fxn
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
