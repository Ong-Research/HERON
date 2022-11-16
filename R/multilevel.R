

#' Making Probe-level Calls
#'
#' \code{makeProbeCalls} returns call information on a probe matrix that has
#' been scored by the
#' @param probe_sample_padj a matrix of adjusted p-values where each row is a
#' feature and column is a sample
#' @param probe_cutoff cutoff to use when calling probes
#' @param pData optional design matrix for condition calling.
#' @param one_hit_filter Indicator to remove probe hits that do not have a
#' matching consecutive probe and if the probe is only found in 1 sample
#'
#' @return list of results
#' sample_probes -> logical matrix of calls on probes
#'
#' @export
#'
#' @examples
#' data(heffron2020_wuhan)
#' probe_meta <- attr(heffron2020_wuhan, "probe_meta")
#' pData <- attr(heffron2020_wuhan, "pData")
#' pval_res <- calcProbePValuesSeqMat(heffron2020_wuhan, probe_meta, pData)
#' calls_res <- makeProbeCalls(pval_res)
makeProbeCalls<-function(
    probe_sample_padj,
    pData,
    probe_cutoff = 0.05,
    one_hit_filter = TRUE
) {

    probe_calls = makeCalls(probe_sample_padj, probe_cutoff, pData);

    if (one_hit_filter) {
        ohp = oneHitProbes(probe_calls$sample);
        #Instead of removing, just set the padj values to 1.0, and remake calls.
        one_hit_padj = probe_sample_padj;
        one_hit_padj[rownames(one_hit_padj) %in% ohp,] = 1.0;
        probe_calls = makeCalls(one_hit_padj, probe_cutoff);
    }

    #Make a list of all of the results
    ans = probe_calls;
    ans$probe_sample_padj = probe_sample_padj;
    if (one_hit_filter) {
        ans$one_hit_padj = one_hit_padj;
        ans$one_hit = ohp;
        ans$one_hit_padj = one_hit_padj;
        ans$orig_padj = probe_sample_padj;
    }
    #Save parameters.
    ans$probe_cutoff = probe_cutoff;
    ans$one_hit_filter = one_hit_filter;
    return(ans);
}

#' Find one hit probes
#'
#' @param sample_probes logical probe matrix where the rows are probe
#' identifiers and the columns are samples
#'
#' @return vector of probes that are one-hits
#'
#' @export
#'
#' @examples
#' hit_mat = data.frame(
#' row.names = c("A;1","A;2","A;3","A;4"),
#' sample1 = c(TRUE, FALSE, FALSE, TRUE),
#' sample2 = c(TRUE, TRUE, FALSE, FALSE),
#' sample3 = c(TRUE, TRUE, FALSE, FALSE)
#' )
#' oneHitProbes(hit_mat)
oneHitProbes<-function(sample_probes) {
    #Probe appears in only 1 sample.
    k1_probes = rownames(sample_probes)[rowSums(sample_probes) == 1];
    probe_hit_support = probeHitSupported(sample_probes);
    supported = rowSums(sample_probes & probe_hit_support) > 0;
    #Mark any k1 probe without consecutive probe support.
    ans = intersect(k1_probes, names(supported)[!supported]);
    return(ans);
}

probeHitSupportedSingle<-function(current_df, tiling, cols) {
    current_tile = tiling[current_df$Protein[1]];
    ans = current_df;
    rownames(ans) = paste0("p",as.character(ans$Pos));

    pos_df = data.frame(
        orig = rownames(current_df),
        pos = current_df$Pos,
        posl = current_df$Pos - current_tile,
        posr = current_df$Pos + current_tile
    );

    pos_df$pos.label = paste0("p", pos_df$pos);
    pos_df$posl.label = paste0("p", pos_df$posl);
    pos_df$posr.label = paste0("p", pos_df$posr);
    rownames(pos_df) = pos_df$pos.label;


    ansl = ans;
    ansl[pos_df$pos.label,cols] =
        ans[pos_df$pos.label,cols] & ans[pos_df$posl.label,cols];
    NAs = is.na(ansl);
    ansl[NAs] = FALSE;
    ansr = ans;
    ansr[pos_df$pos.label,cols] =
        ans[pos_df$pos.label,cols] & ans[pos_df$posr.label,cols];
    NAs = is.na(ansr);
    ansr[NAs] = FALSE;
    ans.or = ans;
    ans.or[pos_df$pos.label,cols] = ansl[pos_df$pos.label,cols] |
        ansr[pos_df$pos.label,cols];
    rownames(ans.or) = pos_df[rownames(ans.or),"orig"];
    return(ans.or);
}

#' Find probe hits with a consecutive probe or another sample
#'
#' @param hit_mat matrix of logical values that indicate a hit with a
#' TRUE value
#'
#' @return matrix of logical values indicate that the TRUE hit is supported by
#' a consecutive probe hit in the sample sample or the within another sample
probeHitSupported<-function(hit_mat) {
    probes = rownames(hit_mat)
    proteins = getProteinLabel(probes)
    positions = getProteinStart(probes)
    tiling = getProteinTiling(probes)
    hit_df = as.data.frame(hit_mat, stringsAsFactors=FALSE);
    cols = seq_len(ncol(hit_df))
    hit_df$Protein = proteins;
    hit_df$Pos = positions;
    hit_df$Order = seq_len(nrow(hit_df))
    hit_df_protein = split(hit_df, hit_df$Protein)
    hit_df_protein = lapply(
        hit_df_protein,
        probeHitSupportedSingle,
        tiling = tiling,
        cols = cols
    )

    ans.dt = data.table::rbindlist(hit_df_protein);
    ans.df = as.data.frame(ans.dt, stringsAsFactors=FALSE);
    ans.df = ans.df[ans.df$Order,];
    rownames(ans.df) = paste0(ans.df$Protein,";",ans.df$Pos);
    ans.df = ans.df[,cols];
    return(ans.df)

}


#' Make calls on an input matrix of p-adjusted values
#'
#' This function takes in a matrix of p-values and using
#' a cutoff, finds the calls on the column/sample sample-level
#' and calculates the number of samples (K) and the frequency
#' of samples (F).  If the pData object is provided with a
#' condition columns defined, then a K and F is calculated for
#' each condition.
#'
#' @param padj_mat adjusted p-value matrix
#' @param padj_cutoff adjusted cutoff
#' @param pData optional pData data.frame with the condition column defined
#'
#'  If the condition column is defined, then K of N will also be reported on
#'  each condition named for the samples.
#'
#' @return list of results
#' sample => calls made on the sample (column)-level
#' k_of_n => data.frame with K of N statistics.
#' @export
#'
#' @examples
#' data(heffron2020_wuhan)
#' probe_meta <- attr(heffron2020_wuhan, "probe_meta")
#' pData <- attr(heffron2020_wuhan, "pData")
#' pval_res <- calcProbePValuesSeqMat(heffron2020_wuhan, probe_meta, pData)
#' calls_res <- makeCalls(pval_res)
makeCalls<-function(padj_mat, padj_cutoff = 0.05, pData) {
    padj_mat[is.na(padj_mat)] = 1.0; #Set all NAs to FDR=1.
    calls = padj_mat < padj_cutoff;
    minFDRs = calcMinFDR(as.matrix(padj_mat),
        additional_stats = FALSE, sort = FALSE);
    k_of_n = minFDRs
    colnames(k_of_n) = paste0("K", seq_len(ncol(minFDRs)),".padj");
    K = rowSums(minFDRs < padj_cutoff)
    F = K / ncol(minFDRs)
    K.padj = rep(1, nrow(padj_mat))
    for (idx in seq_len(nrow(padj_mat))) {
        if (K[idx] > 0) {
            K.padj[idx] = k_of_n[idx,K[idx]]
        }
    }
    k_of_n_prefix = cbind(K, F, K.padj)
    colnames(k_of_n_prefix) = c("K", "F", "K.padj")

    if (!missing(pData) && "condition" %in% colnames(pData)) {
        message("Adding in K of N for conditions");
        condition_tbl = table(pData$condition);
        for (condition in names(condition_tbl)) {
            postCols = pData$ptid[pData$visit == "post" &
                pData$condition == condition];
            K_condition = rowSums(calls[,postCols]);
            F_condition = K_condition / length(postCols);
            klabel = paste0("K.", condition);
            flabel = paste0("F.", condition);
            k_of_n_prefix = cbind(k_of_n_prefix, K_condition)
            colnames(k_of_n_prefix)[ncol(k_of_n_prefix)] = klabel;
            k_of_n_prefix = cbind(k_of_n_prefix, F_condition);
            colnames(k_of_n_prefix)[ncol(k_of_n_prefix)] = flabel;
        }
    }
    k_of_n = cbind(k_of_n_prefix, k_of_n)

    o = order(k_of_n$K, decreasing=TRUE)


    ans = list();
    ans$sample = calls[o,];
    ans$k_of_n = k_of_n[o,];

    return(ans);
}



