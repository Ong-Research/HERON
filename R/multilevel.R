



#' Making Probe-level Calls
#'
#' \code{makeProbeCalls} returns call information on a probe matrix that has
#' been scored by the calcProbePValuesSeqMat or calcProbePValuesProbeMat
#' functions.
#' @param probe_sample_padj a probe matrix of adjusted p-values
#' @param padj_cutoff cutoff to use when calling probes
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
    padj_cutoff = 0.05,
    one_hit_filter = TRUE
) {

    calls = makeCalls(probe_sample_padj, padj_cutoff, pData);

    if (one_hit_filter) {
        ohp = oneHitProbes(calls$sample);
        ## Instead of removing, set the padj values to 1.0, and remake calls.
        one_hit_padj = probe_sample_padj;
        one_hit_padj[rownames(one_hit_padj) %in% ohp,] = 1.0;
        calls = makeCalls(one_hit_padj, padj_cutoff, pData);
    }

    ## Make a list of all of the results
    ans = calls;
    ans$probe_sample_padj = probe_sample_padj;
    if (one_hit_filter) {
        ans$one_hit_padj = one_hit_padj;
        ans$one_hit = ohp;
        ans$orig_padj = probe_sample_padj;
    }
    #Save parameters.
    ans$probe_cutoff = padj_cutoff;
    ans$one_hit_filter = one_hit_filter;
    return(ans);
}

#' Find one hit probes
#'
#' @param sample_probes logical probe matrix from makeCalls
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
    hit_support = probeHitSupported(sample_probes);
    supported = rowSums(sample_probes & hit_support) > 0;
    ## Mark any k1 probe without consecutive probe support.
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
    hit_df_protein = split(hit_df, hit_df$Protein)
    hit_df_protein = lapply(
        hit_df_protein,
        probeHitSupportedSingle,
        tiling = tiling,
        cols = cols
    )
    ans_dt = data.table::rbindlist(hit_df_protein);
    ans_df = as.data.frame(ans_dt, stringsAsFactors=FALSE);
    rownames(ans_df) = paste0(ans_df$Protein,";",ans_df$Pos);
    ans_df = ans_df[probes, cols];  ##Reorder rows to the original matrix
    return(ans_df)

}

#' Find One-hit epitopes
#'
#' @param sample_epitopes logical epitope matrix from makeCalls
#'
#' @return vector of one-hit, one-probe epitopes
#' @export
#'
#' @examples
#' hit_mat = data.frame(
#' row.names = c("A_1_1","A_2_2","A_3_3","A_4_4"),
#' sample1 = c(TRUE, FALSE, FALSE, TRUE),
#' sample2 = c(TRUE, TRUE, FALSE, FALSE),
#' sample3 = c(TRUE, TRUE, FALSE, FALSE)
#' )
#' oneHitEpitopes(hit_mat)
oneHitEpitopes<-function(sample_epitopes) {
    one_probe_epitopes = oneProbeEpitopes(rownames(sample_epitopes))
    one_call_epitopes = rowSums(sample_epitopes) == 1
    one_hit_epitopes = one_probe_epitopes & one_call_epitopes
    ans = rownames(sample_epitopes)[one_hit_epitopes]

    return(ans);
}

#' Make Epitope Calls
#'
#' @param epitope_sample_padj epitope matrix of p-values
#' @param pData optional data.frame of design
#' @param padj_cutoff p-value cutoff to use
#' @param one_hit_filter filter one hit epitopes?
#'
#' @return a list of results
#' @export
#'
#' @examples
#' data(heffron2020_wuhan)
#' probe_meta <- attr(heffron2020_wuhan, "probe_meta")
#' pData <- attr(heffron2020_wuhan, "pData")
#' pr_pval_res <- calcProbePValuesSeqMat(heffron2020_wuhan, probe_meta, pData)
#' pr_calls_res <- makeProbeCalls(pr_pval_res)
#' epi_segments_uniq_res <- findEpitopeSegments(
#' probe_calls = pr_calls_res,
#' segment.method = "unique"
#' );
#' epi_segments_uniq_res <- findEpitopeSegments(
#' probe_calls = pr_calls_res,
#' segment.method = "unique"
#' );
#' epi_pval_uniq <- calcEpitopePValuesMat(
#' probe_pvalues = attr(pr_pval_res, "pvalue"),
#' epitope_ids = epi_segments_uniq_res,
#' metap_method = "wilkinsons_max1"
#' )
#' epi_padj_uniq <- p_adjust_mat(epi_pval_uniq, method="BH")
#' makeEpitopeCalls(epi_padj_uniq)
makeEpitopeCalls<-function(
    epitope_sample_padj,
    pData,
    padj_cutoff = 0.05,
    one_hit_filter = TRUE) {

    calls = makeCalls(epitope_sample_padj, padj_cutoff, pData);

    if (one_hit_filter) {
        ohe = oneHitEpitopes(calls$sample);
        ## Instead of removing, set the padj values to 1.0, and remake calls.
        one_hit_padj = epitope_sample_padj;
        one_hit_padj[rownames(one_hit_padj) %in% ohe,] = 1.0;
        calls = makeCalls(one_hit_padj, padj_cutoff, pData);
    }

    ## Make a list of all of the results
    ans = calls;
    ans$epitope_sample_padj = epitope_sample_padj;
    if (one_hit_filter) {
        ans$one_hit_padj = one_hit_padj;
        ans$one_hit = ohe;
        ans$orig_padj = epitope_sample_padj;
    }
    ## Save parameters.
    ans$epitope_cutoff = padj_cutoff;
    ans$one_hit_filter = one_hit_filter;
    return(ans);

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
    Fr = K / ncol(minFDRs)
    K.padj = rep(1, nrow(padj_mat))
    for (idx in seq_len(nrow(padj_mat))) {
        if (K[idx] > 0) {
            K.padj[idx] = k_of_n[idx,K[idx]]
        }
    }
    k_of_n_prefix = cbind(K, Fr, K.padj)
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



