#' Find one hit probes
#'
#' @param sample_probes logical probe matrix from makeCalls
#'
#' @return vector of probes that are one-hits
#'
#' @export
#'
#' @examples
#' hit_mat <- data.frame(
#' row.names = c("A;1","A;2","A;3","A;4"),
#' sample1 = c(TRUE, FALSE, FALSE, TRUE),
#' sample2 = c(TRUE, TRUE, FALSE, FALSE),
#' sample3 = c(TRUE, TRUE, FALSE, FALSE)
#' )
#' oneHitProbes(hit_mat)
oneHitProbes<-function(sample_probes) {
    #Probe appears in only 1 sample.
    k1_probes <- rownames(sample_probes)[rowSums(sample_probes) == 1]
    hit_support <- probeHitSupported(sample_probes)
    supported <- rowSums(sample_probes & hit_support) > 0
    ## Mark any k1 probe without consecutive probe support.
    ans <- intersect(k1_probes, names(supported)[!supported])
    return(ans)
}

probeHitSupportedSingle<-function(current_df, tiling, cols) {
    current_tile <- tiling[current_df$Protein[1]]
    ans <- current_df
    rownames(ans) <- paste0("p",as.character(ans$Pos))

    pos_df <- data.frame(
        orig = rownames(current_df),
        pos = current_df$Pos,
        posl = current_df$Pos - current_tile,
        posr = current_df$Pos + current_tile
    )

    pos_df$pos.label <- paste0("p", pos_df$pos)
    pos_df$posl.label <- paste0("p", pos_df$posl)
    pos_df$posr.label <- paste0("p", pos_df$posr)
    rownames(pos_df) <- pos_df$pos.label


    ansl <- ans
    ansl[pos_df$pos.label,cols] <-
        ans[pos_df$pos.label,cols] & ans[pos_df$posl.label,cols]
    NAs <- is.na(ansl)
    ansl[NAs] <- FALSE
    ansr <- ans
    ansr[pos_df$pos.label,cols] <-
        ans[pos_df$pos.label,cols] & ans[pos_df$posr.label,cols]
    NAs <- is.na(ansr)
    ansr[NAs] <- FALSE
    ans.or <- ans
    ans.or[pos_df$pos.label,cols] <- ansl[pos_df$pos.label,cols] |
        ansr[pos_df$pos.label,cols]
    rownames(ans.or) <- pos_df[rownames(ans.or),"orig"]
    return(ans.or)
}

#' Find probe hits with a consecutive probe or another sample
#'
#' @param hit_mat matrix of logical values that indicate a hit with a
#' TRUE value
#'
#' @return matrix of logical values indicate that the TRUE hit is supported by
#' a consecutive probe hit in the sample sample or the within another sample
probeHitSupported<-function(hit_mat) {
    probes <- rownames(hit_mat)
    proteins <- getProteinLabel(probes)
    positions <- getProteinStart(probes)
    tiling <- getProteinTiling(probes)
    hit_df <- as.data.frame(hit_mat, stringsAsFactors=FALSE)
    cols <- seq_len(ncol(hit_df))
    hit_df$Protein <- proteins
    hit_df$Pos <- positions
    hit_df_protein <- split(hit_df, hit_df$Protein)
    hit_df_protein <- lapply(
        hit_df_protein,
        probeHitSupportedSingle,
        tiling = tiling,
        cols = cols
    )
    ans_dt <- data.table::rbindlist(hit_df_protein)
    ans_df <- as.data.frame(ans_dt, stringsAsFactors=FALSE)
    rownames(ans_df) <- paste0(ans_df$Protein,";",ans_df$Pos)
    ans_df <- ans_df[probes, cols]  ##Reorder rows to the original matrix
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
    one_probe_epitopes <- oneProbeEpitopes(rownames(sample_epitopes))
    one_call_epitopes <- rowSums(sample_epitopes) == 1
    one_hit_epitopes <- one_probe_epitopes & one_call_epitopes
    ans <- rownames(sample_epitopes)[one_hit_epitopes]
    return(ans)
}

#' Make Epitope Calls
#'
#' @param epi_ds HERONEpitopeDataSet with pvalue assay
#' @param padj_cutoff p-value cutoff to use
#' @param one_hit_filter filter one hit epitopes?
#'
#' @return HERONEpitopeDataSet with calls assay added
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' seq_pval_res <- calcCombPValues(heffron2021_wuhan)
#' pr_pval_res <- convertSequenceDSToProbeDS(seq_pval_res)
#' pr_calls_res <- makeProbeCalls(pr_pval_res)
#' epi_segments_uniq_res <- findEpitopeSegments(
#'     PDS_obj = pr_calls_res,
#'     segment_method = "unique"
#' )
#' epi_padj_uniq <- calcEpitopePValues(
#'     probe_pds = pr_calls_res,
#'     epitope_ids = epi_segments_uniq_res,
#'     metap_method = "wilkinsons_max1"
#' )
#' makeEpitopeCalls(epi_padj_uniq)
makeEpitopeCalls<-function(
        epi_ds,
        padj_cutoff = 0.05,
        one_hit_filter = TRUE) {

    stopifnot(is(epi_ds, "HERONEpitopeDataSet"))
    stopifnot("pvalue" %in% assayNames(epi_ds))
    stopifnot("padj" %in% assayNames(epi_ds))

    res <- makeCalls(se = epi_ds, padj_cutoff = padj_cutoff)
    if (one_hit_filter) {
        ohe <- oneHitEpitopes(assay(res, "calls"))
        padj <- assay(res, "padj")
        padj[ohe,] <- 1.0
        assay(res, "padj") <- padj
        calls <- assay(res, "calls")
        calls[ohe,] <- FALSE
        assay(res, "calls") <- calls
        metadata(res)$one_hit_epitopes <- ohe
        metadata(res)$one_probe_epitopes <- oneProbeEpitopes(rownames(epi_ds))
    }
    return(res)
}

#' Get K of N statistics from an experiment with padj and calls
#'
#' Calculates the number of samples (K), the frequency
#' of samples (F), and the percentage of samples (P) called.
#' If the colData DataFrame contains a condition column with at least two
#' conditions, then a K, F, and P is calculated for each condition and the
#' results are reported as separate columns.
#'
#' @param obj HERON Dataset with a "calls" assay
#'
#' @return DataFrame with K (#calls), F (fraction calls), P (%Calls)
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' seq_pval_res <- calcCombPValues(heffron2021_wuhan)
#' pr_pval_res <- convertSequenceDSToProbeDS(seq_pval_res)
#' pr_calls_res <- makeProbeCalls(pr_pval_res)
#' getKofN(pr_calls_res)
#' @importFrom S4Vectors DataFrame
getKofN<-function(obj) {
    stopifnot(is(obj, "SummarizedExperiment"))
    stopifnot("calls" %in% assayNames(obj))

    calls <- assay(obj, "calls")
    col_data <- colData(obj)
    post_idx <- tolower(col_data$visit) == "post"
    post_cols <- rownames(col_data)[post_idx]

    K_post <- rowSums(calls[,post_cols])
    np <- length(post_cols)
    k_of_n <- DataFrame(
        "K" = K_post,
        "F" = K_post / np,
        "P" = K_post / np * 100.0,
        row.names = rownames(calls)
    )

    if ("condition" %in% colnames(col_data)) {
        uconds <- unique(col_data$condition[post_idx])
        if (length(uconds) > 1) {
            for (cond in uconds) {
                cond_idx <- post_idx & col_data$condition == cond
                cond_cols <- rownames(col_data)[cond_idx]
                nc <- length(cond_cols)
                Kc <- rowSums(calls[,cond_cols])
                k_of_n[,paste0("K_",cond)] <- Kc
                k_of_n[,paste0("F_",cond)] <- Kc / nc
                k_of_n[,paste0("P_",cond)] <- Kc / nc * 100.0
            }
        }
    }
    return(k_of_n)
}

#' Making Probe-level Calls
#'
#' \code{makeProbeCalls} returns call information on a HERONProbeDataSet
#' using the "padj" assay
#'
#' @param pds HERONProbeDataSet with the "padj" assay
#' @param padj_cutoff cutoff to use
#' @param one_hit_filter filter out one-hit probes?
#'
#' @return HERONProbeDataSet with the "calls" assay added
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' pval_seq_res <- calcCombPValues(heffron2021_wuhan)
#' pval_probe_res <- convertSequenceDSToProbeDS(pval_seq_res)
#' calls_res <- makeProbeCalls(pval_probe_res)
makeProbeCalls<-function(pds, padj_cutoff = 0.05, one_hit_filter = TRUE) {
    stopifnot(is(pds, "HERONProbeDataSet"))
    stopifnot("padj" %in% assayNames(pds))
    res <- makeCalls(se = pds, padj_cutoff = padj_cutoff)
    if (one_hit_filter) {
        ohp <- oneHitProbes(assay(res, "calls"))
        ## Set the padj values to 1.0 and set probe call to FALSE.
        padj <- assay(res, "padj")
        padj[rownames(padj) %in% ohp,] <- 1.0
        assay(res, "padj") <- padj
        calls <- assay(res, "calls")
        calls[rownames(calls) %in% ohp,] <- FALSE
        assay(res, "calls") <- calls
        metadata(res)$one_hit_probes <- ohp
    }
    return(res)
}

#' Make Protein-level Calls
#'
#' @param prot_ds HERONProteinDataSet with the "padj" assay
#' @param padj_cutoff cutoff to use
#' @param one_hit_filter use the one-hit filter?
#'
#' @return HERONProteinDataSet with the "calls" assay added
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' seq_pval_res <- calcCombPValues(heffron2021_wuhan)
#' pr_pval_res <- convertSequenceDSToProbeDS(seq_pval_res)
#' pr_calls_res <- makeProbeCalls(pr_pval_res)
#' epi_segments_uniq_res <- findEpitopeSegments(
#'     PDS_obj = pr_calls_res,
#'     segment_method = "unique"
#' )
#' epi_padj_uniq <- calcEpitopePValues(
#'     probe_pds = pr_calls_res,
#'     epitope_ids = epi_segments_uniq_res,
#'     metap_method = "wilkinsons_max1"
#' )
#' prot_padj_uniq <- calcProteinPValues(
#'     epitope_ds = epi_padj_uniq,
#'     metap_method = "tippetts"
#' )
#' prot_calls <- makeProteinCalls(prot_padj_uniq)
makeProteinCalls<-function(prot_ds, padj_cutoff = 0.05, one_hit_filter = FALSE){
    stopifnot(is(prot_ds, "HERONProteinDataSet"))
    stopifnot("padj" %in% assayNames(prot_ds))
    res <- makeCalls(prot_ds, padj_cutoff = padj_cutoff)
    if (one_hit_filter) {
        k_of_n <- getKofN(res)
        k1_protein <- rownames(k_of_n)[k_of_n$K == 1]
        eids <- metadata(prot_ds)$epitope_ids
        one_epi_one_probe_prots <- oneEpiOneProbeProts(eids)
        ohp <- intersect(k1_protein, one_epi_one_probe_prots)
        padj <- assay(res, "padj")
        padj[rownames(padj) %in% ohp,] <- 1.0
        assay(res, "padj") <- padj
        calls <- assay(res, "calls")
        calls[rownames(calls) %in% ohp,] <- FALSE
        assay(res, "calls") <- calls
        metadata(prot_ds)$one_hit_proteins <- ohp
    }
    return(res)
}

oneEpiProts<-function(epitope_ids) {
    eids_protein <- getEpitopeProtein(epitope_ids)
    eids_protein_tbl <- table(eids_protein)
    one_epi_prot <- names(eids_protein_tbl)[eids_protein_tbl == 1]
    return(one_epi_prot)
}

oneEpiOneProbeProts<-function(epitope_ids) {
    one_epi_prot <- oneEpiProts(epitope_ids)
    one_probe_epi <- epitope_ids[oneProbeEpitopes(epitope_ids)]
    one_probe_epi_prot <- getEpitopeProtein(one_probe_epi)
    one_epi_one_probe_prots <- intersect(one_epi_prot, one_probe_epi_prot)
    return(one_epi_one_probe_prots)
}



#' Make calls on an input matrix of p-adjusted values
#'
#' This function takes in a matrix of p-values and using
#' a cutoff, finds the calls on the column/sample sample-level.
#'
#' @param se SummarizedExperiment with the "padj" assay
#' @param padj_cutoff cutoff to use
#'
#' @return SummarizedExperiment with the "calls" assay added
#' @noRd
makeCalls<-function(se, padj_cutoff = 0.05) {
    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot("padj" %in% assayNames(se))
    padj_mat <- assay(se, "padj")
    padj_mat[is.na(padj_mat)] <- 1.0
    calls <- padj_mat < padj_cutoff
    assay(se, "calls") <- calls
    return(se)
}

