#' Calculate Probe-level p-values
#'
#' calculates p-values on the matrix (can be sequence as well), using either a
#' z-score global, t-stat differential, or a combination of both.
#' input a matrix that has samples on the columns and sequence probes on the
#' rows.
#'
#' @param probe_mat numeric matrix, rows as seqs and columns as samples
#' @param colData_in experiment design data.frame
#' @param t_sd_shift sd shift multiplier for diff t-test
#' @param t_abs_shift abs shift to use for the diff t-test
#' @param t_paired do paired t-test
#' @param z_sd_shift sd shift multiplier for global z-test
#' @param use which p-value method(s) to use
#' @param p_adjust_method method for adjusting p-values
#'
#' @return matrix of calculated p-values
#' @noRd
calcProbePValuesProbeMat<-function(
        probe_mat,
        colData_in,
        t_sd_shift = NA,
        t_abs_shift = NA,
        t_paired = FALSE,
        z_sd_shift=0,
        use = "both",
        p_adjust_method = "BH") {

    if (use == "both") { use_t <- TRUE; use_z <- TRUE; use_c <- TRUE }
    else if (use == "t") { use_t <- TRUE; use_z <- FALSE; use_c <- FALSE }
    else if (use == "z") { use_t <- FALSE; use_z <- TRUE; use_c <- FALSE }
    else { stop("Unknown use paramater:" , use) }
    if (missing(colData_in) || is.null(colData_in)) {
        c_mat <- probe_mat
    } else {
        c_mat <- probe_mat[,colData_in$TAG]
    }
    if (use_t) {
        if (t_paired) {
            pvaluet_df <- calcProbePValuesTPaired(
                probe_mat = c_mat, colData_in = colData_in,
                sd_shift = t_sd_shift, abs_shift = t_abs_shift)
        } else {
            pvaluet_df <- calcProbePValuesTUnpaired(
                c_mat, colData_in = colData_in,
                sd_shift = t_sd_shift,
                abs_shift = t_abs_shift)
        }
        praw <- pvaluet_df
    }
    if (use_z) {
        pvaluez_df <- calcProbePValuesZ(c_mat, colData_in, sd_shift = z_sd_shift)
        praw <- pvaluez_df
    }
    if (use_c) {
        praw <- combinePValueMatrix(pvaluet_df, pvaluez_df)
    }
    praw[is.na(praw)] <- 1.0 # Conservatively set NAs to p-value = 1
    padj_df <- p_adjust_mat(praw, method = p_adjust_method)

    ans <- padj_df
    attr(ans, "pvalue") <- praw
    attr(ans, "c_mat") <- c_mat
    attr(ans, "colData") <- colData_in
    if (use_t) { attr(ans, "t") <- pvaluet_df }
    if (use_z) { attr(ans, "z") <- pvaluez_df }
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

#' Calculate p-values using the "exprs" assay
#'
#' @param obj HERONSequenceDataSet or HERONProbeDataSet
#' @param colData_in optional column DataFrame (default: NULL => colData(obj)))
#' @param t_sd_shift standard deviation shift for differential test
#' @param t_abs_shift absolute shift for differential test
#' @param t_paired run paired analysis
#' @param z_sd_shift standard deviation shift for global test
#' @param use use global-test ("z"), differential-test ("t"), or both ("both")
#' @param p_adjust_method method for adjusting p-values
#'
#' @return HERONSequenceDataSet/HERONProbeDataSet with the pvalue assay added
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' seq_pval_res <- calcCombPValues(heffron2021_wuhan)
calcCombPValues<-function(
    obj,
    colData_in = NULL,
    t_sd_shift = NA,
    t_abs_shift = NA,
    t_paired = FALSE,
    z_sd_shift = 0,
    use = "both",
    p_adjust_method = "BH"
) {
    stopifnot(is(obj, "HERONSequenceDataSet") || is(obj, "HERONProbeDataSet"))
    if (is.null(colData_in)) {colData_in <- colData(obj)}
    pval <- calcProbePValuesProbeMat(
        probe_mat = assay(obj, "exprs"),
        colData_in = colData_in,
        t_sd_shift = t_sd_shift,
        t_abs_shift = t_abs_shift,
        t_paired = t_paired,
        z_sd_shift = z_sd_shift,
        use = use,
        p_adjust_method = "none"
    )

    res <- addPValues(obj, pval)
    res <- p_adjust_ds(res)
    return(res)
}

#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assay<-
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' @noRd
addPValues<-function(obj, pval) {
    res <- obj
    if (all(dim(assay(obj, "exprs")) == dim(pval)) &&
        all(colnames(assay(obj,"exprs")) == colnames(pval)) &&
        all(rownames(assay(obj,"exprs")) == rownames(pval))
    ) {
        assay(res, "pvalue") <- pval
    } else {
        # Reduced rows/columns
        exprs_old <- assay(obj, "exprs")
        colData_old <- colData(obj)

        exprs_new <- exprs_old[rownames(pval), colnames(pval)]
        colData_new <- colData_old[(colnames(pval)),]
        if (is(obj, "HERONSequenceDataSet")) {
            res <- HERONSequenceDataSet(
                exprs = exprs_new,
                colData = colData_new
            )
            metadata(res) <- metadata(obj)
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
#' @param colData_in design data.frame
#' @param sd_shift standard deviation shift to use when calculating p-values
#' (default 0)
#' @param all calculate p-values for all samples, otherwise report the p-values
#' for just the post samples
#'
#' @return matrix of "p-values"
#'
#' @examples
#' data(heffron2021_wuhan)
#' colData_wu <- colData(heffron2021_wuhan)
#' pval_res <- calcProbePValuesZ(assay(heffron2021_wuhan), colData_wu)
#' @noRd
calcProbePValuesZ<-function(
        probe_mat,
        colData_in,
        sd_shift = 0,
        all = FALSE
) {

    if (all || missing(colData_in) || is.null(colData_in)) {
        message("No colData or all asked for, calculating on all columns")
        all_cols <- colnames(probe_mat)
        post_cols <- all_cols
    } else {
        pre_cols <- colData_in$TAG[tolower(colData_in$visit)=="pre"]
        post_cols <- colData_in$TAG[tolower(colData_in$visit) == "post"]
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


getPairedMapping<-function(colData) {
    pre_df <- colData[tolower(colData$visit) =="pre",]
    post_df <- colData[tolower(colData$visit) == "post",]

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
#' @param colData_in design data.frame
#' @param sd_shift standard deviation shift to use when calculating p-values.
#' Either sd_shift or abs_shift should be set
#' @param abs_shift absolute shift to use when calculating p-values.
#' @param debug print debugging information
#'
#' @return matrix of p-values on the post columns defined in the colData matrix.
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
#' colData_wu <- colData(heffron2021_wuhan)
#' pre_idx = which(colData_wu$visit == "pre")
#' ## Make some samples paired
#' colData_post = colData_wu[colData_wu$visit == "post",]
#' new_ids = colData_post$Sample_ID[seq_len(5)]
#' colData_wu$ptid[pre_idx[seq_len(5)]] = new_ids
#' exprs <- assay(heffron2021_wuhan, "exprs")
#' pval_res <- calcProbePValuesTPaired(exprs, colData_wu)
calcProbePValuesTPaired <- function(
        probe_mat,
        colData_in,
        sd_shift = NA,
        abs_shift = NA,
        debug = FALSE
) {
    if (!is.na(sd_shift) && !is.na(abs_shift)) {
        stop("Either sd or abs can be set. Not both.")
    }
    mapping <- getPairedMapping(colData_in)
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
#' @param colData_in design data.frame
#' @param sd_shift standard deviation shift to use when calculating p-values
#' Either sd_shift or abs_shift should be set
#' @param abs_shift absolute shift to use when calculating p-values
#'
#' @return matrix of p-values on the post columns defined in the colData matrix
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' colData_wu <- colData(heffron2021_wuhan)
#' pval_res <- calcProbePValuesTUnpaired(assay(heffron2021_wuhan), colData_wu)
calcProbePValuesTUnpaired<-function(
        probe_mat,
        colData_in,
        sd_shift=NA,
        abs_shift=NA
) {
    if (!is.na(sd_shift) && !is.na(abs_shift)) {
        stop("Either sd or abs can be set, not both.")
    }
    pre_cols <- colData_in$TAG[tolower(colData_in$visit) =="pre"]
    post_cols <- colData_in$TAG[tolower(colData_in$visit) == "post"]

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
    colnames(post_tv) <- post_cols

    ans <- matrix(NA, nrow = nrow(probe_mat),ncol=length(post_cols))
    rownames(ans) <- rownames(probe_mat)
    colnames(ans) <- post_cols


    for (post_col in post_cols) {
        ans[,post_col] <-
            stats::pt(q=post_tv[,post_col], df=dfree, lower.tail=FALSE)
    }
    ans <- as.data.frame(ans,stringsAsFactors=FALSE)

    pars <- cbind(pars, post_tv)
    attr(ans, "pars") <- pars
    return(ans)
}

#' Calculate protein-level p-values
#'
#' @param epitope_ds HERONEpitopeDataSet with the "pvalue" assay
#' @param metap_method meta p-value method to use
#' @param p_adjust_method p.adjust method to use
#'
#' @details
#' see calcEpitopePValues for a list of meta p-value methods supported
#' by HERON. the protein should be one that requires at least one of the
#' epitope p-values to be small (e.g. wmax1).
#'
#' @return HERONProteinDataSet with the "pvalue" and "padj" assays
#' @export
#'
#' @seealso [stats::p.adjust()] for p_adjust_parameter.
#' @seealso [calcEpitopePValues()] for meta p-value methods
#'
#' @examples
#' data(heffron2021_wuhan)
#' pval_seq_res <- calcCombPValues(heffron2021_wuhan)
#' pval_pr_res <- convertSequenceDSToProbeDS(pval_seq_res)
#' calls_res <- makeProbeCalls(pval_pr_res)
#' segments_res <- findEpitopeSegments(calls_res, "unique")
#' epval_res <- calcEpitopePValues(calls_res, segments_res)
#' ppval_res <- calcProteinPValues(epval_res)
calcProteinPValues<-function(
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
#' @noRd
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
#' @param probe_pds HERONProbeDataSet with the "pvalue" assay
#' @param epitope_ids vector of epitope ids
#' @param metap_method meta p-value method to use (see below)
#' @param p_adjust_method what p.adjust method to use.
#'
#' @details
#' The meta p-value methods supported by \code{calcEpitopePValues} are:
#' min_bonf*,
#' min*,
#' max*,
#' fischer/sumlog,
#' hmp/harmonicmeanp,
#' wilkinsons_min1/tippets,
#' wilkinsons_min2/wmin2,
#' wilkinsons_min3,
#' wilkinsons_min4,
#' wilkinsons_min5,
#' wilkinsons_max1/wmax1,
#' wilkinsons_max2/wmax2,
#' and cct.
#'
#' When choosing a p-value method, keep in mind that the epitope p-value should
#' be one that requires most of the probe p-values to be small (e.g. *wmax1*)
#' Other p-value methods such as the*cct* and the *hmp* have been shown to be
#' more accurate with p-value that have dependencies.
#'
#' @seealso [stats::p.adjust()] for p_adjust_parameter.
#'
#' @return HERONEpitopeDataSet with "pvalue" and "padj" assays
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' pval_seq_res <- calcCombPValues(heffron2021_wuhan)
#' pval_pr_res <- convertSequenceDSToProbeDS(pval_seq_res)
#' calls_res <- makeProbeCalls(pval_pr_res)
#' segments_res <- findEpitopeSegments(calls_res, "unique")
#' epval_res <- calcEpitopePValues(calls_res, segments_res)
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assays<-
calcEpitopePValues<-function(
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
#'
#' @examples
#' mat = matrix(runif(25), nrow=5)
#' rownames(mat) = paste0("A;",seq_len(nrow(mat)))
#' p_adjust_mat(mat)
#' @noRd
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
#' @return SummarizedExperiment with the "padj" assay added
#' @noRd
p_adjust_ds <- function(obj, method = "BH") {
    stopifnot(is(obj, "SummarizedExperiment"))
    stopifnot("pvalue" %in% assayNames(obj))
    assay(obj, "padj") <- p_adjust_mat(assay(obj, "pvalue"), method = method)
    return(obj)
}

