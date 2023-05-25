
#' getOverlapClusters
#'
#' This function finds all clusters where at least one-probe is called
#' in a sample consecutively across the tiled protein
#' A way of filtering the clusters to concentrate on areas where there is
#' significant signal.
#' @param sample_probes logical matrix of probe calls, where rows are probe ids
#' and columns are samples
#'
#' @return data.frame
getOverlapClusters<-function(sample_probes) {
    calls <- rowSums(sample_probes) > 0
    ncalls <- names(calls)
    blocks <- findBlocksProbeT(ncalls[calls], getProteinTiling(ncalls))
    blocks$Epitope_ID <- rownames(blocks)
    return(blocks)
}


#' Find Epitopes from probe stats and calls.
#'
#' @param PDS_obj HERONProbeDataSet with pvalues and calls in the assay
#' @param segment_method which epitope finding method to use
#' (binary or zscore, applies for hclust or skater)
#' @param segment_score_type which type of scoring to use for probes
#' @param segment_dist_method what kind of distance score method to use
#' @param segment_cutoff for clustering methods, what cutoff to use
#' (either numeric value or 'silhouette')
#'
#' @return a vector of epitope identifiers or segments found
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' probe_meta <- attr(heffron2021_wuhan, "probe_meta")
#' pval_res <- calcProbePValuesSeqMat(heffron2021_wuhan, probe_meta)
#' calls_res <- makeProbeCalls(pval_res)
#' segments_res = findEpitopeSegments(probe_calls = calls_res)
findEpitopeSegmentsPDS<-function(
    PDS_obj,
    segment_method,
    segment_score_type,
    segment_dist_method,
    segment_cutoff
) {
    #Check PDS_obj is a HERONProbeDataSet
    if (!is(PDS_obj, "HERONProbeDataSet")) {
        stop("HERONProbeDataSet object required")
    }
    #Check if pvalue are in the assays
    if (!("pvalue" %in% assayNames(PDS_obj))) {
        stop("Probe p-values not found")
    }
    #Check if calls are in the assays
    if (!("calls" %in% assayNames(PDS_obj))) {
        stop("Probe calls not found")
    }
    calls <- assay(PDS_obj, "calls")

    segments <- c()
    if (segment_method == "unique") {
        segments <- findEpitopeSegmentsUnique(calls)
    } else {
        overlap_cluster_df <- getOverlapClusters(calls)
        if (segment_score_type == "zscore") {
            probe_sample_pvalues <- assay(PDS_obj, "pvalue")
            sample_probe_score <- pvalue_to_zscore(probe_sample_pvalues)
        } else {
            sample_probe_score <- calls
        }

        segments <- getClusterSegments(
            sample_probes = sample_probe_score,
            overlap_cluster_df = overlap_cluster_df,
            method = segment_method,
            dist.method = segment_dist_method,
            cutoff = segment_cutoff
        )
    }

    return(segments)
}


#' Find Epitopes from probe stats and calls.
#'
#' @param probe_pvalues_res results from calcProbePValuesProbeMat
#' @param probe_calls results from makeProbeCalls
#' @param segment.method which epitope finding method to use
#' @param segment.score.type which type of scoring to use for probes
#' (binary or zscore, applies for hclust or skater)
#' @param segment.dist.method what kind of distance score method to use
#' (hamming, euclidean, ..., applies for hclust or skater)
#' @param segment.cutoff for clustering methods, what cutoff to use
#' (either numeric value or 'silhouette')
#'
#' @return a vector of epitope identifiers or segments found.
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' probe_meta <- attr(heffron2021_wuhan, "probe_meta")
#' pData <- attr(heffron2021_wuhan, "pData")
#' pval_res <- calcProbePValuesSeqMat(heffron2021_wuhan, probe_meta, pData)
#' calls_res <- makeProbeCalls(pval_res)
#' segments_res = findEpitopeSegments(probe_calls = calls_res)
findEpitopeSegments<-function(
        probe_pvalues_res,
        probe_calls,
        segment.method = "unique",
        segment.score.type = "binary",
        segment.dist.method = "hamming",
        segment.cutoff = "silhouette"
) {
    if (segment.method == "unique") {
        segments <- findEpitopeSegmentsUnique(probe_calls$sample)
    } else {
        overlap_cluster_df <- getOverlapClusters(probe_calls$sample)
        if (segment.score.type == "zscore") {
            if ("pvalue" %in% names(attributes(probe_pvalues_res))) {
                probe_sample_pvalues <- attr(probe_pvalues_res, "pvalue")
            } else {
                probe_sample_pvalues <- probe_pvalues_res
            }
            if (probe_calls$one_hit_filter)
            {
                ##Make sure the one hit calls are filtered out
                probe_sample_pvalues[rownames(probe_sample_pvalues) %in%
                    probe_calls$to_remove,] <- 1
            }
            sample_probe_score <- pvalue_to_zscore(probe_sample_pvalues)
        } else {
            sample_probe_score <- probe_calls$sample
        }

        segments <- getClusterSegments(
            sample_probes = sample_probe_score,
            overlap_cluster_df = overlap_cluster_df,
            method = segment.method,
            dist.method = segment.dist.method,
            cutoff = segment.cutoff
        )
    }
    return(segments)
}

#' Get the Cluster Segments from a starting list of epitopes
#'
#' @param all_epitopes list of epitopes to further segment
#' @param sample_probes from makeProbeCalls(...), logical matrix of probe calls
#' @param do.plot make supplementary plots?
#' @param method clustering method to use (hclust, skater, all)
#' @param dist.method for clustering methods, what distance metric to use
#' ("hamming", "euclidean")
#' @param cutoff what clustering cutoff to use (silhouette or numeric value)
#' @param overlap_cluster_df result from getOverlapClusters
#'
#' @return list of results
#' overlap_clusters - overlap_cluster_df parameter
#' segments - list of epitope identifiers of the segments
getClusterSegments<-function(
        all_epitopes,
        sample_probes,
        #Step 1 - from the unique set of epitopes, find all epitope that overlap
        do.plot=FALSE,
        method="hclust",
        dist.method = "orig",
        cutoff = "silhouette",
        overlap_cluster_df

) {
    segments <- c()
    for (c_idx in seq_len(nrow(overlap_cluster_df))) {
        cluster_id <- overlap_cluster_df$Epitope_ID[c_idx]
        probes <- getEpitopeIDsToProbeIDs(overlap_cluster_df$Epitope_ID[c_idx])
        #Because of tiling, some probes might be missing,
        sample_probes_sub <-
            sample_probes[probes$PROBE_ID[probes$PROBE_ID %in%
                rownames(sample_probes)],]
        if (nrow(probes) <=2 || nrow(sample_probes_sub) <= 2) {
            #Don't attempt to segment small clusters.
            segments <- c(segments,cluster_id)
        } else {
            if (method == "hclust") {
                segment_ids <- getClusterSegmentsHClust(
                    sample_probes = sample_probes_sub,
                    cluster_id = cluster_id,
                    do.plot=do.plot,
                    cutoff=cutoff,
                    dist.method=dist.method
                )
            } else if (method == "skater") {
                segment_ids <- getClusterSegmentsSkater(
                    sample_probes_sub = sample_probes_sub,
                    cluster_id = cluster_id,
                    dist_method = dist.method,
                    cutoff=cutoff
                )
            }
            segments <- c(segments, segment_ids)
        }
    }
    return(segments)
}


#' Find unique set of epitope regions across all sample calls.
#'
#' @param probe_sample_calls probe hit matrix
#'
#' @return vector of epitope seqments
findEpitopeSegmentsUnique<-function(probe_sample_calls) {
    segments <- c()
    protein_tiling <- getProteinTiling(rownames(probe_sample_calls))

    probe_sample_list <- getSampleProbesList(probe_sample_calls)

    epitope_sample_list <- lapply(
        probe_sample_list,
        findBlocksProbeT,
        protein_tiling = protein_tiling
    )
    segments <- unique(unlist(lapply(epitope_sample_list, rownames)))
    return(segments)
}


getSampleProbesList<-function(probe_calls)
{
    ans <- list()
    for (col_idx in seq_len(ncol(probe_calls))) {
        probes <- rownames(probe_calls)[probe_calls[,col_idx]]
        if (length(probes) > 0) {
            ans[[colnames(probe_calls)[col_idx]]] <- probes
        }
    }
    return(ans)
}

#' Find Blocks of consecutive probes
#'
#' This function will find blocks of consecutive probes within the passed
#' probe parameter
#' @param probes vector of probe identifiers of the format
#' c(Prot1;1, ... Prot1;10)
#' @param protein_tiling tiling of the associated proteins
#' @param proteins associated proteins to probes (cache speed up)
#' @param starts associated starts from probes (cache speed up)
#'
#' @return data.frame with the Protein, Start, Stop, and Number.Of.Probes
#' columns
#' @export
#'
#' @examples
#' findBlocksProbeT(c("A;1", "A;2", "A;3", "B;2", "B;3", "C;10", "A;5", "A;6"))
findBlocksProbeT<-function(
        probes,
        protein_tiling,
        proteins = getProteinLabel(probes),
        starts = getProteinStart(probes)
) {
    if (missing(protein_tiling)) {
        protein_tiling <- 1
    }
    protein.df <- data.frame(
        Protein = proteins,
        Pos = starts,
        stringsAsFactors=FALSE
    )

    protein_list <- split(protein.df, protein.df$Protein)

    ans_l <- lapply(protein_list, findBlocksT, protein_tiling = protein_tiling)
    ans_dt <- data.table::rbindlist(ans_l)
    ans_df <- as.data.frame(ans_dt, stringsAsFactors=FALSE)

    ans_df$ProbeSet.AA.Span <- ans_df$Stop - ans_df$Start + 1
    rownames(ans_df) <- getEpitopeID(ans_df$Protein, ans_df$Start, ans_df$Stop)
    return(ans_df)

}

block_df<-function(protein, start, stop, nprobes) {
    return(
        data.frame(
            Protein = protein,
            Start = start,
            Stop = stop,
            Number.Of.Probes = nprobes,
            stringsAsFactors = FALSE
        )
    )
}

getTiling<-function(protein_tiling, protein) {
    if (missing(protein_tiling) || length(protein_tiling) == 0) {
        tiling <- 1
    } else if (length(protein_tiling) == 1) {
        tiling <- protein_tiling[1]
    } else {
        tiling <- protein_tiling[protein]
    }
    return(tiling)
}


#' Find consecutive probes
#'
#' @param prot_df data.frame with the Protein and Starting position of the
#' probe
#'
#' @param protein_tiling tiling for information for each protein
#'
#' @return data.frame with the Protein, Start, Stop, and Number.Of.Probes
#' columns
#' @export
#'
#' @examples
#' probes = c("A;1","A;2","A;3", "A;5","A;6", "A;8")
#' prot_df = data.frame(
#'     Protein = getProteinLabel(probes),
#'     Pos = getProteinStart(probes)
#' )
#' findBlocksT(prot_df)
findBlocksT<-function(prot_df, protein_tiling) {
    if (is.null(prot_df) || nrow(prot_df) == 0) {
        return(NULL)
    }
    if (nrow(prot_df)  == 1) {
        return(block_df(prot_df$Protein, prot_df$Pos, prot_df$Pos, 1))
    }
    tiling <- getTiling(protein_tiling, prot_df$Protein[1])
    prot_df <- prot_df[order(prot_df$Pos, decreasing=FALSE),]
    ans_list <- list()
    start_idx <- 1
    start_pos <- prot_df$Pos[start_idx]
    for (idx in seq.int(2, nrow(prot_df))) {
        current_idx <- idx
        current_pos <- prot_df$Pos[current_idx]
        prev_idx <- idx-1
        prev_pos <- prot_df$Pos[prev_idx]
        d <- current_pos - prev_pos
        if (d > tiling) {
            nprobes <- length(seq(from=start_pos, to=prev_pos, by=tiling))
            ans_list[[length(ans_list)+1]] <-
                block_df(prot_df$Protein[1], start_pos, prev_pos, nprobes)
            start_idx <- current_idx
            start_pos <- current_pos
        }
    }
    if (start_idx != nrow(prot_df)) {
        end_pos <- prot_df$Pos[nrow(prot_df)]
        nprobes <- length(seq(from=start_pos, to=end_pos, by=tiling))
        ans_list[[length(ans_list)+1]] <-
            block_df(prot_df$Protein[1], start_pos, end_pos, nprobes)
    } else { #Last probe is by itself
        ans_list[[length(ans_list)+1]] <-
            block_df(
                prot_df$Protein[1],
                prot_df$Pos[nrow(prot_df)],
                prot_df$Pos[nrow(prot_df)], 1
            )
    }
    ans.dt <- data.table::rbindlist(ans_list)
    ans.df <- as.data.frame(ans.dt)
    return(ans.df)
}


getHClustSilouette<-function(dist_mat2, hc) {
    n <- max(2,(ncol(dist_mat2)-1))
    k_vec <- seq.int(2, n)
    sil_df <- data.frame(
        K = k_vec,
        SIL = rep(NA, n-1),
        stringsAsFactors = FALSE
    )
    for (k in k_vec) {
        hc_cut <- stats::cutree(hc, k = k)

        silhouette.results <- cluster::silhouette(
            list(clustering=hc_cut),
            stats::as.dist(dist_mat2)
        )
        sil.sum <- summary(silhouette.results)
        sil.mean <- sil.sum$si.summary["Mean"]
        sil_df$SIL[k-1] <- as.vector(sil.mean)
    }
    return(sil_df)
}



#' Find Cluster Segments using hclust
#'
#' @param sample_probes logical or numerical (z-score) matrix of probe calls
#' which is just the probes from the cluster_id column
#' @param cluster_id epitope identifier of cluster to further segment
#' @param do.plot make plots of result
#' @param cutoff cutoff to use when finding the number of clusters
#' @param dist.method distance algorithm to use
#' @param dist_mat2 calculated using getHClustDist(cached result)
#' @param hc calculated using hclust (cached result)
#'
#' @return vector of segmented epitope identifiers
getClusterSegmentsHClust <-function(
        sample_probes,
        cluster_id,
        do.plot = FALSE,
        cutoff = "silhouette",
        dist.method = "hamming",
        dist_mat2 = getHClustDistMat(sample_probes, dist.method=dist.method),
        hc=stats::hclust(stats::as.dist(dist_mat2), method="complete")
) {
    if (do.plot) {plot(hc, cex=0.8, main=cluster_id)}
    sil_df <- NULL
    #Find the number of clusters using the silhouette
    if (cutoff == "silhouette") {
        sil_df <- getHClustSilouette(dist_mat2, hc)
        K <- max(sil_df$K[which(sil_df$SIL == max(sil_df$SIL))])
        #Use K that achieve the maximum mean silo score,
        #breaking ties using maximum K (bias toward more clusters).
        hc_cut <- stats::cutree(hc, k= K)
    } else {
        hc_cut <- stats::cutree(hc, h = cutoff)
    }
    if (is.unsorted(hc_cut)) { #Sanity check
        stop("clusters are not ordered correctly")
    }
    #Build the segments from the clusters.
    nsegments <- max(hc_cut)
    cluster_protein <- getEpitopeProtein(cluster_id)
    cluster_start <- getEpitopeStart(cluster_id)
    cluster_stop <- getEpitopeStop(cluster_id)
    cluster_pos <- getProteinStart(rownames(sample_probes))

    segment_ids <- rep("", nsegments)
    for (segment in seq_len(nsegments)) {
        segment_indices <- which(hc_cut == segment)
        segment_pos <- cluster_pos[segment_indices]
        segment_start <- min(segment_pos)
        segment_end <- max(segment_pos)

        segment_id <- getEpitopeID(cluster_protein, segment_start, segment_end)
        segment_ids[segment] <- segment_id
    }
    attr(segment_ids, "sil_df") <- sil_df
    attr(segment_ids, "hc") <- hc
    attr(segment_ids, "dist_mat") <- attr(dist_mat2, "dist_mat")
    attr(segment_ids, "dist_mat2") <- dist_mat2
    return(segment_ids)

}

#' Calculate a distance matrix for use with hclust
#'
#' @param sample_probes_sub matrix of calls
#' @param dist.method distance method to use (see dist)
#'
#' @return matrix that can be converted to a distance for use with complete
#' hierarchical clustering
getHClustDistMat<-function(sample_probes_sub, dist.method) {

    maxl <- nrow(sample_probes_sub)
    max_end <- nrow(sample_probes_sub)
    #Build build 1st level distance
    dist_mat <- NULL

    if (dist.method == "hamming"){
        dist_mat <- hamming(sample_probes_sub) / ncol(sample_probes_sub)
    } else {
        p <- 2
        if (dist.method == "minkowski") {p <- ncol(sample_probes_sub)}
        dist_mat <- as.matrix(
            stats::dist(
                sample_probes_sub,
                method=dist.method,
                diag=TRUE, upper=TRUE,
                p = p
            )
        )
    }
    #Enforce consecutive probes to be clustered together
    dist_mat2 <- 1-diag(1,max_end, max_end)
    for (idx1 in seq_len((max_end-1))) {
        dist_mat2[idx1,idx1] <- 0

        for (idx2 in seq.int(idx1+1, max_end)) {
            slength <- idx2 - idx1 + 1
            if (slength > maxl) {
                current_dist <- 1
            } else {
                sidx <- seq.int(idx1, idx2)
                ut <- dist_mat[sidx, sidx]
                ut <- ut[upper.tri(ut)]
                current_dist <- max(unlist(ut))
            }
            dist_mat2[idx1,idx2] <- current_dist
            dist_mat2[idx2,idx1] <- current_dist
        }
    }
    attr(dist_mat2, "dist_mat") <- dist_mat
    return(dist_mat2)
}



getSkaterGraph<-function(n) {
    edges <- data.frame(
        A = c(seq_len(n-1)),
        B = c(seq_len(n-1)+1)
    )
    edges <- as.matrix(edges)
    return(edges)
}

getSkaterDist<-function(sample_probes_sub_i, dist_method, p) {
    if (dist_method == "hamming") {
        dist_method_new <- function(data, id) {
            return(
                sum(
                    hamming_dist(
                        rbind(
                            colMeans(data[id, , drop = FALSE]),
                            data[id, , drop = FALSE]
                        )
                    )[seq_len(length(id))]
                )
            )
        }
        sk_dist <- hamming_dist(sample_probes_sub_i)
        attr(sk_dist, "dist_method") <- dist_method_new
    } else {
        sk_dist <- stats::dist(sample_probes_sub_i, method=dist_method, p = p)
        attr(sk_dist, "dist_method") <- dist_method
    }
    attr(sk_dist, "p") <- p
    return(sk_dist)
}

getSkaterSilouette<-function(edges, s_p_sub_i, sk_dist) {
    sil_list <- list()
    sk_res <- NULL
    n <- nrow(edges)
    p <- attr(sk_dist, "p")
    dist_method <- attr(sk_dist, "dist_method")
    for (ncuts in seq_len(n-1)) {
        if (is.null(sk_res)) {
            sk_res <- spdep::skater(edges, s_p_sub_i, ncuts=1,
                method = dist_method, p = p)
        } else {
            sk_res <- spdep::skater(sk_res, s_p_sub_i, ncuts=1,
                method = dist_method, p = p)
        }
        sil_res <- cluster::silhouette(list(clustering=sk_res$groups), sk_dist)
        sil_list[[ncuts]] <-
            data.frame(
                ncuts = ncuts,
                K = ncuts + 1,
                SIL = as.vector(summary(sil_res)$si.summary["Mean"]),
                stringsAsFactors=FALSE
            )

    }
    sil_df <- as.data.frame(data.table::rbindlist(sil_list))
    ## Bias towards smaller regions, i.e. more clusters
    NCUTS <- max(sil_df$ncuts[sil_df$SIL == max(sil_df$SIL)])
    sk_res <- spdep::skater(
        edges,
        s_p_sub_i,
        ncuts = NCUTS, method = dist_method, p=p
    )
    attr(sk_res, "sil_df") <- sil_df
    return(sk_res)
}



#' Get Segmentation Using Skater
#'
#' Acceptable dist.methods are "euclidean", "hamming"
#'    "maximum", "manhattan", "canberra", "binary", "minkowski"
#' Minkowski is hardcoded to have a power (p) of number of columns
#' in sample_probes_sub.
#' @param sample_probes_sub matrix of probe calls
#' @param cluster_id cluster identification
#' @param dist_method distance method to use
#' @param cutoff cutoff to use, silhouette supported.
#'
#' @return vector of segments
getClusterSegmentsSkater<-function(
    sample_probes_sub,
    cluster_id,
    dist_method = "euclidean",
    cutoff = "silhouette"
){

    n <- nrow(sample_probes_sub)
    if (n == 1) {
        stop("Need more than 1 point to cluster")
    }
    edges <- getSkaterGraph(n)
    s_p_sub_i <- toNumericMatrix(sample_probes_sub)
    p <- 2
    if (dist_method == "minkowski") {
        p <- ncol(s_p_sub_i)
    }
    sk_dist <- getSkaterDist(s_p_sub_i, dist_method, p)
    if (cutoff == "silhouette") {
        sk_res <- getSkaterSilouette(edges, s_p_sub_i, sk_dist)
    } else {
        stop("cutoff method Not Implemented: ", cutoff )
    }

    nsegments <- max(sk_res$groups)

    cluster_protein <- getEpitopeProtein(cluster_id)
    cluster_start   <- getEpitopeStart(cluster_id)
    cluster_stop    <- getEpitopeStop(cluster_id)
    cluster_pos     <- getProteinStart(rownames(sample_probes_sub))
    segment_ids     <- rep("", nsegments)

    for (segment in seq_len(nsegments)) {
        segment_indices <- which(sk_res$groups == segment)
        segment_pos <- cluster_pos[segment_indices]
        segment_start <- min(segment_pos)
        segment_end <- max(segment_pos)
        segment_id <- getEpitopeID(cluster_protein, segment_start, segment_end)
        segment_ids[segment] <- segment_id
    }
    attr(segment_ids, "sk_res") <- sk_res
    attr(segment_ids, "dist") <- sk_dist
    return(segment_ids)
}




