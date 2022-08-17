
#' getOverlapClusters
#'
#' This function finds all clusters where at least one-probe is called consecutively across the protein
#' A way of filtering the clusters to concentrate on areas where there is significant signal.
#' @param sample_probes from makeProbeCalls
#' @param protein_tiling tiling for each protein
#'
#' @return data.frame with
#' @export
#' @examples
#' getOverlapClusters
getOverlapClusters<-function(sample_probes, protein_tiling = getProteinTiling(rownames(sample_probes))) {
    calls = rowSums(sample_probes) > 0
    blocks = findBlocksProbeT(names(calls)[calls], getProteinTiling(names(calls)))
    blocks$Epitope_ID = rownames(blocks);
    return(blocks);
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
#' @return
#' @export
#'
#' @examples
findEpitopes<-function(
        probe_pvalues_res,
        probe_calls,
        segment.method = "unique",
        segment.score.type = "binary",
        segment.dist.method = "hamming",
        segment.cutoff = "silhouette"
) {
    if (segment.method == "unique") {
        message("Getting unique epitope calls");
        probe_sample_padj = probe_calls$probe_sample_padj
        probe_cutoff = probe_calls$probe_cutoff;
        probes = probe_calls$probes;
        proteins = probe_calls$proteins
        positions = probe_calls$positions
        protein_tiling = probe_calls$protein_tiling;
        one_hit_filter = probe_calls$one_hit_filter;
        epitope_calls = makeEpitopeCalls(
            probe_sample_padj,
            probe_cutoff = probe_cutoff,
            probes=probes,
            proteins=proteins,
            pos=positions,
            protein_tiling=protein_tiling,
            one_hit_filter=one_hit_filter
        );
        segments = epitope_calls$all_epitopes$Epitope_ID;
    } else {
        message("Epitopes: Segmenting probes");
        overlap_cluster_df = getOverlapClusters(probe_calls$sample_probes);

        if (segment.score.type == "zscore") {
            probe_sample_pvalues = attr(probe_pvalues_res, "pvalue");
            if (probe_calls$one_hit_filter) #TODO investigate whether to find segments, then filter versus filter then find segments
            {
                probe_sample_pvalues[rownames(probe_sample_pvalues) %in% probe_calls$to_remove,] = 1;
            }
            sample_probe_score = pvalue_to_zscore(probe_sample_pvalues)
        } else {
            sample_probe_score = probe_calls$sample_probes
        }

        cluster_res = getClusterSegments(
            sample_probes = sample_probe_score,
            overlap_cluster_df = overlap_cluster_df,
            method = segment.method,
            dist.method = segment.dist.method,
            cutoff = segment.cutoff
        );
        segments = cluster_res$segments;
        #We shouldn't have to worry about one-hit calls here since we are using the probes that have already
        #been filtered.
    }
    return(segments);
}


#Method can be: hclust, skater, all
getClusterSegments<-function(
        all_epitopes,
        sample_probes,
        overlap_cluster_df = getOverlapClusters(all_epitopes), #Step 1 - from the unique set of epitopes, find all epitope that overlap
        do.plot=FALSE,
        method="hclust",
        dist.method = "orig",
        cutoff = "silhouette"
) {



    #Step 2 - segment clusters to try and maximize tradeoff between k of n and # of probe
    segments = c();

    library(cluster) # For silhouette function

    for (cluster_idx in 1:nrow(overlap_cluster_df)) {
        cluster_id = overlap_cluster_df$Epitope_ID[cluster_idx];
        #cluster_name = getClusterShortName(cluster_id);



        probes = getEpitopeIDsToProbeIDs(overlap_cluster_df$Epitope_ID[cluster_idx]);
        #print(probes$PROBE_ID)
        #Because of tiling, some probes might be missing, so just get the ones that exist in the dataset.
        sample_probes_sub = sample_probes[probes$PROBE_ID[probes$PROBE_ID %in% rownames(sample_probes)],];
        #print(sample_probes_sub)
        if (method == "all") {
            segment_ids = getClusterSegmentsAll(sample_probes_sub, cluster_id);
            segments = c(segments, segment_ids);
        } else {
            if (nrow(probes) <=2 || nrow(sample_probes_sub) <= 2) { #Don't attempt to segment small clusters.
                segments = c(segments,cluster_id)
            } else {
                if (method == "hclust") {
                    segment_ids  = getClusterSegmentsHClust(sample_probes_sub, cluster_id, do.plot=do.plot, cutoff=cutoff, dist.method=dist.method)
                } else if (method == "skater") {
                    #message("Calling skater n:",nrow(sample_probes_sub))
                    segment_ids = getClusterSegmentsSkater(sample_probes_sub, cluster_id, dist.method = dist.method, cutoff=cutoff)
                }
                segments = c(segments, segment_ids);
            }
        }
        if (cluster_idx %% 200 == 0) {
            cat(cluster_idx, " / ",nrow(overlap_cluster_df),"\n");
        }


    }
    #Return the final list of segments.

    res = list();
    res$overlap_clusters = overlap_cluster_df
    res$segments = segments;

    return(res);
}
