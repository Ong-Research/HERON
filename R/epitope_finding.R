
#' getOverlapClusters
#'
#' This function finds all clusters where at least one-probe is called
#' in a sample consecutively across the tiled protein
#' A way of filtering the clusters to concentrate on areas where there is
#' significant signal.
#' @param sample_probes from makeProbeCalls(...), logical matrix of probe calls
#' @param protein_tiling tiling for each protein
#'
#' @return data.frame
#' @export
#' @examples
getOverlapClusters<-function(
        sample_probes,
        protein_tiling = getProteinTiling(rownames(sample_probes))
) {
    calls = rowSums(sample_probes) > 0
    ncalls = names(calls)
    blocks = findBlocksProbeT(ncalls[calls], getProteinTiling(ncalls))
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
#' @return a vector of epitope identifiers or segments found.
#' @export
#'
#' @examples
findEpitopeSegments<-function(
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
        epitope_calls = getEpitopeCallsUnique(
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
            if (probe_calls$one_hit_filter)
            {
                #TODO investigate whether to find segments, then filter versus
                #filter then find segments
                probe_sample_pvalues[rownames(probe_sample_pvalues) %in%
                    probe_calls$to_remove,] = 1;
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
        #We shouldn't have to worry about one-hit calls here since we are using
        #the probes that have already been filtered.
    }
    return(segments);
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
#' @export
#'
#' @examples
getClusterSegments<-function(
        all_epitopes,
        sample_probes,
        #Step 1 - from the unique set of epitopes, find all epitope that overlap
        do.plot=FALSE,
        method="hclust",
        dist.method = "orig",
        cutoff = "silhouette",
        overlap_cluster_df = getOverlapClusters(all_epitopes)

) {

    #Step 2 - segment clusters
    segments = c();

    for (c_idx in seq_len(nrow(overlap_cluster_df))) {
        cluster_id = overlap_cluster_df$Epitope_ID[c_idx];
        probes = getEpitopeIDsToProbeIDs(overlap_cluster_df$Epitope_ID[c_idx]);
        #Because of tiling, some probes might be missing,
        #just get the ones that exist in the dataset.
        sample_probes_sub =
            sample_probes[probes$PROBE_ID[probes$PROBE_ID %in%
                rownames(sample_probes)],];

        if (method == "all") {
            #Not sure if we want to support this, can take too long...
            segment_ids = getClusterSegmentsAll(sample_probes_sub, cluster_id);
            segments = c(segments, segment_ids);
        } else {
            if (nrow(probes) <=2 || nrow(sample_probes_sub) <= 2) {
                #Don't attempt to segment small clusters.
                segments = c(segments,cluster_id)
            } else {
                if (method == "hclust") {
                    segment_ids  = getClusterSegmentsHClust(
                        sample_probes = sample_probes_sub,
                        cluster_id = cluster_id,
                        do.plot=do.plot,
                        cutoff=cutoff,
                        dist.method=dist.method
                    );
                } else if (method == "skater") {
                    segment_ids = getClusterSegmentsSkater(
                        sample_probes_sub = sample_probes_sub,
                        cluster_id = cluster_id,
                        dist_method = dist.method,
                        cutoff=cutoff
                    );
                }
                segments = c(segments, segment_ids);
            }
        }
    }
    #Return the final list of segments.

    res = list();
    res$overlap_clusters = overlap_cluster_df;
    res$segments = segments;

    return(res);
}


#' Find unique set of epitope regions across all sample calls.
#'
#' @param probe_calls_res Results from makeProbeCalls
#'
#' @return vector of epitope seqments
#' @export
#'
#' @examples
findEpitopeSegmentsUnique<-function(
        probe_calls_res
        ) {
    segments = c();
    #message("Epitopes: Find sample epitopes")

    probe_calls = probe_calls_res$sample_probes;
    protein_tiling = probe_calls_res$protein_tiling;

    for (col_idx in seq_len(ncol(probe_calls))) {
        #message("col_idx:", col_idx);
        probes = rownames(probe_calls)[probe_calls[,col_idx]]
        if (length(probes) > 0) {
            epitopes = findBlocksProbeT(probes, protein_tiling);
            segments = c(segments, rownames(epitopes));
        }
    }
    segments = unique(segments);
    return(segments);
}


getSampleEpitopes<-function(probe_calls, protein_tiling)
{
    ans = list();
    for (col_idx in seq_len(ncol(probe_calls))) {
        probes = rownames(probe_calls)[,col_idx];
        if (length(probes) > 0) {
            ans[[colnames(probe_calls)[col_idx]]] = probes;
        }
    }
    return(ans);
}




#' Find epitopes using the unique method
#'
#' This function makes epitope level calls using an FDR matrix.
#' For each epitope, we calculate the maximum FDR of the probes within that
#' epitope and use that as the "FDR" for the epitope.  Just means, what
#' would the probe FDR cutoff need to be in
#' order to include this epitope for this sample.
#'
#' @param probe_sample_padj adjusted p-value matrix on the probe-level
#' @param probe_cutoff cutoff to use when calling probes for epitope finding
#' @param one_hit_filter use one hit filter to filter out epitopes
#' @param probes rownames of the probe_sample_padj
#' @param proteins proteins of the probes
#' @param pos starting position of the probe in the associated protein
#' @param protein_tiling tiling of the probe
#'
#' @return list of results
#'    sample_epitopes - epitopes found per sample (column) as list
#'    sample_epitopes_fdr - matrix of estimated FDRs per sample for
#'    sample_epitopes_minfdr = minFDRs;
#'    all_epitopes - unique vector of epitope ids
#'    k_of_n_epitopes
#' @export
#'
#' @examples
getEpitopeCallsUnique<-function(
    probe_sample_padj,
    probe_cutoff = 0.05,
    one_hit_filter = TRUE,
    probes = rownames(probe_sample_padj),
    proteins = getProteinLabel(probes),
    pos = getProteinStart(probes),
    protein_tiling = getProteinTiling(probes)
) {
    message("getEpitopeCallsUnique - start");

    probe_calls  = probe_sample_padj < probe_cutoff;

    message("Sample epitopes");
    sample_epitopes = getSampleEpitopes(probe_calls, protein_tiling);
    all_epitopes = unique(unlist(sample_epitopes))

    if (one_hit_filter) {
        epitope_fdrs =
            calcEpitopePValuesMat(probe_sample_padj, all_epitopes, "maxFDR");
        calls = makeCalls(epitope_fdrs, probe_cutoff);
        k_of_n = calls$k_of_n;
        message("Epitopes: applying one-hit filter");
        #Find all epitopes that have 1 probe and only 1 sample call.
        #Remove those epitopes from the results
        k1_epitopes = rownames(k_of_n)[k_of_n$K <= 1]
        n1_idx = oneProbeEpitopes(all_epitopes)
        n1_epitopes = all_epitopes[n1_idx]
        oh_epitopes = intersect(k1_epitopes,n1_epitopes);
        message("Filter ", length(oh_epitopes),
                " epitopes from ",length(k1_epitopes), " k1 epitopes");
        for (idx in seq_len(sample_epitopes)) {
            to_keep = !(sample_epitopes[[idx]]$Epitope_ID %in% oh_epitopes)
            sample_epitopes[[idx]] = sample_epitopes[[idx]][to_keep,];
        }
        to_keep = !(rownames(epitope_fdrs) %in% oh_epitopes)
        epitope_fdrs = epitope_fdrs[to_keep,]
        to_keep = !(rownames(minFDRs) %in% oh_epitopes)
        minFDRs = minFDRs[to_keep,];
        to_keep = !(all_epitopes$Epitope_ID %in% oh_epitopes)
        all_epitopes = all_epitopes[to_keep,]
        k_of_n = k_of_n[!(rownames(k_of_n) %in% oh_epitopes),]

    }
    message("Epitopes: Build list");
    ans = list(
        sample_epitopes = sample_epitopes,
        sample_epitopes_fdr = epitope_fdrs,
        sample_epitopes_minfdr = minFDRs,
        all_epitopes = all_epitopes,
        k_of_n_epitopes = k_of_n
    );
    return(ans);
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
#' findBlocksProbeT(c("A;1","A;2","A;3", "B;2","B;3", "C;10", "A;5","A;6"))
findBlocksProbeT<-function(
        probes,
        protein_tiling,
        proteins = getProteinLabel(probes),
        starts = getProteinStart(probes)
) {
    if (missing(protein_tiling)) {
        protein_tiling = 1;
    }
    protein.df = data.frame(
        Protein = proteins,
        Pos = starts,
        stringsAsFactors=FALSE
    );

    protein_list = split(protein.df, protein.df$Protein)

    ans_l = lapply(protein_list, findBlocksT, protein_tiling = protein_tiling);
    ans_dt = data.table::rbindlist(ans_l);
    ans_df = as.data.frame(ans_dt, stringsAsFactors=FALSE);

    ans_df$ProbeSet.AA.Span = ans_df$Stop - ans_df$Start + 1;
    rownames(ans_df) = getEpitopeID(ans_df$Protein, ans_df$Start, ans_df$Stop);
    return(ans_df);

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
    );
}

getTiling<-function(protein_tiling, protein) {
    if (missing(protein_tiling) || length(protein_tiling) == 0) {
        tiling = 1;
    } else if (length(protein_tiling) == 1) {
        tiling = protein_tiling[1];
    } else {
        tiling = protein_tiling[protein];
    }
    return(tiling);
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
        return(NULL);
    }
    if (nrow(prot_df)  == 1) {
        return(block_df(prot_df$Protein, prot_df$Pos, prot_df$Pos, 1))
    }
    tiling = getTiling(protein_tiling, prot_df$Protein[1]);
    prot_df = prot_df[order(prot_df$Pos, decreasing=FALSE),];
    ans_list = list();
    start_idx=1;
    start_pos = prot_df$Pos[start_idx];
    for (idx in 2:nrow(prot_df)) {
        current_idx = idx;
        current_pos = prot_df$Pos[current_idx];
        prev_idx = idx-1;
        prev_pos = prot_df$Pos[prev_idx];
        d = current_pos - prev_pos;
        if (d > tiling) {
            nprobes = length(seq(from=start_pos, to=prev_pos, by=tiling));
            ans_list[[length(ans_list)+1]] =
                block_df(prot_df$Protein[1], start_pos, prev_pos, nprobes);
            start_idx = current_idx;
            start_pos = current_pos;
        }
    }
    if (start_idx != nrow(prot_df)) {
        end_pos = prot_df$Pos[nrow(prot_df)]
        nprobes = length(seq(from=start_pos, to=end_pos, by=tiling))
        ans_list[[length(ans_list)+1]] =
            block_df(prot_df$Protein[1], start_pos, end_pos, nprobes)
    } else { #Last probe is by itself
        ans_list[[length(ans_list)+1]] =
            block_df(
                prot_df$Protein[1],
                prot_df$Pos[nrow(prot_df)],
                prot_df$Pos[nrow(prot_df)], 1
            )
    }
    ans.dt = data.table::rbindlist(ans_list);
    ans.df = as.data.frame(ans.dt);
    return(ans.df);
}


getHClustSilouette<-function(dist_mat2, hc) {
    sil_df = NULL;
    for (k in 2:max(2,(ncol(dist_mat2)-1))) {
        hc_cut = stats::cutree(hc, k = k)

        silhouette.results = cluster::silhouette(
            list(clustering=hc_cut),
            stats::as.dist(dist_mat2)
        );
        sil.sum = summary(silhouette.results);
        sil.mean = sil.sum$si.summary["Mean"]
        sil_df = rbind(
            sil_df,
            data.frame(
                K = k,
                SIL = as.vector(sil.mean),
                stringsAsFactors=FALSE
            ))

    }
    return(sil_df);
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
#' @export
#'
#' @examples
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
    sil_df = NULL;
    #Find the number of clusters using the silhouette
    if (cutoff == "silhouette") {
        sil_df = getHClustSilouette(dist_mat2, hc);
        K = max(sil_df$K[which(sil_df$SIL == max(sil_df$SIL))])
        #Use K that achieve the maximum mean silo score,
        #breaking ties using maximum K (bias toward more clusters).
        hc_cut = stats::cutree(hc, k= K);
    } else {
        hc_cut = stats::cutree(hc, h = cutoff);
    }

    if (is.unsorted(hc_cut)) {
        stop("clusters are not ordered correctly")
    }
    #Build the segments from the clusters.
    nsegments = max(hc_cut);

    cluster_protein = getEpitopeProtein(cluster_id);
    cluster_start = getEpitopeStart(cluster_id);
    cluster_stop  = getEpitopeStop(cluster_id);
    cluster_pos = getProteinStart(rownames(sample_probes));

    segment_ids = c();

    for (segment in seq_len(nsegments)) {
        segment_indices = which(hc_cut == segment);
        segment_pos = cluster_pos[segment_indices];
        segment_start = min(segment_pos);
        segment_end   = max(segment_pos);

        segment_id = getEpitopeID(cluster_protein, segment_start, segment_end);
        segment_ids = c(segment_ids, segment_id);
    }
    attr(segment_ids, "sil_df") = sil_df;
    attr(segment_ids, "hc") = hc;
    attr(segment_ids, "dist_mat") = attr(dist_mat2, "dist_mat");
    attr(segment_ids, "dist_mat2") = dist_mat2;
    return(segment_ids);

}

#' Calculate a distance matrix for use with hclust
#'
#' @param sample_probes_sub TODO
#' @param dist.method TODO
#' @param debug TODO
#'
#' @return matrix that can be converted to a distance for use with complete
#' hierarchical clustering
#' @export
#'
#' @examples
getHClustDistMat<-function(sample_probes_sub, dist.method, debug=FALSE) {

    s0 = Sys.time()

    maxl = nrow(sample_probes_sub);
    max_end = nrow(sample_probes_sub);
    if (debug) {message("Step 1");}
    #Build build 1st level distance
    dist_mat = NULL;

    if (dist.method == "hamming"){
        dist_mat = hamming(sample_probes_sub) / ncol(sample_probes_sub);
    } else {
        p = 2;
        if (dist.method == "minkowski") {p=ncol(sample_probes_sub);}
        dist_mat = as.matrix(
            stats::dist(
                sample_probes_sub,
                method=dist.method,
                diag=TRUE, upper=TRUE,
                p = p
            )
        )
    }
    if (debug){
        s1 = Sys.time()
        message("Step1:",s1-s0);
        message("Step 2");
    }
    #Enforce consecutive probes to be clustered together
    dist_mat2 = 1-diag(1,max_end, max_end);
    for (idx1 in seq_len((max_end-1))) {
        dist_mat2[idx1,idx1] = 0;
        for (idx2 in (idx1+1):max_end) {
            slength = idx2 - idx1 + 1;
            if (slength > maxl) {
                message("Max reached:", slength, ">", maxl);
                current_dist = 1;
            } else {
                ut = dist_mat[idx1:idx2,idx1:idx2]
                ut = ut[upper.tri(ut)]
                current_dist = max(unlist(ut));
            }
            dist_mat2[idx1,idx2] = current_dist;
            dist_mat2[idx2,idx1] = current_dist;
        }
    }
    if (debug) {
        s2 = Sys.time();
        message("Step 2:",s2-s1);
        message("Total:",s2-s0);
    }


    attr(dist_mat2, "dist_mat") = dist_mat;
    return(dist_mat2);
}



#' Get all possible cluster segments
#' determines all possible segmentation regions from a cluster/epitopeid
#' Uses sample_probes_sub to find the appropriate tiling for the cluster to use.
#' @param sample_probes_sub TODO
#' @param cluster_id TODO
#'
#' @return vector of epitope identifiers
#' @export
#'
#' @examples
getClusterSegmentsAll<-function(
        sample_probes_sub,
        cluster_id
) {

    cluster_protein = getEpitopeProtein(cluster_id)

    probes = getEpitopeProbeIDs(cluster_id)
    if (length(probes) > 1) {
        probes = probes[probes %in% rownames(sample_probes_sub)]
    }
    cluster_pos = getProteinStart(probes)
    if (sum(is.na(cluster_pos)) > 0) {
        print(cluster_id)
        print(probes)
        print(cluster_pos)
        stop("NA1");
    }
    n = length(cluster_pos);

    segments = c();
    for (idx1 in seq_len(n)) {
        pos1 = cluster_pos[idx1]
        for (idx2 in idx1:n) {
            segment = getEpitopeID(cluster_protein, pos1, cluster_pos[idx2])
            segments = c(segments, segment);
        }
    }

    starts = getEpitopeStart(segments)
    if (sum(is.na(starts)) > 0) {
        print(cluster_id)
        print(probes)
        print(cluster_pos)
        print(segments)
        stop("NA2");

    }


    return(segments);
}


getSkaterGraph<-function(n) {
    edges = data.frame(
        A = seq_len(n-1),
        B = seq_len(n-1)+1
    )
    edges = as.matrix(edges);
    return(edges);
}

getSkaterDist<-function(sample_probes_sub_i, dist_method, p) {
    if (dist_method == "hamming") {
        dist_method = function(data, id) {
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
        sk_dist = hamming_dist(sample_probes_sub_i)
    } else {
        sk_dist = stats::dist(sample_probes_sub_i, method=dist_method, p = p)
    }
    return(sk_dist);
}

getSkaterSilouette<-function(edges, s_p_sub_i, dist_method, p, sk_dist) {
    sil_df = NULL;
    sk_res = NULL;
    n = nrow(edges)
    for (ncuts in seq_len(n-1)) {
        if (is.null(sk_res)) {
            sk_res = spdep::skater(edges, s_p_sub_i, ncuts=1,
                method = dist_method, p = p)
        } else {
            sk_res = spdep::skater(sk_res, s_p_sub_i, ncuts=1,
                method = dist_method, p = p)
        }
        sil_res = cluster::silhouette(list(clustering=sk_res$groups), sk_dist)
        sil_df = rbind(
            sil_df,
            data.frame(
                ncuts = ncuts,
                K = ncuts + 1,
                SIL = as.vector(summary(sil_res)$si.summary["Mean"]),
            )
        )
    }
    NCUTS = max(sil_df$ncuts[sil_df$SIL == max(sil_df$SIL)]);
    sk_res = spdep::skater(
        edges,
        s_p_sub_i,
        ncuts = NCUTS, method = dist_method, p=p
    )
    attr(sk_res, "sil_df") = sil_df;
    return(sk_res);
}



#' Title
#'
#' Acceptable dist.methods are "euclidean", "hamming"
#'    "maximum", "manhattan", "canberra", "binary", "minkowski"
#'  Could also pass in the function, but I haven't tested that yet.
#' Also "mahalanobis" requires a covariance matrix,
#'  which I have implemented/tried yet...
#' Minkowski is hardcoded to have a power (p) of 0.5,
#'  maybe we should give the user
#' an option to modify that?
#' @param sample_probes_sub TODO
#' @param cluster_id TODO
#' @param dist_method TODO
#' @param cutoff TODO
#' @param debug TODO
#'
#' @return vector of segments
#' @export
#'
#' @examples
getClusterSegmentsSkater<-function(
        sample_probes_sub,
        cluster_id,
        dist_method = "euclidean",
        cutoff = "silhouette",
        debug = FALSE) {

    n = nrow(sample_probes_sub);
    if (n == 1) {
        stop("Need more than 1 point to cluster!")
    }
    edges = getSkaterGraph(n)
    s_p_sub_i = toNumericMatrix(sample_probes_sub);
    p = 2;
    if (dist_method == "minkowski") {
        p = ncol(s_p_sub_i)
    }
    sk_dist = getSkaterDist(s_p_sub_i, dist_method, p);
    if (cutoff == "silhouette") {
        sk_res = getSkaterSilouette(edges, s_p_sub_i, dist_method, p, sk_dist);
    } else {
        stop("cutoff method Not Implemented: ", cutoff );
    }

    nsegments = max(sk_res$groups);

    cluster_protein = getEpitopeProtein(cluster_id);
    cluster_start = getEpitopeStart(cluster_id);
    cluster_stop  = getEpitopeStop(cluster_id);
    cluster_pos = getProteinStart(rownames(sample_probes_sub));

    segment_ids = c();

    for (segment in seq_len(nsegments)) {
        segment_indices = which(sk_res$groups == segment);
        segment_pos = cluster_pos[segment_indices];
        segment_start = min(segment_pos);
        segment_end = max(segment_pos);
        segment_id  = getEpitopeID(cluster_protein, segment_start, segment_end);
        segment_ids = c(segment_ids, segment_id);
    }
    attr(segment_ids, "sk_res") = sk_res;
    attr(segment_ids, "dist") = sk_dist;
    return(segment_ids);
}




