

#' Making Probe-level Calls
#'
#' \code{makeProbeCalls} returns call information on a probe matrix that has
#' been scored by the
#' @param probe_sample_padj a matrix of adjusted p-values where each row is a
#' feature and column is a sample
#' @param probe_cutoff
#' @param probes
#' @param proteins
#' @param positions
#' @param protein_tiling
#' @param one_hit_filter
#'
#' @return
#' @export
#'
#' @examples
makeProbeCalls<-function(probe_sample_padj,
                         probe_cutoff = 0.05,
                         probes = rownames(probe_sample_padj),
                         proteins = getProteinLabel(probes),
                         positions = getProteinStart(probes),
                         protein_tiling = getProteinTiling(probes),
                         one_hit_filter = TRUE
) {

    sample_probes = probe_sample_padj < probe_cutoff; # make calls.
    minFDRs = calcMinFDR(as.matrix(probe_sample_padj), additional_stats = FALSE);

    k_of_n = minFDRs;
    colnames(k_of_n) = paste0("K",1:ncol(minFDRs),".padj");
    k_of_n = cbind(rowSums(minFDRs < probe_cutoff), k_of_n);
    colnames(k_of_n)[1] = "K";

    rownames(k_of_n) = rownames(minFDRs);

    #k_of_n = k_of_n[k_of_n$K != 0,]

    if (one_hit_filter) {
        message("Hit support");
        probe_hit_support = probeHitSupported(sample_probes, probes, proteins, positions, protein_tiling);
        message("Applying one-hit filter");
        k1_probes = rownames(k_of_n)[k_of_n$K == 1]; #Probe appears in only 1 sample.

        supported = rowSums(sample_probes & probe_hit_support) > 0;
        to_remove = intersect(k1_probes,names(supported)[!supported]); #Remove any k1 probe without consecutive probe support.

        message("removing ",length(to_remove), " out of ", length(k1_probes)," k1 probes");
        #Instead of removing, just set the padj values to 1.0.
        sample_probes[rownames(sample_probes) %in% to_remove,] = FALSE;
        k_of_n$K[rownames(k_of_n) %in% to_remove]=0;
        k_of_n[rownames(k_of_n) %in% to_remove, 2:ncol(k_of_n)] = 1.0;

        one_hit_padj = probe_sample_padj;
        one_hit_padj[rownames(one_hit_padj) %in% to_remove,] = 1.0;

    }

    #Make a list of all of the results
    ans = list();
    ans$sample_probes = sample_probes;
    ans$minFDRs = minFDRs;
    ans$k_of_n_probes = k_of_n;
    ans$probe_sample_padj = probe_sample_padj;
    if (one_hit_filter) {
        ans$k1_probes = k1_probes;
        ans$probe_hit_support = probe_hit_support;
        ans$supported = names(supported)[supported];
        ans$to_remove = to_remove;
        ans$to_keep = setdiff(rownames(probe_hit_support), to_remove);
        ans$probe_sample_padj = one_hit_padj;
        ans$one_hit_padj = one_hit_padj;
        ans$orig_padj = probe_sample_padj;
    }

    #Save parameters.
    ans$probe_cutoff = probe_cutoff;
    ans$one_hit_filter = one_hit_filter;
    ans$probes = probes;
    ans$proteins = proteins;
    ans$positions = positions;
    ans$protein_tiling = protein_tiling;
    return(ans);
}
