makeMultilevelCalls<-function(
    probe_pvalues_res,
    probe_calls,
    segments,
  epitope.pvalue.method = "maxFDR",
  protein.pvalue.method = "minFDR",
  epitope_cutoff = probe_calls$probe_cutoff,
  protein_cutoff = probe_calls$probe_cutoff
) {
    probe_sample_pvalues = attr(probe_pvalues_res,"pvalue")
    probe_sample_padj = probe_pvalues_res

    if (epitope.pvalue.method == "maxFDR") {
        probe_sample_pvalues = probe_sample_padj;
    }

    epitope_sample_pvalues = calcEpitopePValues(
        epitope_ids = segments,
        probe_pvalues = probe_sample_pvalues,
        method = epitope.pvalue.method
    );

    if (epitope.pvalue.method == "maxFDR") {
        epitope_sample_padj = epitope_sample_pvalues
    } else {
        epitope_sample_padj = p_adjust_mat(epitope_sample_pvalues);
    }

    message("Epitopes: K of N");
    #Now estimate the K of N by first finding the min fdr per K.

    epitope_k_of_n = makeCalls(padj_mat = epitope_sample_padj,cutoff = epitope_cutoff)


    minFDRs = protStat_calcMinFDR(as.matrix(epitope_sample_padj), additional_stats=FALSE);

    k_of_n = minFDRs;
    colnames(k_of_n) = paste0("K",1:ncol(minFDRs),".padj");

    K = rowSums(minFDRs < epitope_cutoff)
    F = K / ncol(minFDRs)

    k_of_n = cbind(K,F,k_of_n);
    colnames(k_of_n)[1:2] = c("K","F"); #K is the count, F is the fraction.
    epitope_k_of_n = k_of_n;

    if (epitope.pvalue.method == "maxFDR" || protein.pvalue.method == "minFDR") {
        message("Proteins: minFDR");
        protein_sample_padj = getProteinPValuesE(epitope_sample_padj, method="min")
    } else {
        message("Proteins: pvalues");
        protein_sample_pvalues = getProteinPValuesE(epitope_sample_pvalues, method=protein.pvalue.method);
        protein_sample_padj = p.adjust.mat(protein_sample_pvalues);
    }
    message("Proteins: K of N");
    minFDRs = protStat_calcMinFDR(as.matrix(protein_sample_padj), additional_stats=FALSE);
    k_of_n = minFDRs;
    colnames(k_of_n) = paste0("K", 1:ncol(minFDRs),".padj");

    K = rowSums(minFDRs < protein_cutoff);
    F = K / ncol(minFDRs);
    k_of_n = cbind(K,F,k_of_n);
    colnames(k_of_n)[1:2] = c("K","F"); #K is the count, F is the fraction.
    protein_k_of_n = k_of_n;


    ans = list();
    #ans$probe_calls = probe_calls
    ans$epitope_calls = list();

    ans$epitope_calls$segments = segments;
    if (epitope.pvalue.method != "maxFDR") {
        ans$epitope_calls$epitope_sample_pvalues = epitope_sample_pvalues;
    }
    ans$epitope_calls$epitope_sample_padj = epitope_sample_padj;
    ans$epitope_calls$epitope_sample_calls = epitope_sample_padj < epitope_cutoff;
    ans$epitope_calls$epitope_k_of_n = epitope_k_of_n;
    ans$protein_calls = list();
    if (epitope.pvalue.method != "maxFDR" && protein.pvalue.method != "minFDR") {
        ans$protein_calls$protein_sample_pvalues = protein_sample_pvalues;
    }
    ans$protein_calls$protein_sample_padj = protein_sample_padj;
    ans$protein_calls$protein_sample_calls = protein_sample_padj < protein_cutoff;
    ans$protein_calls$protein_k_of_n = protein_k_of_n;
    return(ans);

}

#' Making Probe-level Calls
#'
#' \code{makeProbeCalls} returns call information on a probe matrix that has
#' been scored by the
#' @param probe_sample_padj a matrix of adjusted p-values where each row is a
#' feature and column is a sample
#' @param probe_cutoff cutoff to use when calling probes
#' @param one_hit_filter Indicator to remove probe hits that do not have a
#' matching consecutive probe and if the probe is only found in 1 sample
#' @param probes default rownames
#' @param proteins default getProteinLabel(probes), used for speed up for
#' multiple calls
#' @param positions default getProteinStart(probes), used for speed up for
#' multiple calls
#' @param protein_tiling default getProteinTiling(probes), used for speed up
#' for multiple calls
#'
#' @return list of results
#' sample_probes -> logical matrix of calls on probes
#'
#' @export
#'
#' @examples
makeProbeCalls<-function(probe_sample_padj,
                         probe_cutoff = 0.05,
                         one_hit_filter = TRUE,
                         probes = rownames(probe_sample_padj),
                         proteins = getProteinLabel(probes),
                         positions = getProteinStart(probes),
                         protein_tiling = getProteinTiling(probes)
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



#' Find probe hits that are supported by a consecutive probe or another sample
#'
#' @param hit_mat matrix of logical values that indicate a hit with a TRUE value
#' @param probes default rownames, used for speed for across calls that need
#' this information
#' @param proteins default getProteinLabel(probes), used for speed up for
#' multiple getting the protein information from the probes
#' @param positions default getProteinStart(probes), used for speed up for
#' multiple calls
#' @param tiling default getProteinTiling(probes), used for speed up
#' for multiple calls
#'
#' @return matrix of logical values indicate that the TRUE hit is supported by
#' a consecutive probe hit in the sample sample or the within another sample
#' @export
#'
#' @examples
probeHitSupported<-function(hit_mat,
                            probes=rownames(hit_mat),
                            proteins = getProteinLabel(probes),
                            positions = getProteinStart(probes),
                            tiling = getProteinTiling(probes)) {
    hit_df = as.data.frame(hit_mat,stringsAsFactors=FALSE);
    hit_df$Protein = proteins;
    hit_df$Pos = positions;
    hit_df$Order = 1:nrow(hit_df);
    hit_df_protein = split(hit_df, hit_df$Protein);
    cols = 1:(ncol(hit_df)-3)
    hit_df_protein = lapply(
        hit_df_protein,
        function(current_df) {
            current_tile = tiling[current_df$Protein[1]];
            ans = current_df;
            #ans = ans[order(ans$Pos, decreasing=FALSE),]
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
            #if (sum(NAs) > 0) {
            ansl[NAs] = FALSE;
            #}
            ansr = ans;
            ansr[pos_df$pos.label,cols] =
                ans[pos_df$pos.label,cols] & ans[pos_df$posr.label,cols];
            NAs = is.na(ansr);
            #if (sum(NAs) > 0) {
            ansr[NAs] = FALSE;
            #}
            ans.or = ans;
            ans.or[pos_df$pos.label,cols] = ansl[pos_df$pos.label,cols] | ansr[pos_df$pos.label,cols];
            rownames(ans.or) = pos_df[rownames(ans.or),"orig"];
            return(ans.or);
        }
    )

    ans.dt = data.table::rbindlist(hit_df_protein);
    ans.df = as.data.frame(ans.dt, stringsAsFactors=FALSE);
    ans.df = ans.df[ans.df$Order,];
    rownames(ans.df) = paste0(ans.df$Protein,";",ans.df$Pos); #I don't know why this is lost...
    ans.df = ans.df[,cols];
    return(ans.df)

}


makeCalls<-function(padj_mat, cutoff = 0.05) {
    calls = padj_mat < cutoff;
    minFDRs = calcMinFDR(as.matrix(padj_mat), additional_stats = FALSE);
    k_of_n = minFDRs
    colnames(k_of_n) = paste0("K",1:ncol(minFDRs),".padj");
    K = rowSums(minFDRs < epitope_cutoff)
    F = K / ncol(minFDRs)
    K.padj = rep(1, nrow(padj_mat))
    for (idx in 1:nrow(padj_mat)) {
        if (K[idx] > 0) {
          K.padj[idx] = k_of_n[idx,K[idx]]
        }
    }
    ans = cbind(K, F, K.padj, k_of_n)
    colnames(ans)[1:3] = c("K", "F", "K.padj")
    return(ans);
}


makeEpitopeCalls<-function(epitope_sample_padj,
                           epitope_cutoff = 0.05
) {
    ans = makeCalls(epitope_sample_padj, epitope_cutoff);
    return(ans);
}


