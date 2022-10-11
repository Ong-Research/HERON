


#' Obtain Protein Id from Epitope ID
#'
#' Format of EpitopeID is A_B_C, where A is the protein label
#' B is the protein start position of the first probe in the epitope and
#' C is the protein start position of the last probe in the epitope.
#'
#' @param epitope_ids vector of epitope identifier character strings
#'
#' @return vector of protein labels
#' @export
#'
#' @examples getEpitopeProtein("Prot1_1_5")
getEpitopeProtein<-function(epitope_ids) {
    ans = unlist(
        lapply(
            strsplit(epitope_ids,"_"),
            function(l) {
                ans = l[seq_len((length(l)-2))]
                return(paste(ans,collapse="_"))
            }
        )
    )
    return(ans);
}

#' Obtain first probe's protein start position from Epitope ID
#'
#'
#' @param epitope_ids vector of epitope ids
#'
#' @return vector of integers indicating first probe start positions in the
#' epitope(s)
#'
#' @export
#'
#' @examples getEpitopeStart("Prot1_1_5")
getEpitopeStart<-function(epitope_ids) {
    ans = unlist(
        lapply(
            strsplit(epitope_ids,"_"),
            function(l) {
                return(l[length(l)-1])
            }
        )
    )
    return(as.integer(ans));
}

#' Obtain last probe's protein start position from EpitopeID
#'
#' @param epitope_ids vector of epitope ids
#'
#' @return vector of integers indicating the last probe protein start position
#' @export
#'
#' @examples getEpitopeStop("Prot1_1_5")
getEpitopeStop<-function(epitope_ids) {
    ans = unlist(
        lapply(
            strsplit(epitope_ids,"_"),
            function(l) {
                return(l[length(l)])
            }
        )
    )
    return(as.integer(ans));
}


#' Indicate which epitopes are just one probe.
#'
#' @param epitope_ids vector of epitope ids
#'
#' @return vector of logical indicating epitopes that are one probe
#' @export
#'
#' @examples oneProbeEpitopes("A_1_1", "B_1_1","C_1_2")
oneProbeEpitopes<-function(epitope_ids) {
    return(getEpitopeStart(epitope_ids) == getEpitopeStop((epitope_ids)));
}


#' Get the vector of probes from an epitope id
#'
#' @param epitope_id EpitopeID to obtain probes from
#' @param tiling Tiling of the probes across the protein (default 1)
#'
#' @return vector of probe_ids that are contained within the epitope
#' @export
#'
#' @examples
#' getEpitopeProbeIDs("A_1_5")
getEpitopeProbeIDs<-function(epitope_id, tiling=1) {
    protein = getEpitopeProtein(epitope_id);
    start = getEpitopeStart(epitope_id)
    stop = getEpitopeStop(epitope_id);

    ans = paste0(protein,";",seq(from=start, to=stop, by=tiling));
    return(ans);
}


#' Get probe ids from a vector of epitope ids
#'
#' @param epitope_ids vector of epitope identifiers
#' @param tiling tling of probes across proteins
#'
#' @return data.frame of epitope_to_probe mappings
#' @export
#'
#' @examples
#' getEpitopeIDsToProbeIDs(c("A_1_5","C_8_12"))
getEpitopeIDsToProbeIDs<-function(epitope_ids, tiling=1) {
    epitope_to_probe = NULL;

    epitope_to_probe_list = list();

    for (epitope_idx in seq_len(length(epitope_ids))) {
        epitope_id = epitope_ids[epitope_idx];
        if (length(tiling) == length(epitope_ids)) {
            #message("Using custom tiling");
            epi_probes = getEpitopeProbeIDs(epitope_id, tiling[epitope_idx]);
        } else {
            epi_probes = getEpitopeProbeIDs(epitope_id, tiling);
        }
        epitope_to_probe_list[[epitope_idx]] =
            data.frame(
                Epitope_ID = rep(epitope_id, length(epi_probes)),
                PROBE_ID = epi_probes
            );
    }
    epitope_to_probe = data.table::rbindlist(epitope_to_probe_list);
    epitope_to_probe = as.data.frame(epitope_to_probe, stringsAsFactors=FALSE);

    return(epitope_to_probe);
}



#' Create EpitopeID from protein, first and last probes
#'
#' @param protein vector of proteins
#' @param start vector of first probe protein start positions
#' @param stop vector of last probe protein start positions
#'
#' @return vector of epitope ids
#' @export
#'
#' @examples
#' getEpitopeID("A", 1, 2)
getEpitopeID<-function(protein, start, stop) {
    return(paste(protein,start,stop, sep="_"))
}





#' Obtain Sequence Annotations for Epitopes
#'
#' @param epitopes vector of epitope ids
#' @param probe_meta data.frame with the PROBE_ID, SEQ_ID, and PROBE_SEQUENCE
#' columns defined
#' @param debug output debugging information
#'
#' @return data.frame in same order of the epitopes parameter with sequence
#' annotations
#'
#' @export
#'
#' @examples
getSequenceAnnotations<-function(epitopes, probe_meta, debug = FALSE) {
    eproteins = getEpitopeProtein(epitopes)
    estarts = getEpitopeStart(epitopes)
    estops = getEpitopeStop(epitopes)

    umeta = unique(probe_meta[probe_meta$SEQ_ID %in% eproteins,
        c("PROBE_ID", "PROBE_SEQUENCE")])
    umeta$SEQUENCE_LENGTH = nchar(umeta$PROBE_SEQUENCE)
    rownames(umeta) = umeta$PROBE_ID;


    first_probe = paste0(eproteins, ";", estarts )
    last_probe = paste0(eproteins, ";", estops)

    first_length = umeta[first_probe, "SEQUENCE_LENGTH"]
    first_last_pos = estarts + first_length - 1;

    ans_df = data.frame(
        EpitopeID = epitopes,
        Overlap.Sequence.Length = rep(NA, length(epitopes)),
        Full.Sequence.Start = estarts,
        Full.Sequence.Stop = rep(NA, length(epitopes)),
        Full.Sequence.Length = rep(NA, length(epitopes)),
        First.Sequence = umeta[first_probe, "PROBE_SEQUENCE"],
        Last.Sequence = umeta[last_probe, "PROBE_SEQUENCE"],
        Overlap.Sequence = rep("", length(epitopes)),
        Full.Sequence = rep(NA, length(epitopes)),
        stringsAsFactors = FALSE
    );

    for (idx in seq_len(nrow(ans_df))) {
        start = estops[idx]
        stop = first_last_pos[idx]
        if (start <= stop) {
            #overlap sequence is defined, subset the string.
            pstart = start - estarts[idx] + 1
            pstop = stop - estarts[idx] + 1
            ans_df$Overlap.Sequence[idx] = substr(ans_df$First.Sequence[idx],
                pstart, pstop)
        }
        if (estarts[idx] == estops[idx]) {
            # For an epitope of length 1, full sequence is the first sequence
            ans_df$Full.Sequence[idx] = ans_df$First.Sequence[idx]
        }
        else {
            if (start <= stop) {
                #If first and last sequence overlap, then just concatenate.
                #Else, stich the whole sequence together.
                ans_df$Full.Sequence[idx] =
                    catSequences(
                        c(estarts[idx],
                        estops[idx]),
                        c(ans_df$First.Sequence[idx],
                        ans_df$Last.Sequence[idx])
                    )
            }
            else {
                #stitch together full sequence using all of the probes
                probes = paste0(eproteins[idx], ";", estarts[idx]:estops[idx])
                probes = probes[probes %in% rownames(umeta)]
                ans_df$Full.Sequence[idx] =
                    catSequences(
                        getProteinStart(probes),
                        umeta[probes, "PROBE_SEQUENCE"]
                    )
            }
        }
        if (debug && idx%%500 == 0) {
            #Print out every 500 iterations.
            cat(idx, " of ", nrow(ans_df), "\n")
        }
    }

    ans_df$Overlap.Sequence.Length = nchar(ans_df$Overlap.Sequence)
    ans_df$Full.Sequence.Length = nchar(ans_df$Full.Sequence)
    ans_df$Full.Sequence.Stop = ans_df$Full.Sequence.Start +
        ans_df$Full.Sequence.Length - 1
    #Return without the epitope identifer column
    return(ans_df[,-1]);
}


#' Extract aggregate metric statistics from blocks (epitopes)
#'
#' For each block a
#' @param blocks data.frame of with Protein, Start, and Stop columns
#' @param probes list of all probes in an experiment
#' @param metric numeric.vector of values to collect stats on
#' @param label metric label to use in returned values
#' @param proteins getProteinLabel(probes) - cached results
#' @param pos getProteinStart(probes) - cached results
#'
#' @return data.frame with the min, max, and mean value for each block
#' extracted from the metric
#' @export
#'
#' @examples
getBlockMetricStats<-function(blocks, probes, metric, label,
    proteins=getProteinLabel(probes),
    pos = getProteinStart(probes)) {

    keep = proteins %in% blocks$Protein;

    metric.df = data.frame(
        Protein = proteins[keep],
        Pos = pos[keep],
        Metric = metric[keep],
        stringsAsFactors=FALSE
    );

    metric_list = split(metric.df, metric.df$Protein);

    #cat("getMetricStats - start\n");
    min_value = rep(NA, nrow(blocks))
    max_value = rep(NA, nrow(blocks))
    mean_value = rep(NA, nrow(blocks))



    for (idx in seq_len(nrow(blocks))) {
        temp.df = metric_list[[blocks$Protein[idx]]];
        idx_idx = temp.df$Pos >= blocks$Start[idx] &
            temp.df$Pos <= blocks$Stop[idx]
        metric_values = temp.df$Metric[idx_idx];
        min_value[idx] = min(metric_values, na.rm=TRUE);
        max_value[idx] = max(metric_values, na.rm=TRUE);
        mean_value[idx] = mean(metric_values, na.rm=TRUE);
    }

    ans = data.frame(
        min_value = min_value,
        max_value = max_value,
        mean_value = mean_value
    );
    colnames(ans) = c(
        paste0("Min.",label),
        paste0("Max.",label),
        paste0("Mean.",label)
    )
    return(ans);
}

