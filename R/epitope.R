


#' Obtain Protein Id from Epitope ID
#'
#' Format of EpitopeID is A_B_C, where A is the protein label
#' B is the protein start position of the first probe in the epitope and
#' C is the protein start position of the last prboe in the epitope.
#'
#' @param epitope_ids vector of epitope_ids
#'
#' @return vector of protein labels
#' @export
#'
#' @examples
getEpitopeProtein<-function(epitope_ids) {
    ans = unlist(
        lapply(
            strsplit(epitope_ids,"_"),
            function(l) {
                ans = l[1:(length(l)-2)]
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
#' @examples
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
#' @examples
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

#' Get the vector of probes from an epitope id
#'
#' @param epitope_id EpitopeID to obtain probes from
#' @param tiling Tiling of the probes across the protein (default 1)
#'
#' @return vector of probe_ids that are contained within the epitope
#' @export
#'
#' @examples
getEpitopeProbeIDs<-function(epitope_id, tiling=1) {
    protein = getEpitopeProtein(epitope_id);
    start = getEpitopeStart(epitope_id)
    stop = getEpitopeStop(epitope_id);

    ans = paste0(protein,";",seq(from=start, to=stop, by=tiling));
    return(ans);
}

#' Obtain the epitope ids from the vectors of protein, first and last probes
#'
#' @param protein vector of proteins
#' @param start vector of first probe protein start positions
#' @param stop vector of last probe protein start positions
#'
#' @return vector of epitope ids
#' @export
#'
#' @examples
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

    umeta = unique(probe_meta[probe_meta$SEQ_ID %in% eproteins, c("PROBE_ID", "PROBE_SEQUENCE")])
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
        First.Sequence = umeta[first_probe, "PROBE_SEQUENCE"],
        Last.Sequence = umeta[first_probe, "PROBE_SEQUENCE"],
        Overlap.Sequence = rep("", length(epitopes)),
        Full.Sequence = rep(NA, length(epitopes)),
        stringsAsFactors = FALSE
    );

    for (idx in 1:nrow(ans_df)) {
        start = estops[idx]
        stop = first_last_pos[idx]
        if (start <= stop) {
            #overlap sequence is defined.
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
                ans_df$Full.Sequence[idx] =
                    catSequences(c(estarts[idx],
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
            cat(idx, " of ", nrow(ans_df), "\n")
        }
    }

    ans_df$Overlap.Sequence.Length = nchar(ans_df$Overlap.Sequence)
    ans_df$Full.Sequence.Length = nchar(ans_df$Full.Sequence)
    ans_df$Full.Sequence.Stop = ans_df$Full.Sequence.Start +
        ans_df$Full.Sequence.Length - 1

    return(ans_df[,-1]);
}


