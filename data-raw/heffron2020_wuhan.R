## code to prepare `heffron2020_wuhan` dataset goes here

##Probe Meta
sequences_url = "https://dholk.primate.wisc.edu/_webdav/dho/sequencing/Polypeptide%20Microarrays/public/COVID_19/%40files/all_sequences_except_wi.tsv.gz?contentDisposition=attachment"
sequences_path = "all_sequences_except_wi.tsv.gz"

if (!file.exists(sequences_path)) {
    download.file(sequences_url, sequences_path);
}
probe_meta = read.table(sequences_path, sep="\t", header = TRUE);
probe_meta$PROBE_ID = paste0(probe_meta$SEQ_ID, ";", probe_meta$POSITION);

probe_meta_wu1 = probe_meta[grep("Wu1", probe_meta$SEQ_ID),]
probe_meta_wu1 = probe_meta_wu1[grep("orf", probe_meta_wu1$SEQ_ID, invert=TRUE),]

##SeqMat Data
stacked_df_url = "https://dholk.primate.wisc.edu/_webdav/dho/sequencing/Polypeptide%20Microarrays/public/COVID_19/%40files/aggregated_data/df_stacked.tsv.gz?contentDisposition=attachment"
stacked_df_path = "df_stacked.tsv.gz"

if (!file.exists(stacked_df_path)) {
    options(timeout = max(300, getOption("timeout")))

    download.file(stacked_df_url, stacked_df_path);
}

message("Loading matrix")
seq_mat = UW.Adult.Covid.19::loadSeqMat(file_name = stacked_df_path);
sample_meta = attr(seq_mat, "sample_meta")
seq_mat = seq_mat[,-1]

heffron2020_wuhan = seq_mat[rownames(seq_mat) %in% probe_meta_wu1$PROBE_SEQUENCE,];

attr(heffron2020_wuhan, "sample_meta") = sample_meta;
attr(heffron2020_wuhan, "probe_meta") = probe_meta_wu1;

usethis::use_data(heffron2020_wuhan, overwrite = TRUE)
