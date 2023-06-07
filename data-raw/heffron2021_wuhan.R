## code to prepare `heffron2021_wuhan` dataset goes here

##Probe Meta
sequences_url = "https://dholk.primate.wisc.edu/_webdav/dho/sequencing/Polypeptide%20Microarrays/public/COVID_19/%40files/all_sequences_except_wi.tsv.gz?contentDisposition=attachment"
sequences_path = "all_sequences_except_wi.tsv.gz"

if (!file.exists(sequences_path)) {
    download.file(sequences_url, sequences_path);
}
probe_meta <- read.table(sequences_path, sep="\t", header = TRUE);
probe_meta$PROBE_ID <- paste0(probe_meta$SEQ_ID, ";", probe_meta$POSITION);

## Reduce the matrix to make it small enough for examples.
probe_meta_wu1 <- probe_meta[grep("Wu1", probe_meta$SEQ_ID),]
probe_meta_wu1 <-
    probe_meta_wu1[grep("orf", probe_meta_wu1$SEQ_ID, invert=TRUE),]

##SeqMat Data
stacked_df_url <- "https://dholk.primate.wisc.edu/_webdav/dho/sequencing/Polypeptide%20Microarrays/public/COVID_19/%40files/aggregated_data/df_stacked.tsv.gz?contentDisposition=attachment"
stacked_df_path <- "df_stacked.tsv.gz"

if (!file.exists(stacked_df_path)) {
    options(timeout = max(300, getOption("timeout")))
    download.file(stacked_df_url, stacked_df_path);
}

message("Loading matrix")

if (!require(UW.Adult.Covid.19)) {
    devtools::install_github("Ong-Research/UW_Adult_Covid-19/UW.Adult.Covid.19")
}

seq_mat <- UW.Adult.Covid.19::loadSeqMat(file_name = stacked_df_path);
sample_meta = attr(seq_mat, "sample_meta")
seq_mat = seq_mat[,-1]

## Create colData data.frame
create_colData<-function(mat_in) {
    colData <- data.frame(
        Sample_ID = colnames(mat_in),
        ptid = colnames(mat_in),
        visit = "pre",
        condition = "Control",
        stringsAsFactors=FALSE
    );
    pos_samples <- sample_meta$SAMPLE_NAME[sample_meta$COVID_POSITIVE == "YES"]
    colData$condition[colData$Sample_ID %in% pos_samples] = "COVID";
    colData$visit[colData$condition == "COVID"] = "post"
    colData$TAG <- colData$Sample_ID;
    rownames(colData) <- colData$TAG
    return(DataFrame(colData));
}

seq_mat <- seq_mat[rownames(seq_mat) %in% probe_meta_wu1$PROBE_SEQUENCE,]

colData_heffron <- create_colData(seq_mat)

heffron2021_wuhan <- HERON::HERONSequenceDataSet(exprs = seq_mat)
colData(heffron2021_wuhan) <- colData_heffron
metadata(heffron2021_wuhan)$probe_meta <- probe_meta_wu1
head(colData(heffron2021_wuhan))
head(assay(heffron2021_wuhan))
usethis::use_data(heffron2021_wuhan, overwrite = TRUE)

