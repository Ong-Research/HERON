#' SARS CoV-2 Wuhan Peptide Binding Array Data
#'
#' A subset of data from the paper
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8245122/ publication.
#'
#' @format ## `heffron2020_wuhan`
#' A data frame with 1945 rows and 60 columns. Each column is a pre-processed
#' binding signal from a serum sample peptide array set for the SARS-CoV-2.
#' The matrix is a subset of the full matrix and contains sequences from the
#' membrane, envelope, surface (spike), and nucleocapsid proteins.
#'
#' The "probe_meta" attribute is a data frame with 1945 rows and 6 columns.
#' The columns are POSITION - starting position of probe within protein,
#' PROBE_SEQUENCE - amino acid sequence of probe, SEQ_ID - protein identifier
#' SEQ_NAME - name of protein, PROBE_ID - combination of protein identifier
#' and starting position, e.g. prot1;5.
#'
#' The "sample_meta" attribute is a data frame with 60 rows and 2 columns.
#' The columns are SAMPLE_NAME - name of the sample which is a column name in
#' sample matrix and COVID_POSITIVE - whether the sample is positive for
#' COVID.
#'
#' @source <https://github.com/Ong-Research/UW_Adult_Covid-19>
"heffron2020_wuhan"
