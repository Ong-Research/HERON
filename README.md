
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HERON

<!-- badges: start -->
<!-- badges: end -->

The goal of HERON (**H**ierarchical **E**pitope P**R**Otein Bi**N**ding)
is to anayzing peptide binding array data.

## Installation

You can install the released version of HERON from
[github](https://github.com/Ong-Research/HERON) with:

``` r
devtools::install_github("Ong-Research/HERON")
```

And the development version from
[GitHub](https://github.com/Ong-Research/HERON) with:

``` r
# install.packages("devtools")
devtools::install_github("Ong-Research/HERON")
```

# Flowchart of probe p-value calculations

<img src="man/figures/README-flowchart_probe_pvalues-1.png" width="100%" />

## Example

These are examples which shows you how to interact with the code. We
will be using the COVID-19 peptide binding array dataset from <TODO>
this publication.

### Download Data

``` r

##Probe Meta
sequences_url = "https://dholk.primate.wisc.edu/_webdav/dho/sequencing/Polypeptide%20Microarrays/public/COVID_19/%40files/all_sequences_except_wi.tsv.gz?contentDisposition=attachment"
sequences_path = "all_sequences_except_wi.tsv.gz"
download.file(sequences_url, sequences_path);
probe_meta = UW.Adult.Covid.19::loadProbeMeta(sequences_path)

##SeqMat Data
stacked_df_url = "https://dholk.primate.wisc.edu/_webdav/dho/sequencing/Polypeptide%20Microarrays/public/COVID_19/%40files/aggregated_data/df_stacked.tsv.gz?contentDisposition=attachment"
stacked_df_path = "df_stacked.tsv.gz"
options(timeout = max(300, getOption("timeout")))
download.file(stacked_df_url, stacked_df_path);

seq_mat = UW.Adult.Covid.19::loadSeqMat(file_name = stacked_df_path);

sample_meta = attr(seq_mat, "sample_meta");
```

### Process data

``` r


## Quantile normalize
seq_mat_qn = UW.Adult.Covid.19::normalizeQuantile(seq_mat[,-1])

probe_meta = probe_meta[probe_meta$PROBE_SEQUENCE %in% rownames(seq_mat_qn),]
seq_mat_qn = seq_mat_qn[rownames(seq_mat_qn) %in% probe_meta$PROBE_SEQUENCE,]


## Create pData data.frame
create_pData<-function(mat_in) {
  
  
  pData = data.frame(
    Sample_ID = colnames(mat_in),
    ptid = colnames(mat_in),
    visit = "pre",
    condition = "Control",
    stringsAsFactors=FALSE
  );
  pData$condition[pData$Sample_ID %in% sample_meta$SAMPLE_NAME[sample_meta$COVID_POSITIVE == "YES"]] = "COVID";
  pData$visit[pData$condition == "COVID"] = "post"
  pData$TAG = pData$Sample_ID;
  
  return(pData);
}

pData = create_pData(seq_mat_qn)
```

### Calculate probe-level pvalues

``` r
library(HERON)
probe_pvalue_res = HERON::calcProbePValuesSeqMat(seq_mat_qn, probe_meta, pData, t.abs_shift = 1, use="t")
#> differential t-test
#> adjusting using BH
```

<img src="man/figures/README-example_probe_pvalues-1.png" width="100%" />

    #> Generating seq_to_probe
    #> Generating seq_to_probe

<img src="man/figures/README-example_probe_pvalues-2.png" width="100%" />

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
