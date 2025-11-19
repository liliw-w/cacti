# Add FDR (q-values) to per-window top hits across chromosomes

Reads per-chromosome CACTI window result files, computes per-window
top-hit summaries (min p, Bonferroni, ACAT across SNPs in the window),
and assigns q-values to the ACAT top-hit p-values.

## Usage

``` r
cacti_add_fdr(file_all_pval, file_fdr_out)
```

## Arguments

- file_all_pval:

  a character vector containing all pvalue file paths. e.g.,
  "c(chr1.txt.gz, chr2.txt.gz, chr3.txt.gz)"

- file_fdr_out:

  Output filename for combined results with FDR.

## Value

Invisibly returns the combined results with FDR.; writes `file_fdr_out`.

## Details

Expected columns in each file:

- `group` : window / peak-group ID

- `snp` : variant ID

- `pval` : p-value per (snp, group)
