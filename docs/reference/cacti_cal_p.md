# Run CACTI association test for a single chromosome

Run CACTI association test for a single chromosome

## Usage

``` r
cacti_cal_p(
  file_qtl_cis_norm,
  chr,
  file_peak_group,
  file_pheno_cov_residual,
  file_p_peak_group,
  dir_pco = system.file("pco", package = "cacti"),
  min_peaks = 2
)
```

## Arguments

- file_qtl_cis_norm:

  Path to cis QTL summary stats. Must contain columns: phe_id, var_id,
  z.

- chr:

  Chromosome label (e.g., "chr21"). Used to subset peak groups.

- file_peak_group:

  Path to peak-window metadata (groups) file. Produced by
  step_group_peak_window.

- file_pheno_cov_residual:

  Path to residualized phenotype matrix file. Produced by
  step_pheno_cov_residual.

- file_p_peak_group:

  Output path for per-window pvalues. Columns: group, snp, pval.

- dir_pco:

  Directory containing association test helpers:
  ModifiedPCOMerged_acat.R, liu.R, liumod.R, davies.R, qfc.so.

- min_peaks:

  Minimum number of peaks required in a group to run the multivariate
  PCO test (\>= min_peaks -\> PCO; \< min_peaks -\> univariate p).

## Value

Invisibly returns the pvalue tibble; writes `file_p_peak_group`.
