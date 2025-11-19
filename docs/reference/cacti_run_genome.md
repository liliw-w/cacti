# Run CACTI peak-window pipeline genome-wide and add FDR

This function is a convenience wrapper that runs the CACTI peak-window
pipeline across multiple chromosomes and then computes window-level FDR
across all peak windows and chromosomes.

## Usage

``` r
cacti_run_genome(
  window_size,
  file_pheno_meta,
  file_pheno,
  file_cov,
  chrs,
  qtl_files,
  out_prefix,
  dir_pco = system.file("pco", package = "cacti"),
  min_peaks = 2,
  file_fdr_out = NULL
)
```

## Arguments

- window_size:

  Character like `"50kb"` (kb units) or numeric (bp).

- file_pheno_meta:

  Path to input phenotype meta table (with header). Columns are -

  - `phe_chr`: chromosome, with "chr" prefix (e.g. "chr1")

  - `phe_from`: peak start (integer, \< `phe_to`)

  - `phe_to`: peak end (integer)

  - `phe_id`: peak identifier (unique per row)

- file_pheno:

  Path to phenotype expression matrix. Rows = features; first column =
  feature ID; remaining columns = samples.

- file_cov:

  Path to covariate matrix. Rows = covariates; first column = covariate
  ID; remaining columns = samples.

- chrs:

  Character vector of chromosome labels (e.g., `paste0("chr", 1:22)`).

- qtl_files:

  Either:

  - a character vector of the same length as `chrs`, giving the cis-QTL
    file path for each chromosome in order; or

  - a single template string containing the placeholder `"{chr}"`, e.g.,
    `"extdata/test_qtl_sum_stats_{chr}.txt.gz"`. In that case, the
    placeholder is replaced by each element of `chrs`.

- out_prefix:

  Output prefix used to construct all output filenames.

- dir_pco:

  Directory containing association test helpers:
  ModifiedPCOMerged_acat.R, liu.R, liumod.R, davies.R, qfc.so.

- min_peaks:

  Minimum number of peaks required in a window to run the multivariate
  PCO test (\>= min_peaks -\> PCO; \< min_peaks -\> univariate p).

- file_fdr_out:

  Optional output path for the FDR-added window-level file. If `NULL`, a
  default filename is constructed from `out_prefix` and `window_size`.

## Value

Invisibly returns a named list of output paths with elements:

- file_peak_group:

  Path to the window-level group.

- file_peak_group_peaklevel:

  Path to the peak-level group.

- file_pheno_residual:

  Path to the residualized phenotype matrix.

- file_p_peak_group:

  Path to the per-window p-value file for all chromosome.

- file_fdr_out:

  Path to the FDR-added window-level result file.

## Details

For each chromosome in `chrs`, it calls
[`cacti_run_chr()`](https://liliw-w.github.io/cacti/reference/cacti_run_chr.md)
and collects the per-window p-value files. It then calls
[`cacti_add_fdr()`](https://liliw-w.github.io/cacti/reference/cacti_add_fdr.md)
once, aggregating all chromosomes to obtain q-values for the top-hit
p-values in each window.

## Examples
