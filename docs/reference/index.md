# Package index

## Overview

- [`cacti-package`](https://liliw-w.github.io/cacti/reference/cacti.md)
  [`cacti`](https://liliw-w.github.io/cacti/reference/cacti.md) : CACTI:
  Chromatin QTL Mapping with Peak Windows

## Main workflows

- [`cacti_run_chr()`](https://liliw-w.github.io/cacti/reference/cacti_run_chr.md)
  : Run CACTI peak-window pipeline for a single chromosome
- [`cacti_run_genome()`](https://liliw-w.github.io/cacti/reference/cacti_run_genome.md)
  : Run CACTI peak-window pipeline genome-wide and add FDR

## Building blocks

- [`cacti_group_peak_window()`](https://liliw-w.github.io/cacti/reference/cacti_group_peak_window.md)
  : Group adjacent peaks by non-overlapping windows
- [`cacti_pheno_cov_residual()`](https://liliw-w.github.io/cacti/reference/cacti_pheno_cov_residual.md)
  : Residualize peak intensity with covariates
- [`cacti_cal_p()`](https://liliw-w.github.io/cacti/reference/cacti_cal_p.md)
  : Run CACTI association test for a single chromosome
- [`cacti_add_fdr()`](https://liliw-w.github.io/cacti/reference/cacti_add_fdr.md)
  : Add FDR (q-values) to per-window top hits across chromosomes
