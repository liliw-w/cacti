# Residualize peak intensity with covariates

Reads a phenotype matrix (features x samples; first column = feature ID)
and a covariate matrix (covariates x samples; first column = covariate
ID)

## Usage

``` r
cacti_pheno_cov_residual(file_pheno, file_cov, file_pheno_cov_residual)
```

## Arguments

- file_pheno:

  Path to phenotype expression matrix. Rows = features; first column =
  feature ID; remaining columns = samples.

- file_cov:

  Path to covariate matrix. Rows = covariates; first column = covariate
  ID; remaining columns = samples.

- file_pheno_cov_residual:

  Output filename. Rows = features; first column = feature ID; remaining
  columns = samples.

## Value

Invisibly returns the residual matrix (features x samples); writes
`file_pheno_cov_residual`.
