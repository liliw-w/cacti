## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## ----show-pheno-meta, message=FALSE-------------------------------------------
library(cacti)

file_pheno_meta <- system.file(
  "extdata", "test_pheno_meta.bed",
  package = "cacti"
)

head(data.table::fread(file_pheno_meta))

## ----show-pheno, message=FALSE------------------------------------------------
file_pheno <- system.file(
  "extdata", "test_pheno.txt",
  package = "cacti"
)

head(data.table::fread(file_pheno))[, 1:5]

## ----show-cov, message=FALSE--------------------------------------------------
file_cov <- system.file(
  "extdata", "test_covariates.txt",
  package = "cacti"
)

head(data.table::fread(file_cov))[, 1:5]

## ----show-qtl, message=FALSE--------------------------------------------------
qtl_file <- system.file(
  "extdata", "test_qtl_sum_stats_chr5.txt.gz",
  package = "cacti"
)

head(data.table::fread(qtl_file))

## ----example-run-chr5, eval=TRUE----------------------------------------------
library(cacti)

file_pheno_meta <- system.file(
  "extdata", "test_pheno_meta.bed",
  package = "cacti"
)

file_pheno <- system.file(
  "extdata", "test_pheno.txt",
  package = "cacti"
)

file_cov <- system.file(
  "extdata", "test_covariates.txt",
  package = "cacti"
)

qtl_file <- system.file(
  "extdata", "test_qtl_sum_stats_chr5.txt.gz",
  package = "cacti"
)

# Use a temporary prefix to avoid cluttering the working directory
out_prefix <- tempfile("cacti_chr5_")

res_single_chr <- cacti_run_chr(
  window_size     = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno      = file_pheno,
  file_cov        = file_cov,
  chr             = "chr5",
  qtl_file        = qtl_file,
  out_prefix      = out_prefix,
  min_peaks       = 2
)

res_single_chr

## ----example-run-genome, eval=TRUE--------------------------------------------
file_pheno_meta <- system.file(
  "extdata", "test_pheno_meta.bed",
  package = "cacti"
)

file_pheno <- system.file(
  "extdata", "test_pheno.txt",
  package = "cacti"
)

file_cov <- system.file(
  "extdata", "test_covariates.txt",
  package = "cacti"
)

qtl_file <- system.file(
  "extdata", "test_qtl_sum_stats_chr5.txt.gz",
  package = "cacti"
)

out_prefix <- tempfile("cacti_genome_")

res_genome <- cacti_run_genome(
  window_size     = "50kb",
  file_pheno_meta = file_pheno_meta,
  file_pheno      = file_pheno,
  file_cov        = file_cov,
  chrs            = "chr5",                 # toy example: one chromosome
  qtl_files       = qtl_file,
  out_prefix      = out_prefix,
  file_fdr_out    = file.path(tempdir(), "cacti_fdr_chr5.txt.gz")
)

res_genome

