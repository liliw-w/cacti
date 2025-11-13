# write DESCRIPTION
usethis::use_description(fields = list(
  Title       = "CACTI: Chromatin QTL Mapping with Peak Windows",
  Description = "Implements the CACTI peak-window method for multi-peak chromatin QTL mapping.",
  License     = "MIT + file LICENSE",
  Encoding    = "UTF-8",
  LazyData    = "true"
))

usethis::use_author(
  given = "Lili",
  family = "Wang",
  email = "liliw1@uchicago.edu",
  role = c("aut", "cre")
)

usethis::use_package("data.table", type = "Imports")
usethis::use_package("dplyr",      type = "Imports")
usethis::use_package("tidyr",      type = "Imports")
usethis::use_package("tibble",     type = "Imports")
usethis::use_package("stringr",    type = "Imports")
usethis::use_package("ACAT",       type = "Imports")
usethis::use_package("qvalue",     type = "Imports")
usethis::use_package("readr",      type = "Imports")



# add a folder for test data and results
# Example data that ships with the package (for users / vignettes)
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
dir.create("inst/extdata/test_results", recursive = TRUE, showWarnings = FALSE)

# to access the files
pheno_meta_file <- system.file("extdata", "test_pheno_meta.bed", package = "cacti")
pheno_file <- system.file("extdata", "test_pheno.txt", package = "cacti")
cov_file <- system.file("extdata", "test_covariates.txt", package = "cacti")
qtl_file <- system.file("extdata", "test_qtl_sum_stats_chr5.txt.gz", package = "cacti")

system.file("extdata", "test_results/test_pval_window50kb_chr5.txt.gz", package = "cacti")


# write man docs and NAMESPACE
usethis::use_roxygen_md()
devtools::document()


# add vignette
# usethis::use_vignette("cacti_peak_window")
devtools::build_vignettes()


# build package
devtools::build()
devtools::install(build_vignettes = TRUE)



# test package
library(cacti)

?cacti
?cacti_run_chr
?cacti_run_genome
?cacti_group_peak_window

devtools::check()

browseVignettes("cacti")
vignette("cacti_peak_window", package = "cacti")




