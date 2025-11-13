#' CACTI: Chromatin QTL Mapping with Peak Windows
#'
#' The `cacti` package implements the CACTI peak-window method for chromatin QTL
#' mapping. It is designed to:
#'
#' \itemize{
#'   \item group nearby peaks into non-overlapping genomic windows,
#'   \item residualize peak intensities with respect to covariates,
#'   \item run multi-peak association tests using a principal-component omnibus (PCO) test, and
#'   \item aggregate per-chromosome results and compute window-level FDR.
#' }
#'
#' The main user-facing functions are:
#'
#' \itemize{
#'   \item \code{\link{cacti_run_chr}}: run the peak-window pipeline for a single chromosome.
#'   \item \code{\link{cacti_run_genome}}: run the peak-window pipeline across multiple chromosomes and add FDR.
#' }
#'
#' Lower-level building blocks:
#'
#' \itemize{
#'   \item \code{\link{cacti_group_peak_window}}: group peaks into non-overlapping windows.
#'   \item \code{\link{cacti_pheno_cov_residual}}: residualize phenotypes by covariates.
#'   \item \code{\link{cacti_cal_p}}: run CACTI association tests on peak windows.
#'   \item \code{\link{cacti_add_fdr}}: compute window-level FDR across chromosomes.
#' }
#'
#' A small toy dataset from Aracena et al. is bundled under
#' \code{inst/extdata/} and can be accessed using
#' \code{\link[base]{system.file}} for testing and examples.
#'
#' For an overview of the peak-window pipeline and file formats, see the
#' package vignette:
#'
#' \code{browseVignettes("cacti")}
#'
#' @docType package
#' @name cacti
"_PACKAGE"
