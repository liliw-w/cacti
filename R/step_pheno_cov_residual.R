#' Residualize peak intensity with covariates
#'
#' Reads a phenotype matrix (features x samples; first column = feature ID)
#' and a covariate matrix (covariates x samples; first column = covariate ID)
#' 
#' @param file_pheno Path to phenotype expression matrix. 
#' Rows = features; first column = feature ID; remaining columns = samples.
#' @param file_cov Path to covariate matrix. 
#' Rows = covariates; first column = covariate ID; remaining columns = samples.
#' @param file_pheno_cov_residual Output filename. 
#' Rows = features; first column = feature ID; remaining columns = samples.
#'
#' @return Invisibly returns the residual matrix (features x samples); writes `file_pheno_cov_residual`.
#' @export
cacti_pheno_cov_residual <- function(file_pheno, file_cov, file_pheno_cov_residual) {
  suppressPackageStartupMessages(requireNamespace("readr"))
  suppressPackageStartupMessages(requireNamespace("data.table"))
  suppressPackageStartupMessages(requireNamespace("dplyr"))

  # read
  pheno_df <- data.table::fread(file_pheno, header = TRUE)
  cov_df <- data.table::fread(file_cov, header = TRUE)
  
  colnames(pheno_df)[1] <- "ID"
  colnames(cov_df)[1]   <- "ID"
  
  # sample order from input phenotype (to preserve in output)
  sample_order <- setdiff(colnames(pheno_df), "ID")

  # matrices for regression
  pheno_mat <- dplyr::select(pheno_df, -ID) |> as.matrix() |> t()   # samples x features
  colnames(pheno_mat) <- pheno_df$ID                                 # feature IDs become columns

  cov_mat <- dplyr::select(cov_df, -ID) |> as.matrix() |> t()        # samples x covariates
  colnames(cov_mat) <- cov_df$ID

  # align samples (rows) between phenotype and covariates
  if (!all(rownames(pheno_mat) %in% rownames(cov_mat))) {
    stop("Some phenotype samples are missing in covariates.")
  }
  cov_mat <- cov_mat[match(rownames(pheno_mat), rownames(cov_mat)), , drop = FALSE]
  if (any(is.na(rownames(cov_mat)))) stop("Sample alignment failed.")

  # residualize each feature: lm(y ~ covariates)
  extract_residual <- function(y, X) stats::lm(y ~ X)$residuals
  res_s_by_f <- apply(pheno_mat, 2, function(y) extract_residual(y, cov_mat))  # samples x features

  # restore dimnames
  rownames(res_s_by_f) <- rownames(pheno_mat)   # samples
  colnames(res_s_by_f) <- colnames(pheno_mat)   # features

  # transpose back to features x samples, and order samples as original file_pheno
  res_f_by_s <- t(res_s_by_f)
  res_f_by_s <- res_f_by_s[, sample_order, drop = FALSE]

  # build output table in the same format as file_pheno
  out_df <- data.frame(ID = rownames(res_f_by_s), check.names = FALSE)
  out_df <- cbind(out_df, as.data.frame(res_f_by_s, check.names = FALSE))

  # write TSV
  data.table::fwrite(out_df, file = file_pheno_cov_residual, sep = '\t', quote = FALSE, col.names = TRUE)

  invisible(out_df)
}
