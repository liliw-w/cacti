#' Add FDR (q-values) to per-window top hits across chromosomes
#'
#' Reads per-chromosome CACTI window result files, computes per-window
#' top-hit summaries (min p, Bonferroni, ACAT across SNPs in the window),
#' and assigns q-values to the ACAT top-hit p-values.
#'
#' Expected columns in each file:
#' - `group` : window / peak-group ID
#' - `snp`   : variant ID
#' - `pval` : p-value per (snp, group)
#'
#' @param file_all_pval a character vector containing all pvalue file paths.
#'   e.g., "c(chr1.txt.gz, chr2.txt.gz, chr3.txt.gz)"
#' @param file_fdr_out Output filename for combined results with FDR.
#'
#' @return Invisibly returns the combined results with FDR.; writes `file_fdr_out`.
#' @export
cacti_add_fdr <- function(file_all_pval, file_fdr_out) {
  suppressPackageStartupMessages(requireNamespace("data.table"))
  suppressPackageStartupMessages(requireNamespace("dplyr"))
  suppressPackageStartupMessages(requireNamespace("stringr"))
  suppressPackageStartupMessages(requireNamespace("ACAT"))
  suppressPackageStartupMessages(requireNamespace("qvalue"))
  
  # --- parse input string into file paths ---
  if (!length(file_all_pval)) {
    stop("No input files detected.")
  }

  # --- process each file sequentially ---
  dt_fdr <- lapply(
    X = file_all_pval,
    FUN = function(x) {
      # per-group adjustment
      dt <- data.table::fread(x) |>
        dplyr::group_by(group) |>
        dplyr::mutate(n_var_in_cis = dplyr::n()) |>
        dplyr::ungroup()

      # correct p == 1 for ACAT stability
      dt$pval[dt$pval == 1] <- 1 - 1 / dt$n_var_in_cis[dt$pval == 1]

      # summarise per window
      dplyr::group_by(dt, group) |>
        dplyr::summarise(
          best_hit = min(pval, na.rm = TRUE),
          n_var_in_cis = dplyr::n(),
          best_hit_bonf = min(1, best_hit * n_var_in_cis),
          best_hit_acat = ACAT::ACAT(Pvals = pval),
          snp = snp[which.min(pval)],
          .groups = "drop"
        )
    }) |>
    dplyr::bind_rows()

  # --- compute q-values ---
  dt_fdr$q <- signif(qvalue::qvalue(dt_fdr$best_hit_acat)$qvalues, 5)

  data.table::fwrite(dt_fdr, file = file_fdr_out, sep = "\t", quote = FALSE, col.names = TRUE)
  invisible(dt_fdr)
}
