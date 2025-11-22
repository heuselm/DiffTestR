#' Print and Summary Methods for DiffTestR Objects
#'
#' @description S3 methods for printing and summarizing DiffTestR analysis results
#'
#' @name print_methods
NULL

#' Print method for diffExpr objects
#'
#' @param x A diffExpr object (result from testDifferentialAbundance)
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.diffExpr <- function(x, ...) {
  cat("\n")
  cat("====================================\n")
  cat("  Differential Abundance Results\n")
  cat("====================================\n\n")

  # Extract summary stats if available
  if ("summary_stats" %in% names(x)) {
    stats <- x$summary_stats

    cat("Input Data:\n")
    cat(sprintf("  Precursors (input): %d\n", stats$n_precursors_input))
    cat(sprintf("  Precursors (tested): %d\n", stats$n_precursors_tested))
    cat(sprintf("  Proteins (tested): %d\n", stats$n_proteins_tested))

    cat("\nSignificant Results (p < 0.05):\n")
    cat(sprintf("  Precursors: %d (%.1f%%)\n",
                stats$n_significant_precursors,
                100 * stats$n_significant_precursors / stats$n_precursors_tested))
    cat(sprintf("  Proteins: %d (%.1f%%)\n",
                stats$n_significant_proteins,
                100 * stats$n_significant_proteins / stats$n_proteins_tested))

    if (!is.null(stats$n_upregulated) && !is.null(stats$n_downregulated)) {
      cat("\nRegulation (|log2FC| > 1, p < 0.05):\n")
      cat(sprintf("  Up-regulated: %d proteins\n", stats$n_upregulated))
      cat(sprintf("  Down-regulated: %d proteins\n", stats$n_downregulated))
    }

    if (!is.null(stats$median_log2fc_protein)) {
      cat("\nFold Changes:\n")
      cat(sprintf("  Median log2FC (protein): %.3f\n", stats$median_log2fc_protein))
    }

    if (!is.null(stats$median_cv_condition1) && !is.null(stats$median_cv_condition2)) {
      cat("\nData Quality (CV %):\n")
      cat(sprintf("  Condition 1: %.1f%%\n", stats$median_cv_condition1))
      cat(sprintf("  Condition 2: %.1f%%\n", stats$median_cv_condition2))
    }

    if (!is.null(stats$comparison)) {
      cat("\nComparison:\n")
      cat(sprintf("  %s vs %s\n", stats$condition_1, stats$condition_2))
    }

  } else {
    # Fallback if no summary stats
    if ("diffExpr_result_dt" %in% names(x)) {
      dt <- x$diffExpr_result_dt
      cat(sprintf("Precursors tested: %d\n", nrow(dt)))
      cat(sprintf("Proteins tested: %d\n",
                  length(unique(dt$Protein.Group))))

      n_sig <- sum(dt$p_value_BHadj_protein <= 0.05, na.rm = TRUE)
      cat(sprintf("Significant proteins (p<0.05): %d\n", n_sig))
    }
  }

  cat("\n")
  cat("Use summary(x) for more details\n")
  cat("Access results: x$diffExpr_result_dt\n")
  cat("====================================\n\n")

  invisible(x)
}

#' Summary method for diffExpr objects
#'
#' @param object A diffExpr object
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns summary statistics
#'
#' @export
summary.diffExpr <- function(object, ...) {
  print(object)

  if ("summary_stats" %in% names(object)) {
    cat("\n=== Detailed Summary Statistics ===\n\n")

    stats <- object$summary_stats

    # Print all available stats
    for (stat_name in names(stats)) {
      if (!is.null(stats[[stat_name]]) && length(stats[[stat_name]]) == 1) {
        cat(sprintf("%-30s: %s\n", stat_name,
                    format(stats[[stat_name]], digits = 4)))
      }
    }
  }

  invisible(object$summary_stats)
}

#' Calculate summary statistics for differential abundance results
#'
#' @param diffExpr_result_dt Data.table with differential expression results
#' @param input_n_precursors Number of precursors in input data
#' @param condition_1 Name of condition 1
#' @param condition_2 Name of condition 2
#' @param data_matrix Optional: matrix for CV calculation
#' @param study_design Optional: study design for CV per condition
#'
#' @return List of summary statistics
#'
#' @keywords internal
#' @noRd
calculate_summary_stats <- function(diffExpr_result_dt,
                                   input_n_precursors,
                                   condition_1 = NULL,
                                   condition_2 = NULL,
                                   data_matrix = NULL,
                                   study_design = NULL) {

  stats <- list(
    n_precursors_input = input_n_precursors,
    n_precursors_tested = nrow(diffExpr_result_dt),
    n_proteins_tested = length(unique(diffExpr_result_dt$Protein.Group)),
    n_significant_precursors = sum(diffExpr_result_dt$p_value_BHadj <= 0.05,
                                   na.rm = TRUE),
    n_significant_proteins = length(unique(
      diffExpr_result_dt[p_value_BHadj_protein <= 0.05, Protein.Group]
    )),
    median_log2fc_precursor = median(diffExpr_result_dt$log2_fold_change,
                                     na.rm = TRUE),
    median_log2fc_protein = median(diffExpr_result_dt$log2_fold_change_protein,
                                   na.rm = TRUE),
    condition_1 = condition_1,
    condition_2 = condition_2,
    comparison = if (!is.null(condition_1) && !is.null(condition_2)) {
      paste(condition_1, "vs", condition_2)
    } else {
      NULL
    }
  )

  # Calculate regulation counts (|log2FC| > 1, p < 0.05)
  sig_proteins <- diffExpr_result_dt[p_value_BHadj_protein <= 0.05]
  stats$n_upregulated <- length(unique(
    sig_proteins[log2_fold_change_protein > 1, Protein.Group]
  ))
  stats$n_downregulated <- length(unique(
    sig_proteins[log2_fold_change_protein < -1, Protein.Group]
  ))

  # Calculate CV if data provided
  if (!is.null(data_matrix) && !is.null(study_design)) {
    study_design <- as.data.table(study_design)

    # CV for condition 1
    if (!is.null(condition_1)) {
      cond1_samples <- study_design[condition == condition_1, filename]
      if (length(cond1_samples) > 0) {
        cond1_data <- data_matrix[, cond1_samples, drop = FALSE]
        cv1 <- apply(cond1_data, 1, function(x) {
          x_valid <- x[!is.na(x)]
          if (length(x_valid) < 2) return(NA)
          (sd(x_valid) / mean(x_valid)) * 100
        })
        stats$median_cv_condition1 <- median(cv1, na.rm = TRUE)
      }
    }

    # CV for condition 2
    if (!is.null(condition_2)) {
      cond2_samples <- study_design[condition == condition_2, filename]
      if (length(cond2_samples) > 0) {
        cond2_data <- data_matrix[, cond2_samples, drop = FALSE]
        cv2 <- apply(cond2_data, 1, function(x) {
          x_valid <- x[!is.na(x)]
          if (length(x_valid) < 2) return(NA)
          (sd(x_valid) / mean(x_valid)) * 100
        })
        stats$median_cv_condition2 <- median(cv2, na.rm = TRUE)
      }
    }
  }

  return(stats)
}
