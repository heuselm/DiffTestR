#' Assess Proteomics Data Quality
#'
#' @description Comprehensive quality assessment of proteomics data including
#' coefficient of variation analysis, missing value summaries, dynamic range,
#' and replicate correlation metrics. Includes leave-one-out CV analysis to
#' identify outlier replicates.
#'
#' @param data_matrix Numeric matrix of log2-transformed intensities with
#' precursors/peptides in rows and samples in columns
#' @param study_design Data.table or data.frame with columns 'filename',
#' 'condition', and 'replicate'. Rownames or 'filename' column must match
#' column names in data_matrix
#' @param cv_threshold Threshold for flagging outlier replicates based on
#' CV improvement upon removal (default: 0.3 = 30% improvement)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A list with the following components:
#' \itemize{
#'   \item cv_stats: CV statistics per condition
#'   \item cv_loo_analysis: Leave-one-out CV analysis results
#'   \item outlier_replicates: Data.table of flagged outlier replicates
#'   \item missing_stats: Missing value analysis
#'   \item dynamic_range: Min, max, and range of intensities
#'   \item replicate_correlations: Correlation matrix within conditions
#' }
#'
#' @details
#' The leave-one-out CV analysis calculates the coefficient of variation for
#' each condition with all replicates, then recalculates CV after removing
#' each replicate individually. If removing a replicate improves the median
#' CV by more than cv_threshold (default 30%), that replicate is flagged as
#' a potential outlier.
#'
#' @examples
#' \dontrun{
#' data(DiffTestR_example_data_wide)
#' data(DiffTestR_example_study_design)
#'
#' # Convert to matrix
#' data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
#' rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id
#'
#' # Assess quality
#' qc <- assessDataQuality(data_mat, DiffTestR_example_study_design)
#'
#' # View outlier replicates
#' print(qc$outlier_replicates)
#'
#' # Plot CV distributions
#' plot(qc$cv_stats$cv_distribution)
#' }
#'
#' @author Moritz Heusel
#'
#' @import data.table
#' @importFrom stats cor median sd
#'
#' @export
assessDataQuality <- function(data_matrix,
                              study_design,
                              cv_threshold = 0.3,
                              verbose = TRUE) {

  # Validate inputs
  if (!is.matrix(data_matrix)) {
    stop("data_matrix must be a numeric matrix")
  }

  study_design <- as.data.table(study_design)

  # Ensure study design has required columns
  required_cols <- c("filename", "condition", "replicate")
  if (!all(required_cols %in% names(study_design))) {
    stop("study_design must have columns: ", paste(required_cols, collapse = ", "))
  }

  if (verbose) {
    print_analysis_header("Data Quality Assessment")
    report_progress("Starting quality assessment", type = "info")
  }

  # Helper function to extract basename (remove path and extension)
  extract_basename <- function(x) {
    # Remove file extension (common ones: .raw, .d, .tsv, .csv, etc.)
    x <- gsub("\\.(raw|d|tsv|csv|txt)$", "", x, ignore.case = TRUE)
    # Extract basename (last component after path separators)
    # Handle both Windows (\) and Unix (/) paths
    x <- gsub(".*[/\\\\]", "", x)
    return(x)
  }

  # Match column names to study design
  sample_cols <- colnames(data_matrix)
  
  # Extract basenames from column names
  sample_basenames <- extract_basename(sample_cols)
  
  # Extract basenames from study_design$filename (in case they also have paths)
  study_design_basenames <- extract_basename(study_design$filename)
  
  # Create mapping: column name -> study_design filename
  # Match by basename
  col_to_filename <- character(length(sample_cols))
  names(col_to_filename) <- sample_cols
  
  for (i in seq_along(sample_cols)) {
    # Find matching study_design entry by basename
    match_idx <- which(study_design_basenames == sample_basenames[i])
    
    if (length(match_idx) == 0) {
      stop(sprintf(
        "Column '%s' (basename: '%s') does not match any filename in study_design",
        sample_cols[i], sample_basenames[i]
      ))
    } else if (length(match_idx) > 1) {
      warning(sprintf(
        "Column '%s' (basename: '%s') matches multiple entries in study_design. Using first match.",
        sample_cols[i], sample_basenames[i]
      ))
      match_idx <- match_idx[1]
    }
    
    col_to_filename[i] <- study_design$filename[match_idx]
  }
  
  # Check that all study_design entries are matched
  unmatched <- setdiff(study_design$filename, col_to_filename)
  if (length(unmatched) > 0) {
    warning(sprintf(
      "Some study_design filenames not found in data_matrix: %s",
      paste(unmatched, collapse = ", ")
    ))
  }
  
  # Rename columns in data_matrix to use study_design$filename (shorter names)
  # Only rename if they're different
  if (!identical(sample_cols, col_to_filename)) {
    if (verbose) {
      report_progress("Renaming columns from paths to study_design filenames", type = "info")
    }
    colnames(data_matrix) <- col_to_filename
    sample_cols <- col_to_filename  # Update for use in rest of function
  }

  #---------------------------------------------------------------------------
  # 1. CV Analysis per Condition
  #---------------------------------------------------------------------------
  if (verbose) {
    print_section_header("Coefficient of Variation Analysis")
  }

  cv_stats_list <- list()
  conditions <- unique(study_design$condition)

  for (cond in conditions) {
    # Get samples for this condition
    cond_samples <- study_design[condition == cond, filename]
    cond_data <- data_matrix[, cond_samples, drop = FALSE]

    # Calculate CV per precursor
    cv_values <- apply(cond_data, 1, function(x) {
      x_valid <- x[!is.na(x)]
      if (length(x_valid) < 2) return(NA)
      (sd(x_valid) / mean(x_valid)) * 100
    })

    cv_stats_list[[cond]] <- list(
      condition = cond,
      n_samples = length(cond_samples),
      median_cv = median(cv_values, na.rm = TRUE),
      mean_cv = mean(cv_values, na.rm = TRUE),
      sd_cv = sd(cv_values, na.rm = TRUE),
      cv_values = cv_values
    )

    if (verbose) {
      report_progress(sprintf("Condition '%s': Median CV = %.1f%%",
                             cond, cv_stats_list[[cond]]$median_cv),
                     type = "info")
    }
  }

  #---------------------------------------------------------------------------
  # 2. Leave-One-Out CV Analysis
  #---------------------------------------------------------------------------
  if (verbose) {
    print_section_header("Leave-One-Out Replicate Analysis")
  }

  loo_results <- list()
  outlier_flags <- list()

  for (cond in conditions) {
    cond_samples <- study_design[condition == cond, filename]
    n_reps <- length(cond_samples)

    if (n_reps < 3) {
      if (verbose) {
        report_progress(sprintf("Condition '%s': < 3 replicates, skipping LOO analysis",
                               cond), type = "warning")
      }
      next
    }

    cond_data <- data_matrix[, cond_samples, drop = FALSE]

    # Calculate baseline CV (all replicates)
    baseline_cv <- apply(cond_data, 1, function(x) {
      x_valid <- x[!is.na(x)]
      if (length(x_valid) < 2) return(NA)
      (sd(x_valid) / mean(x_valid)) * 100
    })
    baseline_median_cv <- median(baseline_cv, na.rm = TRUE)

    # Leave-one-out analysis
    loo_cv_results <- data.table()

    for (i in seq_along(cond_samples)) {
      # Remove one replicate
      loo_samples <- cond_samples[-i]
      loo_data <- cond_data[, loo_samples, drop = FALSE]

      # Calculate CV without this replicate
      loo_cv <- apply(loo_data, 1, function(x) {
        x_valid <- x[!is.na(x)]
        if (length(x_valid) < 2) return(NA)
        (sd(x_valid) / mean(x_valid)) * 100
      })

      loo_median_cv <- median(loo_cv, na.rm = TRUE)

      # Calculate improvement
      cv_improvement <- (baseline_median_cv - loo_median_cv) / baseline_median_cv

      loo_cv_results <- rbind(loo_cv_results, data.table(
        condition = cond,
        replicate = study_design[filename == cond_samples[i], replicate],
        filename = cond_samples[i],
        baseline_cv = baseline_median_cv,
        loo_cv = loo_median_cv,
        cv_improvement = cv_improvement,
        is_outlier = cv_improvement > cv_threshold
      ))
    }

    loo_results[[cond]] <- loo_cv_results

    # Report outliers
    outliers <- loo_cv_results[is_outlier == TRUE]
    if (nrow(outliers) > 0) {
      for (j in seq_len(nrow(outliers))) {
        if (verbose) {
          report_progress(sprintf("Outlier detected: %s (CV improves by %.1f%% when removed)",
                                 outliers[j, filename],
                                 outliers[j, cv_improvement] * 100),
                         type = "warning")
        }
      }
      outlier_flags[[cond]] <- outliers
    } else {
      if (verbose) {
        report_progress(sprintf("Condition '%s': No outlier replicates detected", cond),
                       type = "success")
      }
    }
  }

  #---------------------------------------------------------------------------
  # 3. Missing Value Analysis
  #---------------------------------------------------------------------------
  if (verbose) {
    print_section_header("Missing Value Analysis")
  }

  total_values <- length(data_matrix)
  n_missing <- sum(is.na(data_matrix))
  pct_missing <- 100 * n_missing / total_values

  missing_per_sample <- colSums(is.na(data_matrix)) / nrow(data_matrix) * 100
  missing_per_precursor <- rowSums(is.na(data_matrix)) / ncol(data_matrix) * 100

  missing_stats <- list(
    total_values = total_values,
    n_missing = n_missing,
    pct_missing = pct_missing,
    missing_per_sample = missing_per_sample,
    missing_per_precursor = missing_per_precursor,
    samples_with_high_missing = names(missing_per_sample[missing_per_sample > 50])
  )

  if (verbose) {
    report_progress(sprintf("Overall missing values: %.1f%%", pct_missing),
                   type = "info")
    if (length(missing_stats$samples_with_high_missing) > 0) {
      report_progress(sprintf("Samples with >50%% missing: %s",
                             paste(missing_stats$samples_with_high_missing, collapse = ", ")),
                     type = "warning")
    }
  }

  #---------------------------------------------------------------------------
  # 4. Dynamic Range
  #---------------------------------------------------------------------------
  if (verbose) {
    print_section_header("Dynamic Range")
  }

  dynamic_range <- list(
    min = min(data_matrix, na.rm = TRUE),
    max = max(data_matrix, na.rm = TRUE),
    range = max(data_matrix, na.rm = TRUE) - min(data_matrix, na.rm = TRUE)
  )

  if (verbose) {
    report_progress(sprintf("Dynamic range: %.2f to %.2f (range: %.2f log2 units)",
                           dynamic_range$min, dynamic_range$max,
                           dynamic_range$range),
                   type = "info")
  }

  #---------------------------------------------------------------------------
  # 5. Replicate Correlation
  #---------------------------------------------------------------------------
  if (verbose) {
    print_section_header("Replicate Correlation")
  }

  replicate_cors <- list()

  for (cond in conditions) {
    cond_samples <- study_design[condition == cond, filename]
    cond_data <- data_matrix[, cond_samples, drop = FALSE]

    # Calculate correlation matrix (pairwise complete observations)
    cor_mat <- cor(cond_data, use = "pairwise.complete.obs")

    # Get off-diagonal correlations
    cor_values <- cor_mat[lower.tri(cor_mat)]

    replicate_cors[[cond]] <- list(
      cor_matrix = cor_mat,
      median_cor = median(cor_values, na.rm = TRUE),
      mean_cor = mean(cor_values, na.rm = TRUE),
      min_cor = min(cor_values, na.rm = TRUE),
      max_cor = max(cor_values, na.rm = TRUE)
    )

    if (verbose) {
      report_progress(sprintf("Condition '%s': Median correlation = %.3f",
                             cond, replicate_cors[[cond]]$median_cor),
                     type = "info")
    }
  }

  #---------------------------------------------------------------------------
  # Compile results
  #---------------------------------------------------------------------------
  if (verbose) {
    report_progress("Quality assessment complete", type = "success")
  }

  # Combine all outlier flags
  outlier_dt <- rbindlist(outlier_flags, fill = TRUE)

  result <- list(
    cv_stats = cv_stats_list,
    cv_loo_analysis = loo_results,
    outlier_replicates = outlier_dt,
    missing_stats = missing_stats,
    dynamic_range = dynamic_range,
    replicate_correlations = replicate_cors
  )

  class(result) <- c("DiffTestR_QC", "list")
  return(result)
}

#' Print method for DiffTestR_QC objects
#'
#' @param x A DiffTestR_QC object
#' @param ... Additional arguments (not used)
#'
#' @export
print.DiffTestR_QC <- function(x, ...) {
  cat("\n=== DiffTestR Data Quality Report ===\n\n")

  cat("Coefficient of Variation:\n")
  for (cond_name in names(x$cv_stats)) {
    stats <- x$cv_stats[[cond_name]]
    cat(sprintf("  %s: Median CV = %.1f%%, Mean CV = %.1f%%\n",
                cond_name, stats$median_cv, stats$mean_cv))
  }

  cat("\nOutlier Replicates:\n")
  if (nrow(x$outlier_replicates) > 0) {
    print(x$outlier_replicates[, .(condition, filename, cv_improvement)])
  } else {
    cat("  No outliers detected\n")
  }

  cat("\nMissing Values:\n")
  cat(sprintf("  Overall: %.1f%% missing\n", x$missing_stats$pct_missing))

  cat("\nDynamic Range:\n")
  cat(sprintf("  %.2f to %.2f log2 units (range: %.2f)\n",
              x$dynamic_range$min, x$dynamic_range$max, x$dynamic_range$range))

  cat("\nReplicate Correlation:\n")
  for (cond_name in names(x$replicate_correlations)) {
    cors <- x$replicate_correlations[[cond_name]]
    cat(sprintf("  %s: Median r = %.3f\n", cond_name, cors$median_cor))
  }

  cat("\n")
  invisible(x)
}
