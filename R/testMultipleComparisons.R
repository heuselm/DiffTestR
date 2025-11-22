#' Test Multiple Differential Abundance Comparisons
#'
#' @description Perform differential abundance testing for multiple pairwise comparisons.
#' This function either accepts explicit comparison pairs or auto-generates all pairwise
#' comparisons from unique conditions in the study design. Results are aggregated and
#' returned as a structured multiDiffExpr object.
#'
#' @param input_dt Input data table either in tsv/txt format or already in R as data.table or data.frame with the following columns:
#' \itemize{
#' \item Protein.Group: Semicolon-separated Uniprot IDs (or similar, as long as it matches)
#' \item Precursor.Id: Unique Precursor Id for which the quantitative values are contained
#' \item "filename": The file names of the MS raw data, must be identical to the entries in study_design$filename
#' }
#' Note: The data will be log2-transformed internally.
#'
#' @param protein_group_annotation Protein annotation table with columns Protein.Group and Protein.Names (and others if desired)
#' that will be used to annotate the results. By default it is assumed to be a subset of and and an attempt will be made
#' to extract it from the input_dt.
#' @param study_design Study design in tab-separated .txt with mandatory columns:
#' \itemize{
#' \item filename: Must match quantitative data-containing column headers in the input_dt
#' \item condition: String, biological condition (e.g. "treated" and "untreated")
#' \item replicate: Replicate number (integer). Minimally 3 replicates are needed per condition for this type of analysis.
#' }
#' @param comparisons Optional list of comparisons to perform. Each element should be a named list with
#' \code{condition_1} and \code{condition_2}. If NULL (default), all pairwise comparisons will be auto-generated.
#' Example: \code{list(list(condition_1 = "A", condition_2 = "B"), list(condition_1 = "A", condition_2 = "C"))}
#' @param normalize_data Whether or not data is scaled/normalized before differential testing. In some cases
#' it might be preferable not to scale the datasets, e.g. when comparing pulldowns vs. input samples! Defaults to TRUE.
#' @param normalization_function Normalization function to use that transforms a matrix of quantities where columns are
#' samples and rows are analytes. Defaults to limma:normalizeQuantiles, but can be replaced with any such function. You may want
#' to try limma::normalizeVSN or limma::normalizeMedianValues.
#' @param min_n_obs Minimum number of observations per precursor (number of runs it was identified in)
#' in order to keep in in the analysis
#' @param imp_percentile Percentile of the total distribution of values on which the random
#' distribution for sampling will be centered
#' @param imp_sd standard deviation of the normal distribution from which values are sampled to impute missing values
#' @param output_dir Output directory for results. If NULL, uses current working directory. Each comparison
#' will be saved in a subdirectory named "condition1_vs_condition2"
#' @param plot_pdf Document processing steps in a string of pdf graphs for each comparison
#' @param write_tsv_tables Write out final quant table with differential expression testing results
#' @param target_protein Optional string with protein identifier to highlight in volcano plots
#' @param stop_on_error Logical, whether to stop on first error (TRUE) or continue with remaining
#' comparisons (FALSE, default). When FALSE, errors are logged in the error_log element.
#' @param verbose Logical, whether to print progress messages during analysis (default: TRUE)
#'
#' @return A multiDiffExpr object (list) containing:
#' \itemize{
#' \item individual_results: List of diffExpr objects, one per comparison
#' \item combined_results_dt: Aggregated data.table with results from all comparisons
#' \item comparisons: Specification of comparisons performed
#' \item n_comparisons: Total number of comparisons attempted
#' \item n_successful: Number of comparisons completed successfully
#' \item error_log: Data.table with any errors encountered (NULL if all successful)
#' \item output_dir: Directory where results were saved
#' }
#'
#' @author Moritz Heusel
#'
#' @examples
#' \dontrun{
#' # Auto-generate all pairwise comparisons
#' result <- testMultipleComparisons(
#'   input_dt = "path/to/DIANN_matrix.tsv",
#'   study_design = "path/to/Study_design.tsv",
#'   output_dir = "multi_comparison_results"
#' )
#'
#' # Explicit comparisons
#' my_comparisons <- list(
#'   list(condition_1 = "Treatment_A", condition_2 = "Control"),
#'   list(condition_1 = "Treatment_B", condition_2 = "Control")
#' )
#' result <- testMultipleComparisons(
#'   input_dt = data_dt,
#'   study_design = design_dt,
#'   comparisons = my_comparisons
#' )
#' }
#'
#' @import data.table
#' @export
testMultipleComparisons <- function(input_dt = "path/to/DIANN_matrix.tsv",
                                    protein_group_annotation = NULL,
                                    study_design = "path/to/Study_design_filled.tsv",
                                    comparisons = NULL,

                                    # toggle normalization & -function
                                    normalize_data = TRUE,
                                    normalization_function = limma::normalizeQuantiles,

                                    # filtering options
                                    min_n_obs = 4,

                                    # imputation of missing values options
                                    imp_percentile = 0.001,
                                    imp_sd = 0.2,

                                    # output options
                                    output_dir = NULL,
                                    plot_pdf = TRUE,
                                    write_tsv_tables = TRUE,

                                    # target protein highlight
                                    target_protein = "O08760",

                                    # error handling
                                    stop_on_error = FALSE,

                                    # progress reporting
                                    verbose = TRUE) {

  #####################################################################################################
  ## INITIALIZATION
  #####################################################################################################

  # Set seed to ensure reproducibility
  set.seed(123)

  # Print analysis header
  if (verbose) {
    print_analysis_header("DiffTestR: Multiple Comparisons Analysis")
  }

  # Load study design to determine comparisons
  if (verbose) {
    report_progress("Initializing multi-comparison analysis", step = 1, total_steps = 4, verbose = verbose)
  }

  if (is.character(study_design)) {
    design_dt <- data.table::fread(study_design)
  } else {
    design_dt <- data.table::as.data.table(study_design)
  }

  # Validate study design
  if (!(all(c("filename", "condition", "replicate") %in% names(design_dt)))) {
    stop("study design must contain columns filename, condition and replicate")
  }

  #####################################################################################################
  ## GENERATE OR VALIDATE COMPARISONS
  #####################################################################################################

  if (is.null(comparisons)) {
    # Auto-generate all pairwise comparisons
    unique_conditions <- unique(design_dt$condition)

    if (length(unique_conditions) < 2) {
      stop("At least 2 unique conditions are required for differential testing")
    }

    if (verbose) {
      report_progress(sprintf("  Auto-generating pairwise comparisons from %d conditions",
                             length(unique_conditions)),
                     type = "info", verbose = verbose)
    }

    # Generate all pairs
    comparisons <- list()
    for (i in 1:(length(unique_conditions) - 1)) {
      for (j in (i + 1):length(unique_conditions)) {
        comparisons[[length(comparisons) + 1]] <- list(
          condition_1 = unique_conditions[i],
          condition_2 = unique_conditions[j]
        )
      }
    }
  } else {
    # Validate user-provided comparisons
    if (!is.list(comparisons)) {
      stop("comparisons must be a list of comparison specifications")
    }

    for (i in seq_along(comparisons)) {
      comp <- comparisons[[i]]
      if (!all(c("condition_1", "condition_2") %in% names(comp))) {
        stop(sprintf("Comparison %d must contain condition_1 and condition_2", i))
      }

      # Check that conditions exist in study design
      available_conditions <- unique(design_dt$condition)
      if (!(comp$condition_1 %in% available_conditions)) {
        stop(sprintf("condition_1 '%s' not found in study design", comp$condition_1))
      }
      if (!(comp$condition_2 %in% available_conditions)) {
        stop(sprintf("condition_2 '%s' not found in study design", comp$condition_2))
      }
    }
  }

  n_comparisons <- length(comparisons)

  if (verbose) {
    report_progress(sprintf("  Total comparisons to perform: %d", n_comparisons),
                   type = "info", verbose = verbose)
    for (i in seq_along(comparisons)) {
      comp <- comparisons[[i]]
      report_progress(sprintf("    %d. %s vs %s", i, comp$condition_1, comp$condition_2),
                     type = "info", verbose = verbose)
    }
  }

  #####################################################################################################
  ## SETUP OUTPUT DIRECTORY
  #####################################################################################################

  if (is.null(output_dir)) {
    output_dir <- getwd()
    if (verbose) {
      report_progress(sprintf("  Using current directory for output: %s", output_dir),
                     type = "info", verbose = verbose)
    }
  } else {
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      if (verbose) {
        report_progress(sprintf("  Created output directory: %s", output_dir),
                       type = "info", verbose = verbose)
      }
    }
  }

  #####################################################################################################
  ## PERFORM COMPARISONS
  #####################################################################################################

  if (verbose) {
    report_progress("Performing differential abundance testing", step = 2, total_steps = 4, verbose = verbose)
  }

  individual_results <- list()
  error_log <- list()
  n_successful <- 0

  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    comparison_name <- sprintf("%s_vs_%s", comp$condition_1, comp$condition_2)

    if (verbose) {
      cat("\n")
      report_progress(sprintf("Comparison %d/%d: %s", i, n_comparisons, comparison_name),
                     type = "step", verbose = verbose)
    }

    # Create subdirectory for this comparison
    comparison_dir <- file.path(output_dir, comparison_name)
    if (!dir.exists(comparison_dir)) {
      dir.create(comparison_dir, recursive = TRUE)
    }

    # Store current working directory
    original_wd <- getwd()

    tryCatch({
      # Change to comparison directory for output
      setwd(comparison_dir)

      # Run differential abundance testing
      result <- testDifferentialAbundance(
        input_dt = input_dt,
        protein_group_annotation = protein_group_annotation,
        study_design = study_design,
        normalize_data = normalize_data,
        normalization_function = normalization_function,
        condition_1 = comp$condition_1,
        condition_2 = comp$condition_2,
        min_n_obs = min_n_obs,
        imp_percentile = imp_percentile,
        imp_sd = imp_sd,
        plot_pdf = plot_pdf,
        write_tsv_tables = write_tsv_tables,
        target_protein = target_protein,
        verbose = verbose
      )

      # Add comparison metadata
      result$comparison_name <- comparison_name
      result$comparison_dir <- comparison_dir

      # Store result
      individual_results[[comparison_name]] <- result
      n_successful <- n_successful + 1

      if (verbose) {
        report_progress(sprintf("  Completed successfully (saved to %s)", comparison_dir),
                       type = "success", verbose = verbose)
      }

    }, error = function(e) {
      # Log error
      error_log[[length(error_log) + 1]] <- list(
        comparison = comparison_name,
        condition_1 = comp$condition_1,
        condition_2 = comp$condition_2,
        error_message = as.character(e$message),
        error_call = as.character(e$call)
      )

      if (verbose) {
        report_progress(sprintf("  FAILED: %s", e$message),
                       type = "error", verbose = verbose)
      }

      if (stop_on_error) {
        # Restore working directory before stopping
        setwd(original_wd)
        stop(sprintf("Error in comparison %s: %s", comparison_name, e$message))
      }

    }, finally = {
      # Always restore working directory
      setwd(original_wd)
    })
  }

  #####################################################################################################
  ## AGGREGATE RESULTS
  #####################################################################################################

  if (verbose) {
    cat("\n")
    report_progress("Aggregating results", step = 3, total_steps = 4, verbose = verbose)
  }

  # Convert error log to data.table
  error_log_dt <- NULL
  if (length(error_log) > 0) {
    error_log_dt <- data.table::rbindlist(error_log)
  }

  # Combine all result data.tables
  combined_results_dt <- NULL
  if (n_successful > 0) {
    result_list <- lapply(individual_results, function(x) {
      dt <- data.table::copy(x$diffExpr_result_dt)
      dt[, comparison := x$comparison_name]
      return(dt)
    })
    combined_results_dt <- data.table::rbindlist(result_list, fill = TRUE)

    if (verbose) {
      report_progress(sprintf("  Combined %d comparisons into single data.table", n_successful),
                     type = "info", verbose = verbose)
      report_progress(sprintf("  Total rows: %d", nrow(combined_results_dt)),
                     type = "info", verbose = verbose)
    }
  }

  #####################################################################################################
  ## FINALIZE AND RETURN RESULTS
  #####################################################################################################

  if (verbose) {
    report_progress("Finalizing multi-comparison results", step = 4, total_steps = 4, verbose = verbose)
  }

  # Prepare overview statistics
  overview_stats <- list(
    n_comparisons = n_comparisons,
    n_successful = n_successful,
    n_failed = n_comparisons - n_successful,
    comparisons_performed = sapply(comparisons, function(x) {
      sprintf("%s vs %s", x$condition_1, x$condition_2)
    })
  )

  # Create result object
  result_list <- list(
    individual_results = individual_results,
    combined_results_dt = combined_results_dt,
    comparisons = comparisons,
    overview_stats = overview_stats,
    n_comparisons = n_comparisons,
    n_successful = n_successful,
    error_log = error_log_dt,
    output_dir = output_dir,
    study_design = design_dt
  )

  # Assign S3 class
  class(result_list) <- c("multiDiffExpr", "list")

  # Print summary
  if (verbose) {
    cat("\n")
    report_progress("Multi-comparison analysis complete!", type = "success", verbose = verbose)
    report_progress(sprintf("  Successfully completed: %d/%d comparisons", n_successful, n_comparisons),
                   type = "info", verbose = verbose)
    if (n_comparisons > n_successful) {
      report_progress(sprintf("  Failed: %d/%d comparisons (see error_log for details)",
                             n_comparisons - n_successful, n_comparisons),
                     type = "warning", verbose = verbose)
    }
    report_progress(sprintf("  Results saved to: %s", output_dir),
                   type = "info", verbose = verbose)
    cat("\n")
  }

  return(result_list)
}
