#' S3 Methods for qdata Objects
#'
#' @description S3 class methods for proteomics quantification data (qdata)
#' objects created by importFromDIANN() and importFromTables()
#'
#' @name qdata-methods
NULL

#' Print method for qdata objects
#'
#' @param x A qdata object
#' @param ... Additional arguments (not used)
#'
#' @export
print.qdata <- function(x, ...) {
  cat("\n=== DiffTestR qdata Object ===\n\n")

  cat("Source file:", x$source_file, "\n")

  cat("\nDimensions:\n")
  cat(sprintf("  Precursors: %d\n", nrow(x$data_wide)))
  cat(sprintf("  Samples: %d\n", ncol(x$matrix_raw)))

  cat("\nStudy Design:\n")
  if (is.data.frame(x$study_design) || is.data.table(x$study_design)) {
    print(x$study_design)
  } else {
    cat("  Study design file:", x$study_design, "\n")
  }

  cat("\nConditions:\n")
  if (is.data.frame(x$ann_col) || is.data.table(x$ann_col)) {
    cond_table <- table(x$ann_col$condition)
    print(cond_table)
  }

  cat("\nAnnotation columns available:\n")
  cat(" ", paste(names(x$ann_row), collapse = ", "), "\n")

  cat("\nComponents:\n")
  cat("  $data_long   - Long format data\n")
  cat("  $data_wide   - Wide format data\n")
  cat("  $matrix_raw  - Raw intensity matrix\n")
  cat("  $ann_col     - Sample annotations\n")
  cat("  $ann_row     - Precursor annotations\n")
  cat("  $study_design - Study design\n")
  cat("  $source_file - Source file path\n")

  cat("\n")
  invisible(x)
}

#' Summary method for qdata objects
#'
#' @param object A qdata object
#' @param ... Additional arguments (not used)
#'
#' @export
summary.qdata <- function(object, ...) {
  cat("\n=== DiffTestR qdata Summary ===\n\n")

  cat("Data Dimensions:\n")
  cat(sprintf("  Precursors: %d\n", nrow(object$data_wide)))
  cat(sprintf("  Proteins: %d\n", length(unique(object$ann_row$Protein.Group))))
  cat(sprintf("  Samples: %d\n", ncol(object$matrix_raw)))

  cat("\nMissing Values:\n")
  n_missing <- sum(is.na(object$matrix_raw))
  total <- length(object$matrix_raw)
  pct_missing <- 100 * n_missing / total
  cat(sprintf("  Total: %d / %d (%.1f%%)\n", n_missing, total, pct_missing))

  cat("\nIntensity Range (raw):\n")
  cat(sprintf("  Min: %.2e\n", min(object$matrix_raw, na.rm = TRUE)))
  cat(sprintf("  Max: %.2e\n", max(object$matrix_raw, na.rm = TRUE)))
  cat(sprintf("  Median: %.2e\n", median(object$matrix_raw, na.rm = TRUE)))

  cat("\nConditions:\n")
  if (is.data.frame(object$ann_col) || is.data.table(object$ann_col)) {
    cond_counts <- as.data.table(table(object$ann_col$condition))
    setnames(cond_counts, c("Condition", "N_Samples"))
    print(cond_counts)
  }

  cat("\n")
  invisible(object)
}

#' Validate qdata object structure
#'
#' @param x A qdata object
#'
#' @return TRUE if valid, stops with error if invalid
#'
#' @export
validate_qdata <- function(x) {
  if (!inherits(x, "qdata")) {
    stop("Object must be of class 'qdata'")
  }

  required_components <- c("data_long", "data_wide", "matrix_raw",
                           "ann_col", "ann_row", "study_design",
                           "source_file")

  missing <- setdiff(required_components, names(x))
  if (length(missing) > 0) {
    stop("Missing required components: ", paste(missing, collapse = ", "))
  }

  # Check dimensions match
  if (nrow(x$data_wide) != nrow(x$matrix_raw)) {
    stop("Dimension mismatch: data_wide rows != matrix_raw rows")
  }

  if (ncol(x$matrix_raw) != nrow(x$ann_col)) {
    stop("Dimension mismatch: matrix_raw cols != ann_col rows")
  }

  if (nrow(x$data_wide) != nrow(x$ann_row)) {
    stop("Dimension mismatch: data_wide rows != ann_row rows")
  }

  return(TRUE)
}
