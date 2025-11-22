#' Progress Reporting Utilities for DiffTestR
#'
#' @description Internal utility functions for progress reporting during analysis
#'
#' @keywords internal
#' @noRd
NULL

#' Report progress message
#'
#' @param msg Message to display
#' @param step Current step number
#' @param total_steps Total number of steps
#' @param verbose Whether to display message
#' @param type Type of message: "step", "info", "success", "warning", "error"
#'
#' @keywords internal
#' @noRd
report_progress <- function(msg, step = NULL, total_steps = NULL,
                           verbose = TRUE, type = "step") {
  if (!verbose) return(invisible(NULL))

  # Build prefix
  prefix <- switch(type,
    "step" = {
      if (!is.null(step) && !is.null(total_steps)) {
        sprintf("[%d/%d]", step, total_steps)
      } else {
        "[INFO]"
      }
    },
    "info" = "[INFO]",
    "success" = "[OK]",
    "warning" = "[WARNING]",
    "error" = "[ERROR]",
    ""
  )

  # Format message
  full_msg <- if (nchar(prefix) > 0) {
    sprintf("%s %s\n", prefix, msg)
  } else {
    sprintf("%s\n", msg)
  }

  # Color support (basic)
  if (type == "success") {
    cat(full_msg)  # Could add color codes if desired
  } else if (type == "warning") {
    cat(full_msg)
  } else if (type == "error") {
    cat(full_msg)
  } else {
    cat(full_msg)
  }
}

#' Create a simple text-based progress bar
#'
#' @param current Current iteration
#' @param total Total iterations
#' @param width Width of progress bar in characters
#' @param prefix Optional prefix message
#' @param verbose Whether to display
#'
#' @keywords internal
#' @noRd
show_progress_bar <- function(current, total, width = 50,
                              prefix = "", verbose = TRUE) {
  if (!verbose) return(invisible(NULL))

  percent <- current / total
  filled <- round(width * percent)
  empty <- width - filled

  bar <- paste0(
    "\r", prefix,
    " [", paste(rep("=", filled), collapse = ""),
    ifelse(filled < width, ">", ""),
    paste(rep(" ", max(0, empty - 1)), collapse = ""),
    "] ",
    sprintf("%3.0f%%", percent * 100)
  )

  cat(bar)
  if (current == total) cat("\n")
  flush.console()
}

#' Print analysis header
#'
#' @param title Analysis title
#' @param verbose Whether to display
#'
#' @keywords internal
#' @noRd
print_analysis_header <- function(title = "DiffTestR Analysis Pipeline",
                                  verbose = TRUE) {
  if (!verbose) return(invisible(NULL))

  width <- max(nchar(title) + 4, 60)
  border <- paste(rep("=", width), collapse = "")

  cat("\n", border, "\n", sep = "")
  cat(" ", title, "\n", sep = "")
  cat(border, "\n\n", sep = "")
}

#' Print analysis summary section
#'
#' @param title Section title
#' @param verbose Whether to display
#'
#' @keywords internal
#' @noRd
print_section_header <- function(title, verbose = TRUE) {
  if (!verbose) return(invisible(NULL))
  cat("\n--- ", title, " ---\n", sep = "")
}
