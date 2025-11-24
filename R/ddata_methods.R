#' S3 Methods for ddata Objects
#'
#' @description S3 class methods for differential expression data (ddata)
#' objects created by testDifferentialAbundance()
#'
#' @name ddata-methods
NULL

#' Print method for ddata objects
#'
#' @param x A ddata object
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.ddata <- function(x, ...) {
  cat("\n")
  cat("====================================\n")
  cat("  DiffTestR ddata Object\n")
  cat("  Differential Abundance Results\n")
  cat("====================================\n\n")

  # Basic information
  cat("Source:", x$data_source, "\n")
  cat("Comparison:", x$comparison, "\n")
  cat(sprintf("  Condition 1: %s\n", x$conditions$condition_1))
  cat(sprintf("  Condition 2: %s\n", x$conditions$condition_2))

  # Extract summary stats if available
  if ("summary_stats" %in% names(x) && !is.null(x$summary_stats)) {
    stats <- x$summary_stats

    cat("\nInput Data:\n")
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
      cat(sprintf("  Higher in %s: %d proteins\n",
                  x$conditions$condition_1, stats$n_upregulated))
      cat(sprintf("  Higher in %s: %d proteins\n",
                  x$conditions$condition_2, stats$n_downregulated))
    }

    if (!is.null(stats$median_log2fc_protein)) {
      cat("\nFold Changes:\n")
      cat(sprintf("  Median log2FC (protein): %.3f\n", stats$median_log2fc_protein))
    }

  } else {
    # Fallback if no summary stats
    if ("diffExpr_result_dt" %in% names(x)) {
      dt <- x$diffExpr_result_dt
      cat(sprintf("\nPrecursors tested: %d\n", nrow(dt)))
      cat(sprintf("Proteins tested: %d\n",
                  length(unique(dt$Protein.Group))))

      n_sig <- sum(dt$p_value_BHadj_protein <= 0.05, na.rm = TRUE)
      cat(sprintf("Significant proteins (p<0.05): %d\n", n_sig))
    }
  }

  cat("\nComponents:\n")
  cat("  $data_source                      - Source file path\n")
  cat("  $comparison                       - Comparison name\n")
  cat("  $conditions                       - Condition names\n")
  cat("  $diffExpr_result_dt               - Full results table\n")
  cat("  $mat_quant_log2_qnorm_imp_minObs  - Processed matrix\n")
  cat("  $mat_quant_log2_qnorm             - Normalized matrix\n")
  cat("  $mat_quant_log2                   - Log2 matrix\n")
  cat("  $mat_quant                        - Raw matrix\n")
  cat("  $study_design                     - Study design\n")
  cat("  $summary_stats                    - Summary statistics\n")
  cat("  $candidates_condition1            - Proteins higher in condition 1\n")
  cat("  $candidates_condition2            - Proteins higher in condition 2\n")

  cat("\n")
  cat("Use summary(x) for more details\n")
  cat("Use plot(x) for volcano plot\n")
  cat("Use plot(x, interactive = TRUE) for interactive volcano plot\n")
  cat("====================================\n\n")

  invisible(x)
}

#' Summary method for ddata objects
#'
#' @param object A ddata object
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns summary statistics
#'
#' @export
summary.ddata <- function(object, ...) {
  cat("\n")
  cat("================================================\n")
  cat("  DiffTestR ddata Summary\n")
  cat("  Differential Abundance Analysis\n")
  cat("================================================\n\n")

  # Basic information
  cat("Data Source:\n")
  cat(sprintf("  %s\n", object$data_source))

  cat("\nComparison:\n")
  cat(sprintf("  %s\n", object$comparison))
  cat(sprintf("  Condition 1: %s\n", object$conditions$condition_1))
  cat(sprintf("  Condition 2: %s\n", object$conditions$condition_2))

  if ("summary_stats" %in% names(object) && !is.null(object$summary_stats)) {
    cat("\n=== Detailed Summary Statistics ===\n\n")

    stats <- object$summary_stats

    cat("Input and Processing:\n")
    cat(sprintf("  Precursors (input):        %d\n", stats$n_precursors_input))
    cat(sprintf("  Precursors (tested):       %d\n", stats$n_precursors_tested))
    cat(sprintf("  Proteins (tested):         %d\n", stats$n_proteins_tested))

    cat("\nSignificance Testing (p < 0.05):\n")
    cat(sprintf("  Significant precursors:    %d (%.1f%%)\n",
                stats$n_significant_precursors,
                100 * stats$n_significant_precursors / stats$n_precursors_tested))
    cat(sprintf("  Significant proteins:      %d (%.1f%%)\n",
                stats$n_significant_proteins,
                100 * stats$n_significant_proteins / stats$n_proteins_tested))

    cat("\nRegulation (|log2FC| > 1, p < 0.05):\n")
    if (!is.null(stats$n_upregulated)) {
      cat(sprintf("  Up-regulated proteins:     %d\n", stats$n_upregulated))
    }
    if (!is.null(stats$n_downregulated)) {
      cat(sprintf("  Down-regulated proteins:   %d\n", stats$n_downregulated))
    }

    cat("\nFold Change Statistics:\n")
    if (!is.null(stats$median_log2fc_precursor)) {
      cat(sprintf("  Median log2FC (precursor): %.3f\n", stats$median_log2fc_precursor))
    }
    if (!is.null(stats$median_log2fc_protein)) {
      cat(sprintf("  Median log2FC (protein):   %.3f\n", stats$median_log2fc_protein))
    }

    if (!is.null(stats$median_cv_condition1) && !is.null(stats$median_cv_condition2)) {
      cat("\nData Quality (Coefficient of Variation):\n")
      cat(sprintf("  Median CV %s:  %.1f%%\n",
                  object$conditions$condition_1, stats$median_cv_condition1))
      cat(sprintf("  Median CV %s:  %.1f%%\n",
                  object$conditions$condition_2, stats$median_cv_condition2))
    }
  }

  # Candidate proteins
  cat("\nCandidate Proteins:\n")
  if ("candidates_condition1" %in% names(object) && nrow(object$candidates_condition1) > 0) {
    cat(sprintf("  Higher in %s: %d proteins\n",
                object$conditions$condition_1,
                length(unique(object$candidates_condition1$Protein.Group))))

    # Show top 5
    top5 <- unique(object$candidates_condition1[order(-log2_fold_change_protein)][1:min(5, .N)])
    if (nrow(top5) > 0) {
      cat("\n  Top 5 (by log2FC):\n")
      for (i in 1:nrow(top5)) {
        cat(sprintf("    %d. %s (log2FC: %.2f, p: %.2e)\n",
                    i,
                    top5$Protein.Names[i],
                    top5$log2_fold_change_protein[i],
                    top5$p_value_BHadj_protein[i]))
      }
    }
  } else {
    cat(sprintf("  Higher in %s: 0 proteins\n", object$conditions$condition_1))
  }

  cat("\n")
  if ("candidates_condition2" %in% names(object) && nrow(object$candidates_condition2) > 0) {
    cat(sprintf("  Higher in %s: %d proteins\n",
                object$conditions$condition_2,
                length(unique(object$candidates_condition2$Protein.Group))))

    # Show top 5
    top5 <- unique(object$candidates_condition2[order(log2_fold_change_protein)][1:min(5, .N)])
    if (nrow(top5) > 0) {
      cat("\n  Top 5 (by log2FC):\n")
      for (i in 1:nrow(top5)) {
        cat(sprintf("    %d. %s (log2FC: %.2f, p: %.2e)\n",
                    i,
                    top5$Protein.Names[i],
                    top5$log2_fold_change_protein[i],
                    top5$p_value_BHadj_protein[i]))
      }
    }
  } else {
    cat(sprintf("  Higher in %s: 0 proteins\n", object$conditions$condition_2))
  }

  cat("\n================================================\n\n")

  invisible(object$summary_stats)
}

#' Plot method for ddata objects
#'
#' @param x A ddata object
#' @param interactive Logical, whether to create an interactive plotly plot
#'   (default: FALSE for static ggplot2)
#' @param which Type of plot: "volcano" (default) or "ma"
#' @param p_cutoff P-value cutoff for significance (default: 0.05)
#' @param fc_cutoff Log2 fold change cutoff for significance (default: 1)
#' @param label_top Number of top proteins to label (default: 10, set to 0 for no labels)
#' @param point_size Size of points in the plot (default: 2)
#' @param alpha Transparency of points (default: 0.6)
#' @param ... Additional arguments passed to ggplot2 or plotly
#'
#' @return A ggplot2 object (if interactive=FALSE) or plotly object (if interactive=TRUE)
#'
#' @examples
#' \dontrun{
#' # Static volcano plot
#' plot(ddata)
#'
#' # Interactive volcano plot
#' plot(ddata, interactive = TRUE)
#'
#' # MA plot
#' plot(ddata, which = "ma")
#'
#' # Custom thresholds
#' plot(ddata, p_cutoff = 0.01, fc_cutoff = 2)
#' }
#'
#' @import ggplot2 data.table
#' @importFrom plotly ggplotly
#' @importFrom ggrepel geom_text_repel
#'
#' @export
plot.ddata <- function(x,
                       interactive = FALSE,
                       which = "volcano",
                       p_cutoff = 0.05,
                       fc_cutoff = 1,
                       label_top = 10,
                       point_size = 2,
                       alpha = 0.6,
                       ...) {

  # Validate ddata object
  validate_ddata(x)

  # Get the result data
  dt <- as.data.table(x$diffExpr_result_dt)

  # Get unique protein-level data
  plot_data <- unique(dt[, .(
    Protein.Group,
    Protein.Names,
    log2_fold_change_protein,
    p_value_BHadj_protein,
    n_precursors
  )])

  # Calculate -log10(p-value) for volcano plot
  plot_data[, neg_log10_p := -log10(p_value_BHadj_protein)]

  # Handle Inf values (p-value = 0)
  max_finite_p <- max(plot_data[is.finite(neg_log10_p), neg_log10_p], na.rm = TRUE)
  plot_data[is.infinite(neg_log10_p), neg_log10_p := max_finite_p * 1.2]

  # Determine significance category
  plot_data[, significance := "Not Significant"]
  plot_data[p_value_BHadj_protein <= p_cutoff & log2_fold_change_protein >= fc_cutoff,
            significance := sprintf("Up in %s", x$conditions$condition_1)]
  plot_data[p_value_BHadj_protein <= p_cutoff & log2_fold_change_protein <= -fc_cutoff,
            significance := sprintf("Up in %s", x$conditions$condition_2)]
  plot_data[p_value_BHadj_protein <= p_cutoff &
              abs(log2_fold_change_protein) < fc_cutoff,
            significance := "Significant (|log2FC| < threshold)"]

  # Define colors
  up_in_cond1 <- sprintf("Up in %s", x$conditions$condition_1)
  up_in_cond2 <- sprintf("Up in %s", x$conditions$condition_2)

  colors <- c(
    "Not Significant" = "grey70",
    "Significant (|log2FC| < threshold)" = "grey40"
  )
  colors[up_in_cond1] <- "#E41A1C"
  colors[up_in_cond2] <- "#377EB8"

  if (which == "volcano") {
    # Create volcano plot
    p <- ggplot(plot_data, aes(
      x = log2_fold_change_protein,
      y = neg_log10_p,
      color = significance,
      text = paste0(
        "Protein: ", Protein.Names, "\n",
        "log2FC: ", round(log2_fold_change_protein, 3), "\n",
        "p-value: ", format(p_value_BHadj_protein, digits = 3, scientific = TRUE), "\n",
        "Precursors: ", n_precursors
      )
    )) +
      geom_point(size = point_size, alpha = alpha) +
      geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
      geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      scale_color_manual(values = colors, name = "Regulation") +
      labs(
        title = sprintf("Volcano Plot: %s", x$comparison),
        subtitle = sprintf("Cutoffs: p < %.3f, |log2FC| > %g", p_cutoff, fc_cutoff),
        x = sprintf("log2 Fold Change (%s / %s)",
                    x$conditions$condition_1, x$conditions$condition_2),
        y = "-log10(Adjusted p-value)"
      ) +
      theme_bw() +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()
      )

    # Add labels for top proteins if requested
    if (label_top > 0) {
      # Select top proteins by significance (lowest p-value and highest |log2FC|)
      sig_proteins <- plot_data[significance %in% c(
        sprintf("Up in %s", x$conditions$condition_1),
        sprintf("Up in %s", x$conditions$condition_2)
      )]

      if (nrow(sig_proteins) > 0) {
        # Rank by p-value and fold change
        sig_proteins[, rank_score := neg_log10_p * abs(log2_fold_change_protein)]
        top_proteins <- sig_proteins[order(-rank_score)][1:min(label_top, .N)]

        p <- p + ggrepel::geom_text_repel(
          data = top_proteins,
          aes(label = Protein.Names),
          size = 3,
          max.overlaps = 20,
          box.padding = 0.5,
          point.padding = 0.3,
          segment.color = "grey50",
          segment.size = 0.3,
          show.legend = FALSE
        )
      }
    }

  } else if (which == "ma") {
    # Calculate mean abundance for MA plot
    # Get mean abundance from the normalized matrix
    mean_abundance <- rowMeans(x$mat_quant_log2_qnorm_imp_minObs, na.rm = TRUE)
    ma_data <- data.table(
      Precursor.Id = rownames(x$mat_quant_log2_qnorm_imp_minObs),
      mean_abundance = mean_abundance
    )

    # Merge with plot_data via diffExpr_result_dt
    ma_data <- merge(ma_data, dt[, .(Precursor.Id, Protein.Group)], by = "Precursor.Id")
    ma_data <- ma_data[, .(mean_abundance = mean(mean_abundance)), by = Protein.Group]
    plot_data <- merge(plot_data, ma_data, by = "Protein.Group")

    # Create MA plot
    p <- ggplot(plot_data, aes(
      x = mean_abundance,
      y = log2_fold_change_protein,
      color = significance,
      text = paste0(
        "Protein: ", Protein.Names, "\n",
        "log2FC: ", round(log2_fold_change_protein, 3), "\n",
        "p-value: ", format(p_value_BHadj_protein, digits = 3, scientific = TRUE), "\n",
        "Mean abundance: ", round(mean_abundance, 2)
      )
    )) +
      geom_point(size = point_size, alpha = alpha) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
      scale_color_manual(values = colors, name = "Regulation") +
      labs(
        title = sprintf("MA Plot: %s", x$comparison),
        subtitle = sprintf("Cutoffs: p < %.3f, |log2FC| > %g", p_cutoff, fc_cutoff),
        x = "Mean log2 Abundance",
        y = sprintf("log2 Fold Change (%s / %s)",
                    x$conditions$condition_1, x$conditions$condition_2)
      ) +
      theme_bw() +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()
      )

    # Add labels for top proteins if requested
    if (label_top > 0) {
      sig_proteins <- plot_data[significance %in% c(
        sprintf("Up in %s", x$conditions$condition_1),
        sprintf("Up in %s", x$conditions$condition_2)
      )]

      if (nrow(sig_proteins) > 0) {
        sig_proteins[, rank_score := -log10(p_value_BHadj_protein) * abs(log2_fold_change_protein)]
        top_proteins <- sig_proteins[order(-rank_score)][1:min(label_top, .N)]

        p <- p + ggrepel::geom_text_repel(
          data = top_proteins,
          aes(label = Protein.Names),
          size = 3,
          max.overlaps = 20,
          box.padding = 0.5,
          point.padding = 0.3,
          segment.color = "grey50",
          segment.size = 0.3,
          show.legend = FALSE
        )
      }
    }

  } else {
    stop("Unknown plot type. Use 'volcano' or 'ma'.")
  }

  # Convert to interactive plot if requested
  if (interactive) {
    p <- plotly::ggplotly(p, tooltip = "text")
  }

  return(p)
}

#' Validate ddata object structure
#'
#' @param x A ddata object
#'
#' @return TRUE if valid, stops with error if invalid
#'
#' @export
validate_ddata <- function(x) {
  if (!inherits(x, "ddata")) {
    stop("Object must be of class 'ddata'")
  }

  required_components <- c(
    "data_source",
    "comparison",
    "conditions",
    "diffExpr_result_dt",
    "mat_quant_log2_qnorm_imp_minObs",
    "mat_quant_log2_qnorm",
    "mat_quant_log2",
    "mat_quant",
    "study_design"
  )

  missing <- setdiff(required_components, names(x))
  if (length(missing) > 0) {
    stop("Missing required components: ", paste(missing, collapse = ", "))
  }

  # Check conditions structure
  if (!is.list(x$conditions)) {
    stop("'conditions' must be a list")
  }
  if (!all(c("condition_1", "condition_2") %in% names(x$conditions))) {
    stop("'conditions' must contain 'condition_1' and 'condition_2'")
  }

  # Check diffExpr_result_dt has required columns
  required_cols <- c(
    "Protein.Group",
    "Precursor.Id",
    "log2_fold_change_protein",
    "p_value_BHadj_protein"
  )
  missing_cols <- setdiff(required_cols, names(x$diffExpr_result_dt))
  if (length(missing_cols) > 0) {
    stop("diffExpr_result_dt missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  return(TRUE)
}
