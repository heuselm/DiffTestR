#' Test Differential Protein Expression in MS proteomics data
#'
#' @description Test Differential Protein Expression in MS proteomics data
#' starting from the precursor level. This function now primarily accepts qdata
#' objects from importFromDIANN() or importFromTables().
#'
#' @param data Either a qdata object (recommended) OR a path to a data file
#' (will attempt to import). If a data.table/data.frame is provided, a
#' deprecation warning is shown.
#' @param normalize_data Whether or not data is scaled/normalized before
#' differential testing. In some cases it might be preferable not to scale
#' the datasets, e.g. when comparing pulldowns vs. input samples!
#' Defaults to TRUE.
#' @param normalization_function Normalization function to use that transforms
#' a matrix of quantities where columns are samples and rows are analytes.
#' Defaults to limma::normalizeQuantiles, but can be replaced with any such
#' function. Try limma::normalizeVSN or limma::normalizeMedianValues.
#' @param condition_1 First condition for the differential comparison. By
#' default it is guessed from unique conditions in the qdata object.
#' @param condition_2 Second condition for the differential comparison. By
#' default it is guessed from unique conditions in the qdata object.
#' @param min_n_obs Minimum number of observations per precursor (number of
#' runs it was identified in) to keep it in the analysis (default: 4)
#' @param imp_percentile Percentile of the total distribution of values on
#' which the random distribution for sampling will be centered (default: 0.001)
#' @param imp_sd Standard deviation of the normal distribution from which
#' values are sampled to impute missing values (default: 0.2)
#' @param plot_pdf Document processing steps in a series of pdf graphs
#' (default: TRUE)
#' @param write_tsv_tables Write out final quant table with differential
#' expression testing results (default: TRUE)
#' @param target_protein Optional string with protein identifier to highlight
#' in volcano plots
#' @param verbose Logical, whether to print progress messages during analysis
#' (default: TRUE)
#'
#' @return A ddata object (S3 class) containing differential expression results:
#' \itemize{
#' \item data_source: Source file path
#' \item comparison: String describing the comparison (e.g., "A/B")
#' \item conditions: Named list with condition_1 and condition_2
#' \item diffExpr_result_dt: Result table with intensities and statistics
#' \item mat_quant_log2_qnorm_imp_minObs: Processed matrix (log2, normalized,
#'   imputed, filtered)
#' \item mat_quant_log2_qnorm: Log2 normalized matrix
#' \item mat_quant_log2: Log2 transformed matrix
#' \item mat_quant: Raw intensity matrix
#' \item study_design: Study design table
#' \item candidates_condition1: Proteins higher in condition 1
#' \item candidates_condition2: Proteins higher in condition 2
#' \item summary_stats: Summary statistics of the analysis
#' }
#'
#' @details
#' This function performs the following steps:
#' 1. Log2 transformation
#' 2. Optional normalization (quantile normalization by default)
#' 3. Filtering based on minimum observations
#' 4. Missing value imputation
#' 5. Precursor-level t-tests
#' 6. Protein-level rollup and multiple testing correction
#'
#' The returned ddata object has plot() and summary() methods for
#' visualization and result exploration.
#'
#' @examples
#' \dontrun{
#' # Recommended workflow
#' qdata <- importFromTables("data.tsv", study_design = "design.txt")
#' ddata <- testDifferentialAbundance(qdata)
#'
#' # Interactive volcano plot
#' plot(ddata, interactive = TRUE)
#'
#' # Static PDF volcano plot
#' plot(ddata, interactive = FALSE)
#'
#' # Summary
#' summary(ddata)
#' }
#'
#' @author Moritz Heusel
#'
#' @import data.table ggplot2 corrplot pheatmap ggrepel limma
#' @importFrom reshape2 melt
#'
#' @export

testDifferentialAbundance <- function(data,
                                      # toggle normalization & -function
                                      normalize_data = TRUE,
                                      normalization_function = limma::normalizeQuantiles,

                                      # select conditions to be compared
                                      condition_1 = NULL,
                                      condition_2 = NULL,

                                      # filtering options
                                      min_n_obs = 4,

                                      # imputation of missing values options
                                      imp_percentile = 0.001,
                                      imp_sd = 0.2,

                                      # output options
                                      plot_pdf = TRUE,
                                      write_tsv_tables = TRUE,

                                      # target protein highlight
                                      target_protein = NULL,

                                      # progress reporting
                                      verbose = TRUE)
{
  #####################################################################################################

  ## Set seed to ensure reproducibility
  set.seed(123)

  # Print analysis header
  if (verbose) {
    print_analysis_header("DiffTestR: Differential Abundance Analysis")
  }

  ## Handle different input types
  if (verbose) {
    report_progress("Preparing data", step = 1, total_steps = 8, verbose = verbose)
  }

  # Check if qdata object
  if (inherits(data, "qdata")) {
    # Extract from qdata object
    if (verbose) {
      report_progress("  Using qdata object", type = "info", verbose = verbose)
    }

    data_source <- data$source_file
    data.s.wide <- data$data_wide
    annotation <- as.data.table(data$ann_col)
    protein_group_annotation <- unique(data$ann_row[, .(Protein.Group, Protein.Names)])
    input_n_precursors <- nrow(data$data_wide)

    # Prepare matrix and column names
    data.s.wide.quant <- data$matrix_raw
    annotation_col <- as.data.frame(annotation)
    row.names(annotation_col) <- colnames(data.s.wide.quant)

  } else if (is.character(data)) {
    # Try to import the file
    warning("Providing a file path directly is deprecated. Please use importFromTables() or importFromDIANN() first, then pass the qdata object.")

    if (verbose) {
      report_progress("  Attempting to import file with importFromTables()",
                     type = "warning", verbose = verbose)
    }

    # Try to find study_design file in same directory
    study_design_path <- file.path(dirname(data), "study_design.txt")
    if (!file.exists(study_design_path)) {
      study_design_path <- file.path(dirname(data), "Study_design.txt")
    }
    if (!file.exists(study_design_path)) {
      stop("Cannot find study_design.txt in the same directory as the data file. Please import data using importFromTables() or importFromDIANN() first.")
    }

    # Import and recursively call this function
    qdata_obj <- importFromTables(path_to_tsv = data,
                                   study_design = study_design_path)
    return(testDifferentialAbundance(
      data = qdata_obj,
      normalize_data = normalize_data,
      normalization_function = normalization_function,
      condition_1 = condition_1,
      condition_2 = condition_2,
      min_n_obs = min_n_obs,
      imp_percentile = imp_percentile,
      imp_sd = imp_sd,
      plot_pdf = plot_pdf,
      write_tsv_tables = write_tsv_tables,
      target_protein = target_protein,
      verbose = verbose
    ))

  } else {
    # Legacy mode - data.table or data.frame
    stop("Direct data.table/data.frame input is no longer supported. Please use importFromTables() or importFromDIANN() to create a qdata object first.")
  }

  # Auto-detect conditions if not provided
  if (is.null(condition_1)) {
    condition_1 <- unique(annotation$condition)[1]
    if (verbose) {
      report_progress(sprintf("  Auto-detected condition_1: %s", condition_1),
                     type = "info", verbose = verbose)
    }
  }

  if (is.null(condition_2)) {
    condition_2 <- unique(annotation$condition)[2]
    if (verbose) {
      report_progress(sprintf("  Auto-detected condition_2: %s", condition_2),
                     type = "info", verbose = verbose)
    }
  }

  # log transform and remove inf
  if (verbose) {
    report_progress("Log2 transformation", step = 2, total_steps = 8, verbose = verbose)
  }

  data.s.wide.quant.log2 = log2(data.s.wide.quant)
  data.s.wide.quant.log2[is.infinite(data.s.wide.quant.log2)] <- 0

  # Plot abundance distributions before and after normalization (if applied)
  boxplot(data.s.wide.quant.log2, las = 2, main = "Input intensities, log2-transformed")

  if (normalize_data == FALSE) {
    normalization_function = function(x){return(x)}
    if (verbose) {
      report_progress("Normalization: SKIPPED (normalize_data = FALSE)", step = 3,
                     total_steps = 8, verbose = verbose)
    }
  } else {
    if (verbose) {
      report_progress("Quantile normalization", step = 3, total_steps = 8,
                     verbose = verbose)
    }
  }

  data.s.wide.quant.log2.qnorm = normalization_function(data.s.wide.quant.log2)
  boxplot(data.s.wide.quant.log2.qnorm, las = 2, main = "log2-transformed and quantile-normalized intensities")
  data.s.wide.log2.qnorm = cbind(data.s.wide[,1:2],data.s.wide.quant.log2.qnorm)

  # Plot normalization effect in density plot
  data.s.long.prenorm = reshape2::melt(data.s.wide.quant.log2)
  data.s.long.prenorm$normalization = "before (input)"
  data.s.long.qNorm = reshape2::melt(data.s.wide.quant.log2.qnorm)
  data.s.long.qNorm$normalization = "after normalization"
  data.s.long = rbind(data.s.long.prenorm, data.s.long.qNorm)

  ggplot(data.s.long, aes(value, group = Var2, color = normalization)) + geom_density() + facet_wrap(~normalization, ncol = 1) +
    ggtitle(paste("Impact of normalization, normalize_data = ",normalize_data)) + theme_bw() + xlab("log2 Precursor.Quantity")
  if(plot_pdf){
    ggsave("Normalization_impact_abundance_density.pdf", height = 5, width = 5)
  }

  ## Impute missing values
  # Visualize input data before imputation
  data.s.wide.quant.log2.qnorm.noNa = copy(data.s.wide.quant.log2.qnorm)
  data.s.wide.quant.log2.qnorm.noNa[is.na(data.s.wide.quant.log2.qnorm.noNa)] = 0
  pheatmap::pheatmap(data.s.wide.quant.log2.qnorm.noNa[order(-rowSums(data.s.wide.quant.log2.qnorm.noNa)),], show_rownames = F,
           main = "Quantile-normalized log2 precursor-level intensities, non-imputed", cluster_rows = F)
  if(plot_pdf){
    dev.copy(pdf, "Heatmap_unfiltered_pre_imputation.pdf", width = ncol(data.s.wide.quant.log2.qnorm.noNa)/2, height = 8)
    dev.off()
  }
  # too big for hclust, thus cluster_rows = F

  # Check correlation to remove outliers / select data subset
  dev.off()
  corrplot::corrplot(cor(data.s.wide.quant.log2.qnorm.noNa),
           method = "color",
           is.corr = F,
           order = "hclust",
           mar = c(4,4,4,4),
           tl.cex = 1-(ncol(data.s.wide.quant.log2.qnorm.noNa)/100),
           title = "Sample-sample correlation before min_obs filtering")
  if(plot_pdf){
    dev.copy(pdf, "Correlation_sample_sample_unfiltered.pdf")
    dev.off()
  }

  # Filter dataset based on minimum n_obs to avoid imputation around noise..
  if (verbose) {
    report_progress(sprintf("Filtering (min_n_obs >= %d)", min_n_obs), step = 4,
                   total_steps = 8, verbose = verbose)
  }

  n_obs = apply(data.s.wide.quant.log2.qnorm, MARGIN = 1, function(x) length(x)-sum(is.na(x)))
  hist(n_obs, n = 100)
  abline(v = min_n_obs, col = "red", lty = 2)
  data.s.wide.quant.log2.qnorm.noNa.minObs = data.s.wide.quant.log2.qnorm.noNa[n_obs >= min_n_obs,]

  if (verbose) {
    report_progress(sprintf("  Retained %d/%d precursors after filtering",
                           nrow(data.s.wide.quant.log2.qnorm.noNa.minObs),
                           nrow(data.s.wide.quant.log2.qnorm)),
                   type = "info", verbose = verbose)
  }

  nrow(data.s.wide.quant.log2.qnorm.noNa.minObs)

  plot(n_obs[order(-rowSums(data.s.wide.quant.log2.qnorm.noNa))], main = "n_observations over abundance rank")
  abline(h = min_n_obs, col = "red", lty = 2)

  if(plot_pdf){
    par(mfrow = c(2,1))
    pdf("Filtering_min_obs.pdf")
    hist(n_obs, n = 100)
    abline(v = min_n_obs, col = "red", lty = 2)
    plot(n_obs[order(-rowSums(data.s.wide.quant.log2.qnorm.noNa))],
         main = "n_observations over abundance rank")
    abline(h = min_n_obs, col = "red", lty = 2)
    dev.off()
  }

  # Plot heatmap after min obs filtering, before imputation
  pheatmap::pheatmap(data.s.wide.quant.log2.qnorm.noNa.minObs[order(-rowSums(data.s.wide.quant.log2.qnorm.noNa.minObs)),], show_rownames = F,
           ann_col = annotation_col,
           main = "Quantile-normalized log2 precursor-level intensities, min obs. filtered, before imputation",
           cluster_rows = F,
           cluster_cols = F)
  if(plot_pdf){
    dev.copy(pdf, "Heatmap_min_obs_filtered_pre_imputation.pdf", width = ncol(data.s.wide.quant.log2.qnorm.noNa)/2, height = 8)
    dev.off()
  }

  # Sample-sample correlation after filtering to min obs
  dev.off()
  corrplot::corrplot(cor(data.s.wide.quant.log2.qnorm.noNa.minObs),
           method = "color",
           is.corr = F,
           order = "hclust",
           mar = c(4,4,4,4),
           tl.cex = 1-(ncol(data.s.wide.quant.log2.qnorm.noNa)/100),
           title = "Sample-sample correlation after min_obs filtering")
  if(plot_pdf){
    dev.copy(pdf, "Correlation_sample_sample_min_obs_filtered_pre_imputation.pdf")
    dev.off()
  }

    ## Impute missing values
  if (verbose) {
    report_progress("Imputing missing values", step = 5, total_steps = 8,
                   verbose = verbose)
  }

  # distribution of all values before imputation
  hist(data.s.wide.quant.log2.qnorm.noNa.minObs)

  # impute missing values
  imp_percentile = quantile(unique(data.s.wide.quant.log2.qnorm.noNa.minObs[data.s.wide.quant.log2.qnorm.noNa.minObs>0]),
                        na.rm = T, imp_percentile)
  data.s.wide.quant.log2.qnorm.noNa.minObs.imp = copy(data.s.wide.quant.log2.qnorm.noNa.minObs)

  n_imputed <- sum(data.s.wide.quant.log2.qnorm.noNa.minObs.imp == 0)

  data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs.imp == 0] =
    rnorm(length(data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs.imp == 0]),
          mean = imp_percentile,  sd = imp_sd)

  if (verbose) {
    report_progress(sprintf("  Imputed %d missing values", n_imputed),
                   type = "info", verbose = verbose)
  }

  # visualize distribution of imputed values
  hist(data.s.wide.quant.log2.qnorm.noNa.minObs)
  hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs == 0], add = T, col = "red")
  hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp)
  hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs == 0], add = T, col = "red")

  if(plot_pdf){
    par(mfrow = c(2,1))
    pdf("Imputation_histograms.pdf")
    hist(data.s.wide.quant.log2.qnorm.noNa.minObs)
    hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs == 0], add = T, col = "red")
    hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp)
    hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs == 0], add = T, col = "red")
    dev.off()
  }

  # heatmap after filtering to min obs and imputation
  visualize_imputation_heatmap = cbind(data.s.wide.quant.log2.qnorm.noNa.minObs, data.s.wide.quant.log2.qnorm.noNa.minObs.imp)
  colnames(visualize_imputation_heatmap)[(ncol(data.s.wide.quant.log2.qnorm.noNa.minObs)+1):ncol(visualize_imputation_heatmap)] =
    paste(colnames(data.s.wide.quant.log2.qnorm.noNa.minObs.imp), "imputed")

  # To avoid ram overflow the plottable number of precursors is limited to 20000
  if(nrow(visualize_imputation_heatmap)>20000){
    visualize_imputation_heatmap =
      visualize_imputation_heatmap[sample(1:nrow(visualize_imputation_heatmap), 20000, replace = F),]
  }

  pheatmap::pheatmap(visualize_imputation_heatmap[order(-rowSums(visualize_imputation_heatmap)),],
           main = "Quantile-normalized log2 precursor-level intensities, min obs.\nnon-imputed      vs       imputed",
           show_rownames = F,
           gaps_col = ncol(data.s.wide.quant.log2.qnorm.noNa.minObs.imp),
           cluster_rows = F,
           cluster_cols = F)
  dev.copy(pdf, "Imputation_heatmap.pdf", width = ncol(data.s.wide.quant.log2.qnorm.noNa.minObs.imp))
  dev.off()

  ## Precursor level differential Abundance testing t.tests

  # Assemble differential expression testing result table
  res = copy(data.s.wide[,1:2])

  # Add raw data to the result table
  data.s.wide.quant.dt = as.data.table(data.s.wide.quant, keep.rownames = T)
  names(data.s.wide.quant.dt)[-1] = paste(names(data.s.wide.quant.dt)[-1], "raw_intensities")
  res = merge(res, data.s.wide.quant.dt,
              by.x = "Precursor.Id", by.y = "rn", all.x = F)

  # Add the filtered, normalized, imputed quantdata
  data.s.wide.quant.log2.qnorm.noNa.minObs.imp.dt = as.data.table(data.s.wide.quant.log2.qnorm.noNa.minObs.imp, keep.rownames = T)
  names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp.dt)[-1] = paste(names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp.dt)[-1], "log2_qnorm_imp_intensities")

  res = merge(res, data.s.wide.quant.log2.qnorm.noNa.minObs.imp.dt,
              by.x = "Precursor.Id", by.y = "rn", all.x = F)

  ## perform t.tests

  # get conditions
  if (is.null(condition_1)){
    condition_1 = unique(annotation$condition[1])
  }
  condition_1
  if (is.null(condition_1)){
    condition_2 = unique(annotation$condition[2])
  }
  condition_2

  # Keep only precursors that made it through the filtering
  res = res[Precursor.Id %in% row.names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp)]
  res[, comparison:=paste0(condition_1,"/",condition_2)]

  # preallocate result vectors
  prec = res$Precursor.Id
  n_prec = length(prec)
  pvals = numeric(length = n_prec)
  log2fcs = numeric(length = n_prec)

  # presplit matrix
  data.s.wide.quant.log2.qnorm.noNa.minObs.imp =
    data.s.wide.quant.log2.qnorm.noNa.minObs.imp[order(row.names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp), method = "radix"),]

  # double-check ordering
  if (!(all(row.names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp) == res$Precursor.Id))){
    setorder(res, "Precursor.Id") # try to fix ordering
    if(!(all(row.names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp) == res$Precursor.Id))){
      stop("Index error, check ordering of precursors")
    }
  }
  data.s.wide.quant.log2.qnorm.noNa.minObs.imp.cond1 = data.s.wide.quant.log2.qnorm.noNa.minObs.imp[, grep(condition_1, data$study_design$condition)]
  data.s.wide.quant.log2.qnorm.noNa.minObs.imp.cond2 = data.s.wide.quant.log2.qnorm.noNa.minObs.imp[, grep(condition_2, data$study_design$condition)]

  # Start progress bar
  if (verbose) {
    report_progress(sprintf("Precursor-level statistical testing (%d tests)", n_prec),
                   step = 6, total_steps = 8, verbose = verbose)
  }

  total = n_prec
  pb <- txtProgressBar(title = "Calculating precursor-level statistics", min = 0,
                       max = total, width = 300, style = 3)
  for (i in seq_along(prec)){
    if (verbose || plot_pdf) {  # Only show progress bar if verbose or making PDFs
      setTxtProgressBar(pb, i)
    }
    if (row.names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp.cond1)[i] != prec[i]){
      stop("precursor mismatch")
    }
    v1 = data.s.wide.quant.log2.qnorm.noNa.minObs.imp.cond1[i,]
    v2 = data.s.wide.quant.log2.qnorm.noNa.minObs.imp.cond2[i,]

    if(all(v1 == v2)){
      pvals[i] = 1
      log2fcs[i] = 0
    } else{
      t_test_result =t.test(v1, v2)
      pvals[i] = t_test_result$p.value
      log2fcs[i] = mean(v1) - mean(v2)
    }
  }
  close(pb)
  res$p_value = pvals
  res$log2_fold_change = log2fcs

  # Add precursor-level multiple testing correction
  res[, p_value_BHadj:=p.adjust(p_value, method = "BH")]

  # Visualize intermediate, precursor-level results in volcano plot
  if (!is.null(target_protein)) {
    res[, target_prot:=grepl(target_protein, Protein.Group), Precursor.Id]
  } else {
    res[, target_prot:=FALSE]
  }

  res_v = ggplot(res, aes(x = log2_fold_change, y = -log10(p_value), col = target_prot)) + geom_point() +
    scale_color_manual(values = c("darkgrey", "red")) +
    theme_bw() + ggtitle(paste0(condition_1, "/", condition_2, ", precursor-level statistics")) +
    geom_hline(yintercept = -log10(0.01), lty = 2) + geom_hline(yintercept = 0) +
    geom_vline(xintercept = c(-log2(2), log2(2)), lty = 2) + geom_vline(xintercept = 0)
  res_v
  if(plot_pdf){
    ggsave("Volcano_plot_precursorlevel.pdf", plot = res_v, height = 6, width = 6)
  }

  # Summarize to protein level
  if (verbose) {
    report_progress("Protein-level rollup and multiple testing correction", step = 7,
                   total_steps = 8, verbose = verbose)
  }

  res[, log2_fold_change_protein:=mean(log2_fold_change), Protein.Group]
  res[, n_precursors:=length(unique(Precursor.Id)), Protein.Group]
  res[, p_value_protein:=mean(p_value), Protein.Group]

  # do multiple testing correction
  pcorr = unique(res[, .(Protein.Group, p_value_protein)])
  pcorr[, p_value_BHadj_protein:=p.adjust(p_value_protein, method = "BH")]
  # Add corrected pvalues
  res = merge(res, pcorr[, .(Protein.Group, p_value_BHadj_protein)], by = "Protein.Group")

  if (verbose) {
    n_sig_proteins <- length(unique(res[p_value_BHadj_protein <= 0.05, Protein.Group]))
    report_progress(sprintf("  Found %d significant proteins (p < 0.05)",
                           n_sig_proteins),
                   type = "info", verbose = verbose)
  }

  # Add protein annotation from protein_group_annotation
  res = merge(res, protein_group_annotation, by = "Protein.Group")

  res_v_p = ggplot(res, aes(x = log2_fold_change_protein, y = -log10(p_value_BHadj_protein)))+
    geom_point() +
    scale_color_manual(values = c("darkgrey", "red")) +
    theme_bw() + ggtitle(paste0(condition_1, "/", condition_2, ", protein-level statistics")) +
    geom_hline(yintercept = -log10(0.01), lty = 2) +
    geom_vline(xintercept = c(-log2(2), log2(2)), lty = 2) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0)
  res_v_p
  if(plot_pdf){
    ggsave("Volcano_plot_proteinlevel.pdf", plot = res_v_p, height = 6, width = 6)
  }

  # Plot labels for up- and down-regulated proteins
  res_v_p + geom_text_repel(data = res[p_value_BHadj_protein<=0.01 & log2_fold_change >= 1],
                            aes(label = Protein.Names), col = "red") +
    ggtitle(paste("Proteins higher abundant in condition_1: ", condition_1))
  ggsave(paste0("Volcano_plot_candidates_up_in_",condition_1,".pdf"))

  res_v_p + geom_text_repel(data = res[p_value_BHadj_protein<=0.01 & log2_fold_change <= -1],
                            aes(label = Protein.Names), col = "blue") +
    ggtitle(paste("Proteins higher abundant in condition_2: ", condition_2))
  ggsave(paste0("Volcano_plot_candidates_up_in_",condition_2,".pdf"))

  # Write out final result table
  write.table(res, "DiffTest_result.tsv", sep = "\t", quote = F, row.names = F)

  # return results
  if (verbose) {
    report_progress("Finalizing results", step = 8, total_steps = 8, verbose = verbose)
  }

  names(data.s.long) = c("Precursor.Id", "filename", "Precursor.Quantity.log2", "normalization")
  data.s.long = merge(data.s.long, unique(data.s.wide[, .(Protein.Group, Precursor.Id)]), by = "Precursor.Id", all.x = T)
  annotation[, colname_in_matrix:=paste(filename,condition,replicate)]

  # Calculate summary statistics
  summary_stats <- calculate_summary_stats(
    diffExpr_result_dt = res,
    input_n_precursors = input_n_precursors,
    condition_1 = condition_1,
    condition_2 = condition_2,
    data_matrix = data.s.wide.quant.log2.qnorm.noNa.minObs,
    study_design = annotation
  )

  result_list = list("data_source" = data_source,
             "comparison" = paste0(condition_1, "/", condition_2),
             "conditions" = list(condition_1 = condition_1, condition_2 = condition_2),
             "diffExpr_result_dt" = res,
             "mat_quant_log2_qnorm_imp_minObs" = data.s.wide.quant.log2.qnorm.noNa.minObs.imp,
             "mat_quant_log2_qnorm" = data.s.wide.quant.log2.qnorm,
             "mat_quant_log2" = data.s.wide.quant.log2,
             "mat_quant" = data.s.wide.quant,
             "study_design" = annotation,
             "summary_stats" = summary_stats,
             "candidates_condition1" = res[p_value_BHadj_protein<=0.01 & log2_fold_change >= 1],
             "candidates_condition2" = res[p_value_BHadj_protein<=0.01 & log2_fold_change <= -1])

  # Add S3 class - now using ddata instead of diffExpr
  class(result_list) <- c("ddata", "list")

  if (verbose) {
    report_progress("Analysis complete!", type = "success", verbose = verbose)
    cat("\n")
  }

  return(result_list)
}

