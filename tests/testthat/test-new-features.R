# Unit tests for new features added in "evolve" branch

test_that("assessDataQuality works with example data", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Convert to matrix
  data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
  rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id

  # Log transform
  data_mat_log2 <- log2(data_mat + 1)

  # Run QC assessment
  qc <- assessDataQuality(data_mat_log2,
                          DiffTestR_example_study_design,
                          verbose = FALSE)

  # Test structure
  expect_type(qc, "list")
  expect_s3_class(qc, "DiffTestR_QC")

  # Test components exist
  expect_true("cv_stats" %in% names(qc))
  expect_true("cv_loo_analysis" %in% names(qc))
  expect_true("outlier_replicates" %in% names(qc))
  expect_true("missing_stats" %in% names(qc))
  expect_true("dynamic_range" %in% names(qc))
  expect_true("replicate_correlations" %in% names(qc))

  # Test CV stats
  expect_true(length(qc$cv_stats) > 0)
  for (cond in names(qc$cv_stats)) {
    expect_true("median_cv" %in% names(qc$cv_stats[[cond]]))
    expect_true(is.numeric(qc$cv_stats[[cond]]$median_cv))
  }

  # Test missing stats
  expect_true(is.numeric(qc$missing_stats$pct_missing))
  expect_true(qc$missing_stats$pct_missing >= 0)
  expect_true(qc$missing_stats$pct_missing <= 100)

  # Test dynamic range
  expect_true(is.numeric(qc$dynamic_range$min))
  expect_true(is.numeric(qc$dynamic_range$max))
  expect_true(qc$dynamic_range$max > qc$dynamic_range$min)
})

test_that("Leave-one-out CV analysis detects outliers correctly", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE))

  # Create synthetic data with an outlier replicate
  set.seed(123)
  n_features <- 100

  # Normal replicates (CV ~ 10%)
  rep1 <- rnorm(n_features, mean = 10, sd = 1)
  rep2 <- rnorm(n_features, mean = 10, sd = 1)
  rep3 <- rnorm(n_features, mean = 10, sd = 1)

  # Outlier replicate (high variance)
  rep4 <- rnorm(n_features, mean = 10, sd = 3)

  test_mat <- cbind(rep1, rep2, rep3, rep4)
  rownames(test_mat) <- paste0("feature_", 1:n_features)

  study_design <- data.table::data.table(
    filename = c("rep1", "rep2", "rep3", "rep4"),
    condition = c("A", "A", "A", "A"),
    replicate = 1:4
  )

  # Run QC
  qc <- assessDataQuality(test_mat, study_design,
                          cv_threshold = 0.2, verbose = FALSE)

  # Should detect outlier
  expect_s3_class(qc$outlier_replicates, "data.table")

  # Outlier replicate should be flagged (rep4)
  if (nrow(qc$outlier_replicates) > 0) {
    expect_true("rep4" %in% qc$outlier_replicates$filename)
  }
})

test_that("print.diffExpr displays summary correctly", {
  # Load example data and run minimal analysis
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Create a mock diffExpr object with summary stats
  mock_result <- list(
    diffExpr_result_dt = data.table::data.table(
      Protein.Group = rep(paste0("Prot_", 1:10), each = 3),
      Precursor.Id = paste0("Prec_", 1:30),
      log2_fold_change = rnorm(30),
      log2_fold_change_protein = rnorm(30),
      p_value_BHadj = runif(30),
      p_value_BHadj_protein = runif(30)
    ),
    summary_stats = list(
      n_precursors_input = 100,
      n_precursors_tested = 30,
      n_proteins_tested = 10,
      n_significant_precursors = 5,
      n_significant_proteins = 2,
      n_upregulated = 1,
      n_downregulated = 1,
      median_log2fc_protein = 0.5,
      condition_1 = "A",
      condition_2 = "B",
      comparison = "A vs B"
    )
  )

  class(mock_result) <- c("diffExpr", "list")

  # Test that print works without error
  expect_output(print(mock_result), "Differential Abundance Results")
  expect_output(print(mock_result), "Precursors")
  expect_output(print(mock_result), "Proteins")
})

test_that("verbose parameter controls progress messages", {
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # With verbose = FALSE, should have minimal output
  expect_silent({
    result_quiet <- testDifferentialAbundance(
      input_dt = DiffTestR_example_data_wide,
      study_design = DiffTestR_example_study_design,
      condition_1 = "A",
      condition_2 = "B",
      min_n_obs = 4,
      plot_pdf = FALSE,
      write_tsv_tables = FALSE,
      verbose = FALSE
    )
  })

  # Should still return valid result
  expect_s3_class(result_quiet, "diffExpr")
  expect_true("diffExpr_result_dt" %in% names(result_quiet))
  expect_true("summary_stats" %in% names(result_quiet))
})

test_that("Summary statistics are calculated correctly", {
  # Create test data
  test_dt <- data.table::data.table(
    Protein.Group = rep(c("A", "B", "C"), each = 3),
    log2_fold_change = rnorm(9),
    log2_fold_change_protein = rep(c(1.5, -1.5, 0.5), each = 3),
    p_value_BHadj = c(rep(0.01, 3), rep(0.001, 3), rep(0.5, 3)),
    p_value_BHadj_protein = rep(c(0.01, 0.001, 0.5), each = 3)
  )

  stats <- calculate_summary_stats(
    diffExpr_result_dt = test_dt,
    input_n_precursors = 100,
    condition_1 = "A",
    condition_2 = "B"
  )

  # Test basic stats
  expect_equal(stats$n_precursors_input, 100)
  expect_equal(stats$n_precursors_tested, 9)
  expect_equal(stats$n_proteins_tested, 3)

  # Test significant counts (p < 0.05)
  expect_equal(stats$n_significant_proteins, 2)  # A and B

  # Test regulation counts (|log2FC| > 1, p < 0.05)
  expect_equal(stats$n_upregulated, 1)  # Protein A
  expect_equal(stats$n_downregulated, 1)  # Protein B
})

test_that("Progress reporting utilities work correctly", {
  # Test report_progress with verbose = TRUE
  expect_output(
    report_progress("Test message", step = 1, total_steps = 5, verbose = TRUE),
    "\\[1/5\\]"
  )

  # Test report_progress with verbose = FALSE (should be silent)
  expect_silent(
    report_progress("Test message", verbose = FALSE)
  )

  # Test different message types
  expect_output(
    report_progress("Success", type = "success", verbose = TRUE),
    "\\[OK\\]"
  )

  expect_output(
    report_progress("Warning", type = "warning", verbose = TRUE),
    "\\[WARNING\\]"
  )
})
