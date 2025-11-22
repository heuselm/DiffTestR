# Unit tests for testDifferentialAbundance function

test_that("testDifferentialAbundance returns correct structure", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Run function
  set.seed(123)
  result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = TRUE,
    normalization_function = limma::normalizeQuantiles,
    imp_percentile = 0.001,
    imp_sd = 0.2,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  # Test structure
  expect_type(result, "list")
  expect_named(result, c("diffExpr_result_dt", "mat_quant_log2_qnorm_imp_minObs",
                         "mat_quant_log2_qnorm_imp", "mat_quant_log2_qnorm",
                         "mat_quant_log2", "mat_quant", "study_design",
                         "input_dt"))

  # Test data.table output
  expect_s3_class(result$diffExpr_result_dt, "data.table")

  # Test key columns exist
  expected_cols <- c("Precursor.Id", "Protein.Group", "log2_fold_change",
                     "p_value", "p_value_BHadj", "log2_fold_change_protein",
                     "p_value_BHadj_protein")
  expect_true(all(expected_cols %in% names(result$diffExpr_result_dt)))
})

test_that("testDifferentialAbundance produces consistent results", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Run function twice with same seed
  set.seed(123)
  result1 <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = TRUE,
    normalization_function = limma::normalizeQuantiles,
    imp_percentile = 0.001,
    imp_sd = 0.2,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  set.seed(123)
  result2 <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = TRUE,
    normalization_function = limma::normalizeQuantiles,
    imp_percentile = 0.001,
    imp_sd = 0.2,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  # Results should be identical with same seed
  expect_equal(result1$diffExpr_result_dt$log2_fold_change,
               result2$diffExpr_result_dt$log2_fold_change)
  expect_equal(result1$diffExpr_result_dt$p_value,
               result2$diffExpr_result_dt$p_value)
})

test_that("testDifferentialAbundance regression test against reference", {
  skip_if_not(file.exists("reference_results.rds"),
              "Reference results not generated yet")

  # Load reference results
  ref <- readRDS("reference_results.rds")

  # Load example data and run analysis
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  set.seed(123)
  result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = TRUE,
    normalization_function = limma::normalizeQuantiles,
    imp_percentile = 0.001,
    imp_sd = 0.2,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  # Test key statistics match reference
  expect_equal(nrow(result$diffExpr_result_dt), ref$n_precursors)
  expect_equal(length(unique(result$diffExpr_result_dt$Protein.Group)),
               ref$n_proteins)

  n_sig_precursors <- sum(result$diffExpr_result_dt$p_value_BHadj <= 0.05, na.rm = TRUE)
  expect_equal(n_sig_precursors, ref$n_significant_precursors)

  # Median fold changes should be very close
  expect_equal(median(result$diffExpr_result_dt$log2_fold_change, na.rm = TRUE),
               ref$median_log2fc_precursor, tolerance = 0.01)
  expect_equal(median(result$diffExpr_result_dt$log2_fold_change_protein, na.rm = TRUE),
               ref$median_log2fc_protein, tolerance = 0.01)
})

test_that("testDifferentialAbundance handles missing values correctly", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Count NAs in input
  input_nas <- sum(is.na(DiffTestR_example_data_wide[, -c(1:3)]))

  # Run function
  set.seed(123)
  result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = TRUE,
    normalization_function = limma::normalizeQuantiles,
    imp_percentile = 0.001,
    imp_sd = 0.2,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  # After imputation, there should be no NAs in the min obs matrix
  imputed_nas <- sum(is.na(result$mat_quant_log2_qnorm_imp_minObs))
  expect_equal(imputed_nas, 0)

  # Input should have some NAs
  expect_true(input_nas > 0)
})

test_that("testDifferentialAbundance min_n_obs parameter works", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Run with different min_n_obs values
  set.seed(123)
  result_strict <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 6,  # More strict
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  set.seed(123)
  result_lenient <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 3,  # Less strict
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  # Lenient should have more precursors passing filter
  expect_true(nrow(result_lenient$diffExpr_result_dt) >=
              nrow(result_strict$diffExpr_result_dt))
})

test_that("testDifferentialAbundance with/without normalization differs", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  set.seed(123)
  result_norm <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = TRUE,
    normalization_function = limma::normalizeQuantiles,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  set.seed(123)
  result_no_norm <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  # Results should differ
  expect_false(isTRUE(all.equal(
    result_norm$diffExpr_result_dt$log2_fold_change,
    result_no_norm$diffExpr_result_dt$log2_fold_change
  )))
})

test_that("testDifferentialAbundance validates input correctly", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Test with mismatched filenames
  bad_study_design <- copy(DiffTestR_example_study_design)
  bad_study_design$filename <- paste0("WRONG_", bad_study_design$filename)

  expect_error(
    testDifferentialAbundance(
      input_dt = DiffTestR_example_data_wide,
      study_design = bad_study_design,
      condition_1 = "A",
      condition_2 = "B",
      plot_pdf = FALSE,
      write_tsv_tables = FALSE
    ),
    "filename in study design must match column name in input_dt"
  )
})
