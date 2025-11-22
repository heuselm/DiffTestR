# Unit tests for testMultipleComparisons function

test_that("testMultipleComparisons returns correct structure", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Create temporary output directory
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_multi_comp_1")

  # Run with explicit comparisons
  set.seed(123)
  comparisons <- list(
    list(condition_1 = "A", condition_2 = "B")
  )

  result <- testMultipleComparisons(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    comparisons = comparisons,
    min_n_obs = 4,
    normalize_data = TRUE,
    normalization_function = limma::normalizeQuantiles,
    imp_percentile = 0.001,
    imp_sd = 0.2,
    output_dir = output_dir,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE,
    verbose = FALSE
  )

  # Test structure
  expect_s3_class(result, "multiDiffExpr")
  expect_type(result, "list")

  # Check that key components exist
  expected_names <- c("individual_results", "combined_results_dt", "comparisons",
                      "overview_stats", "n_comparisons", "n_successful",
                      "error_log", "output_dir", "study_design")
  expect_true(all(expected_names %in% names(result)))

  # Test individual results
  expect_type(result$individual_results, "list")
  expect_equal(length(result$individual_results), 1)
  expect_true("A_vs_B" %in% names(result$individual_results))

  # Each individual result should be a diffExpr object
  expect_s3_class(result$individual_results$A_vs_B, "diffExpr")

  # Test combined results
  expect_s3_class(result$combined_results_dt, "data.table")
  expect_true("comparison" %in% names(result$combined_results_dt))

  # Test counts
  expect_equal(result$n_comparisons, 1)
  expect_equal(result$n_successful, 1)

  # Test error log (should be NULL for successful run)
  expect_null(result$error_log)

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("testMultipleComparisons auto-generates pairwise comparisons", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Create a study design with 3 conditions
  design_3cond <- copy(DiffTestR_example_study_design)
  # Add a third condition by duplicating some samples
  design_extra <- copy(design_3cond[condition == "A"])
  design_extra[, condition := "C"]
  design_3cond <- rbind(design_3cond, design_extra)

  # Create temporary output directory
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_multi_comp_auto")

  # Run with auto-generated comparisons (comparisons = NULL)
  set.seed(123)
  result <- testMultipleComparisons(
    input_dt = DiffTestR_example_data_wide,
    study_design = design_3cond,
    comparisons = NULL,  # Auto-generate
    min_n_obs = 4,
    normalize_data = FALSE,
    output_dir = output_dir,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE,
    verbose = FALSE
  )

  # With 3 conditions (A, B, C), should generate 3 pairwise comparisons
  expect_equal(result$n_comparisons, 3)

  # Check that all pairwise combinations are present
  comparison_names <- names(result$individual_results)
  expect_true("A_vs_B" %in% comparison_names || "B_vs_A" %in% comparison_names)
  expect_true("A_vs_C" %in% comparison_names || "C_vs_A" %in% comparison_names)
  expect_true("B_vs_C" %in% comparison_names || "C_vs_B" %in% comparison_names)

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("testMultipleComparisons handles errors gracefully with stop_on_error=FALSE", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Create temporary output directory
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_multi_comp_error")

  # Create comparisons where one will fail (non-existent condition)
  comparisons <- list(
    list(condition_1 = "A", condition_2 = "B"),  # Valid
    list(condition_1 = "A", condition_2 = "INVALID")  # Invalid - should fail
  )

  # This should throw an error during validation
  expect_error(
    testMultipleComparisons(
      input_dt = DiffTestR_example_data_wide,
      study_design = DiffTestR_example_study_design,
      comparisons = comparisons,
      output_dir = output_dir,
      plot_pdf = FALSE,
      write_tsv_tables = FALSE,
      stop_on_error = FALSE,
      verbose = FALSE
    ),
    "condition_2 'INVALID' not found in study design"
  )

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("testMultipleComparisons creates subdirectories correctly", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Create temporary output directory
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_multi_comp_dirs")

  # Run with explicit comparisons
  set.seed(123)
  comparisons <- list(
    list(condition_1 = "A", condition_2 = "B")
  )

  result <- testMultipleComparisons(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    comparisons = comparisons,
    min_n_obs = 4,
    output_dir = output_dir,
    plot_pdf = FALSE,
    write_tsv_tables = TRUE,  # Should create files
    verbose = FALSE
  )

  # Check that output directory was created
  expect_true(dir.exists(output_dir))

  # Check that subdirectory was created
  comparison_dir <- file.path(output_dir, "A_vs_B")
  expect_true(dir.exists(comparison_dir))

  # Check that result file was created (because write_tsv_tables = TRUE)
  result_file <- file.path(comparison_dir, "DiffTest_result.tsv")
  expect_true(file.exists(result_file))

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("testMultipleComparisons combined results contain all comparisons", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Create a study design with 3 conditions
  design_3cond <- copy(DiffTestR_example_study_design)
  design_extra <- copy(design_3cond[condition == "A"])
  design_extra[, condition := "C"]
  design_3cond <- rbind(design_3cond, design_extra)

  # Create temporary output directory
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_multi_comp_combined")

  # Run with explicit comparisons
  set.seed(123)
  comparisons <- list(
    list(condition_1 = "A", condition_2 = "B"),
    list(condition_1 = "A", condition_2 = "C")
  )

  result <- testMultipleComparisons(
    input_dt = DiffTestR_example_data_wide,
    study_design = design_3cond,
    comparisons = comparisons,
    min_n_obs = 4,
    normalize_data = FALSE,
    output_dir = output_dir,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE,
    verbose = FALSE
  )

  # Check combined results
  expect_s3_class(result$combined_results_dt, "data.table")

  # Should have comparison column
  expect_true("comparison" %in% names(result$combined_results_dt))

  # Should have both comparisons represented
  unique_comparisons <- unique(result$combined_results_dt$comparison)
  expect_equal(length(unique_comparisons), 2)
  expect_true("A_vs_B" %in% unique_comparisons)
  expect_true("A_vs_C" %in% unique_comparisons)

  # Each comparison should have the same number of rows (same precursors tested)
  comp_counts <- result$combined_results_dt[, .N, by = comparison]
  expect_equal(length(unique(comp_counts$N)), 1)  # All counts should be equal

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("testMultipleComparisons validates input correctly", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Test with invalid study design (missing required columns)
  bad_study_design <- copy(DiffTestR_example_study_design)
  bad_study_design[, condition := NULL]

  expect_error(
    testMultipleComparisons(
      input_dt = DiffTestR_example_data_wide,
      study_design = bad_study_design,
      plot_pdf = FALSE,
      write_tsv_tables = FALSE,
      verbose = FALSE
    ),
    "study design must contain columns filename, condition and replicate"
  )

  # Test with only one condition (should fail)
  single_condition_design <- DiffTestR_example_study_design[condition == "A"]

  expect_error(
    testMultipleComparisons(
      input_dt = DiffTestR_example_data_wide,
      study_design = single_condition_design,
      comparisons = NULL,
      plot_pdf = FALSE,
      write_tsv_tables = FALSE,
      verbose = FALSE
    ),
    "At least 2 unique conditions are required"
  )

  # Test with invalid comparisons format
  expect_error(
    testMultipleComparisons(
      input_dt = DiffTestR_example_data_wide,
      study_design = DiffTestR_example_study_design,
      comparisons = "not a list",
      plot_pdf = FALSE,
      write_tsv_tables = FALSE,
      verbose = FALSE
    ),
    "comparisons must be a list"
  )

  # Test with incomplete comparison specification
  bad_comparisons <- list(
    list(condition_1 = "A")  # Missing condition_2
  )

  expect_error(
    testMultipleComparisons(
      input_dt = DiffTestR_example_data_wide,
      study_design = DiffTestR_example_study_design,
      comparisons = bad_comparisons,
      plot_pdf = FALSE,
      write_tsv_tables = FALSE,
      verbose = FALSE
    ),
    "must contain condition_1 and condition_2"
  )
})

test_that("testMultipleComparisons print and summary methods work", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Create temporary output directory
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_multi_comp_print")

  # Run analysis
  set.seed(123)
  comparisons <- list(
    list(condition_1 = "A", condition_2 = "B")
  )

  result <- testMultipleComparisons(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    comparisons = comparisons,
    min_n_obs = 4,
    output_dir = output_dir,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE,
    verbose = FALSE
  )

  # Test print method (should not error)
  expect_output(print(result), "Multiple Differential Abundance Comparisons")
  expect_output(print(result), "Total comparisons: 1")
  expect_output(print(result), "Successful: 1")

  # Test summary method (should not error)
  expect_output(summary(result), "Multiple Differential Abundance Comparisons")
  expect_output(summary(result), "Detailed Summary by Comparison")
  expect_output(summary(result), "A_vs_B")

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("testMultipleComparisons handles single comparison", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Create temporary output directory
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_multi_comp_single")

  # Run with single comparison
  set.seed(123)
  comparisons <- list(
    list(condition_1 = "A", condition_2 = "B")
  )

  result <- testMultipleComparisons(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    comparisons = comparisons,
    min_n_obs = 4,
    output_dir = output_dir,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE,
    verbose = FALSE
  )

  # Should work just fine with a single comparison
  expect_equal(result$n_comparisons, 1)
  expect_equal(result$n_successful, 1)
  expect_equal(length(result$individual_results), 1)

  # Combined results should have same structure
  expect_s3_class(result$combined_results_dt, "data.table")
  expect_equal(length(unique(result$combined_results_dt$comparison)), 1)

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("testMultipleComparisons respects verbose parameter", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Create temporary output directory
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_multi_comp_verbose")

  # Run with verbose = FALSE (should produce minimal output)
  set.seed(123)
  comparisons <- list(
    list(condition_1 = "A", condition_2 = "B")
  )

  # Capture output
  output_silent <- capture.output(
    result_silent <- testMultipleComparisons(
      input_dt = DiffTestR_example_data_wide,
      study_design = DiffTestR_example_study_design,
      comparisons = comparisons,
      min_n_obs = 4,
      output_dir = output_dir,
      plot_pdf = FALSE,
      write_tsv_tables = FALSE,
      verbose = FALSE
    ),
    type = "output"
  )

  # Should produce minimal or no output
  expect_true(length(output_silent) < 10)

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})

test_that("testMultipleComparisons preserves working directory", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Store original working directory
  original_wd <- getwd()

  # Create temporary output directory
  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_multi_comp_wd")

  # Run analysis
  set.seed(123)
  comparisons <- list(
    list(condition_1 = "A", condition_2 = "B")
  )

  result <- testMultipleComparisons(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    comparisons = comparisons,
    min_n_obs = 4,
    output_dir = output_dir,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE,
    verbose = FALSE
  )

  # Working directory should be unchanged
  expect_equal(getwd(), original_wd)

  # Cleanup
  unlink(output_dir, recursive = TRUE)
})
