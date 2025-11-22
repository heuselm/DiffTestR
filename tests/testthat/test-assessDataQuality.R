# Unit tests for assessDataQuality function

test_that("assessDataQuality returns correct structure", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Convert to matrix format
  data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
  rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id

  # Run function
  result <- assessDataQuality(
    data_matrix = data_mat,
    study_design = DiffTestR_example_study_design,
    verbose = FALSE
  )

  # Test structure
  expect_type(result, "list")
  expect_s3_class(result, "DiffTestR_QC")

  # Test components exist
  expected_components <- c("cv_stats", "cv_loo_analysis", "outlier_replicates",
                          "missing_stats", "dynamic_range", "replicate_correlations")
  expect_true(all(expected_components %in% names(result)))
})

test_that("assessDataQuality CV statistics are valid", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
  rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id

  result <- assessDataQuality(data_mat, DiffTestR_example_study_design, verbose = FALSE)

  # Check CV stats for each condition
  conditions <- unique(DiffTestR_example_study_design$condition)
  expect_equal(length(result$cv_stats), length(conditions))

  for (cond in conditions) {
    cv_stat <- result$cv_stats[[cond]]

    # CV values should be numeric and positive
    expect_type(cv_stat$median_cv, "double")
    expect_true(cv_stat$median_cv > 0)
    expect_true(cv_stat$mean_cv > 0)

    # Standard deviation should be non-negative
    expect_true(cv_stat$sd_cv >= 0)

    # Number of samples should match study design
    n_samples <- sum(DiffTestR_example_study_design$condition == cond)
    expect_equal(cv_stat$n_samples, n_samples)
  }
})

test_that("assessDataQuality leave-one-out analysis works", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
  rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id

  result <- assessDataQuality(data_mat, DiffTestR_example_study_design, verbose = FALSE)

  # Should have LOO results for conditions with >= 3 replicates
  conditions <- unique(DiffTestR_example_study_design$condition)

  for (cond in conditions) {
    n_reps <- sum(DiffTestR_example_study_design$condition == cond)

    if (n_reps >= 3) {
      expect_true(cond %in% names(result$cv_loo_analysis))

      loo_dt <- result$cv_loo_analysis[[cond]]

      # Should have one row per replicate
      expect_equal(nrow(loo_dt), n_reps)

      # Required columns
      expect_true(all(c("condition", "replicate", "filename", "baseline_cv",
                       "loo_cv", "cv_improvement", "is_outlier") %in% names(loo_dt)))

      # CV improvement should be numeric
      expect_type(loo_dt$cv_improvement, "double")

      # is_outlier should be logical
      expect_type(loo_dt$is_outlier, "logical")
    }
  }
})

test_that("assessDataQuality outlier detection works with threshold", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
  rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id

  # Test with very strict threshold (should flag more outliers)
  result_strict <- assessDataQuality(
    data_mat,
    DiffTestR_example_study_design,
    cv_threshold = 0.01,  # 1% improvement flags outlier
    verbose = FALSE
  )

  # Test with lenient threshold (should flag fewer outliers)
  result_lenient <- assessDataQuality(
    data_mat,
    DiffTestR_example_study_design,
    cv_threshold = 0.9,  # 90% improvement needed
    verbose = FALSE
  )

  # Outlier replicates should be data.table
  expect_s3_class(result_strict$outlier_replicates, "data.table")
  expect_s3_class(result_lenient$outlier_replicates, "data.table")

  # Strict threshold should flag same or more outliers than lenient
  expect_true(nrow(result_strict$outlier_replicates) >=
              nrow(result_lenient$outlier_replicates))
})

test_that("assessDataQuality missing value analysis is correct", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
  rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id

  result <- assessDataQuality(data_mat, DiffTestR_example_study_design, verbose = FALSE)

  # Calculate expected values
  expected_n_missing <- sum(is.na(data_mat))
  expected_total <- length(data_mat)
  expected_pct <- 100 * expected_n_missing / expected_total

  # Check missing stats
  expect_equal(result$missing_stats$n_missing, expected_n_missing)
  expect_equal(result$missing_stats$total_values, expected_total)
  expect_equal(result$missing_stats$pct_missing, expected_pct)

  # Missing per sample should have entry for each sample
  expect_equal(length(result$missing_stats$missing_per_sample), ncol(data_mat))

  # Missing per precursor should have entry for each precursor
  expect_equal(length(result$missing_stats$missing_per_precursor), nrow(data_mat))

  # All percentages should be between 0 and 100
  expect_true(all(result$missing_stats$missing_per_sample >= 0 &
                  result$missing_stats$missing_per_sample <= 100))
  expect_true(all(result$missing_stats$missing_per_precursor >= 0 &
                  result$missing_stats$missing_per_precursor <= 100))
})

test_that("assessDataQuality dynamic range is calculated correctly", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
  rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id

  result <- assessDataQuality(data_mat, DiffTestR_example_study_design, verbose = FALSE)

  # Calculate expected values
  expected_min <- min(data_mat, na.rm = TRUE)
  expected_max <- max(data_mat, na.rm = TRUE)
  expected_range <- expected_max - expected_min

  # Check dynamic range
  expect_equal(result$dynamic_range$min, expected_min)
  expect_equal(result$dynamic_range$max, expected_max)
  expect_equal(result$dynamic_range$range, expected_range)

  # Range should be positive
  expect_true(result$dynamic_range$range > 0)
})

test_that("assessDataQuality replicate correlation analysis works", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
  rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id

  result <- assessDataQuality(data_mat, DiffTestR_example_study_design, verbose = FALSE)

  # Should have correlations for each condition
  conditions <- unique(DiffTestR_example_study_design$condition)
  expect_equal(length(result$replicate_correlations), length(conditions))

  for (cond in conditions) {
    cor_result <- result$replicate_correlations[[cond]]

    # Should have correlation matrix
    expect_true(is.matrix(cor_result$cor_matrix))

    # Correlation values should be between -1 and 1
    expect_true(cor_result$median_cor >= -1 && cor_result$median_cor <= 1)
    expect_true(cor_result$mean_cor >= -1 && cor_result$mean_cor <= 1)
    expect_true(cor_result$min_cor >= -1 && cor_result$min_cor <= 1)
    expect_true(cor_result$max_cor >= -1 && cor_result$max_cor <= 1)

    # Matrix should be square
    n_samples <- sum(DiffTestR_example_study_design$condition == cond)
    expect_equal(nrow(cor_result$cor_matrix), n_samples)
    expect_equal(ncol(cor_result$cor_matrix), n_samples)
  }
})

test_that("assessDataQuality validates input correctly", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
  rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id

  # Test with non-matrix input
  expect_error(
    assessDataQuality(
      data_matrix = DiffTestR_example_data_wide,
      study_design = DiffTestR_example_study_design,
      verbose = FALSE
    ),
    "data_matrix must be a numeric matrix"
  )

  # Test with missing columns in study design
  bad_study_design <- copy(DiffTestR_example_study_design)
  bad_study_design[, condition := NULL]

  expect_error(
    assessDataQuality(
      data_matrix = data_mat,
      study_design = bad_study_design,
      verbose = FALSE
    ),
    "study_design must have columns"
  )

  # Test with mismatched column names
  colnames(data_mat) <- paste0("WRONG_", colnames(data_mat))

  expect_error(
    assessDataQuality(
      data_matrix = data_mat,
      study_design = DiffTestR_example_study_design,
      verbose = FALSE
    ),
    "Not all column names in data_matrix match study_design"
  )
})

test_that("assessDataQuality verbose parameter works", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
  rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id

  # Test verbose = FALSE (should not print)
  expect_silent(
    assessDataQuality(data_mat, DiffTestR_example_study_design, verbose = FALSE)
  )

  # Test verbose = TRUE (should print output)
  expect_output(
    assessDataQuality(data_mat, DiffTestR_example_study_design, verbose = TRUE),
    "Data Quality Assessment"
  )
})

test_that("assessDataQuality print method works", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  data_mat <- as.matrix(DiffTestR_example_data_wide[, 4:11])
  rownames(data_mat) <- DiffTestR_example_data_wide$Precursor.Id

  result <- assessDataQuality(data_mat, DiffTestR_example_study_design, verbose = FALSE)

  # Test print output
  expect_output(print(result), "DiffTestR Data Quality Report")
  expect_output(print(result), "Coefficient of Variation")
  expect_output(print(result), "Missing Values")
  expect_output(print(result), "Dynamic Range")
  expect_output(print(result), "Replicate Correlation")
})

test_that("assessDataQuality handles edge case with few replicates", {
  # Create minimal study design with only 2 replicates per condition
  minimal_study_design <- data.table(
    filename = c("sample1", "sample2", "sample3", "sample4"),
    condition = c("A", "A", "B", "B"),
    replicate = c(1, 2, 1, 2)
  )

  # Create minimal data matrix
  set.seed(123)
  minimal_mat <- matrix(
    rnorm(100 * 4, mean = 20, sd = 2),
    nrow = 100,
    ncol = 4
  )
  colnames(minimal_mat) <- minimal_study_design$filename
  rownames(minimal_mat) <- paste0("precursor_", 1:100)

  # Should run without errors but skip LOO analysis
  result <- assessDataQuality(
    minimal_mat,
    minimal_study_design,
    verbose = FALSE
  )

  # LOO analysis should be empty (need >= 3 replicates)
  expect_equal(length(result$cv_loo_analysis), 0)

  # Other analyses should still work
  expect_true(!is.null(result$cv_stats))
  expect_true(!is.null(result$missing_stats))
  expect_true(!is.null(result$replicate_correlations))
})
