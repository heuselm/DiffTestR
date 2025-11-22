# Unit tests for plotDifferentialAbundanceOverview function

test_that("plotDifferentialAbundanceOverview returns correct structure", {
  # Load example data and run differential abundance test
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  set.seed(123)
  diff_result <- testDifferentialAbundance(
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

  # Add comparison column as expected by plotDifferentialAbundanceOverview
  diff_result$diffExpr_result_dt[, comparison := "A_vs_B"]

  # Suppress plot output during testing
  pdf(NULL)
  result <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = diff_result$diffExpr_result_dt,
    label_prefix = tempfile("test_plot"),
    browsable_html = FALSE,
    heatmap = TRUE
  )
  dev.off()

  # Test structure
  expect_type(result, "list")
  expect_named(result, c("diffExpr_result_dt", "volcanos", "heatmap", "FC_matrix"))

  # Test components
  expect_s3_class(result$diffExpr_result_dt, "data.table")
  expect_s3_class(result$volcanos, "gg")
  expect_s3_class(result$volcanos, "ggplot")
  expect_true(is.matrix(result$FC_matrix))
})

test_that("plotDifferentialAbundanceOverview processes protein-level results correctly", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  set.seed(123)
  diff_result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  diff_result$diffExpr_result_dt[, comparison := "A_vs_B"]

  pdf(NULL)
  result <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = diff_result$diffExpr_result_dt,
    label_prefix = tempfile("test_plot"),
    browsable_html = FALSE,
    heatmap = FALSE
  )
  dev.off()

  # Should have unique protein-level results
  protein_dt <- result$diffExpr_result_dt
  expect_true("Protein.Group" %in% names(protein_dt))
  expect_true("log2_fold_change_protein" %in% names(protein_dt))
  expect_true("p_value_BHadj_protein" %in% names(protein_dt))

  # Should have one row per protein per comparison
  n_unique <- nrow(unique(protein_dt[, .(Protein.Group, comparison)]))
  expect_equal(nrow(protein_dt), n_unique)
})

test_that("plotDifferentialAbundanceOverview significance thresholds work", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  set.seed(123)
  diff_result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  diff_result$diffExpr_result_dt[, comparison := "A_vs_B"]

  # Test with different significance thresholds
  pdf(NULL)
  result_strict <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = diff_result$diffExpr_result_dt,
    significance_threshold_p = 0.01,  # More strict
    significance_threshold_fc = 3,     # 3-fold change
    label_prefix = tempfile("test_strict"),
    browsable_html = FALSE,
    heatmap = TRUE
  )

  result_lenient <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = diff_result$diffExpr_result_dt,
    significance_threshold_p = 0.1,   # More lenient
    significance_threshold_fc = 1.5,  # 1.5-fold change
    label_prefix = tempfile("test_lenient"),
    browsable_html = FALSE,
    heatmap = TRUE
  )
  dev.off()

  # Lenient should have more or equal significant proteins in FC matrix
  expect_true(nrow(result_lenient$FC_matrix) >= nrow(result_strict$FC_matrix))
})

test_that("plotDifferentialAbundanceOverview handles multiple comparisons", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Create multiple comparisons by duplicating and modifying
  set.seed(123)
  diff_result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  # Create two comparisons
  comp1 <- copy(diff_result$diffExpr_result_dt)
  comp1[, comparison := "A_vs_B"]

  comp2 <- copy(diff_result$diffExpr_result_dt)
  comp2[, comparison := "B_vs_A"]
  comp2[, log2_fold_change_protein := -log2_fold_change_protein]

  multi_comp <- rbind(comp1, comp2)

  pdf(NULL)
  result <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = multi_comp,
    label_prefix = tempfile("test_multi"),
    browsable_html = FALSE,
    heatmap = TRUE
  )
  dev.off()

  # Should handle multiple comparisons
  expect_true(length(unique(result$diffExpr_result_dt$comparison)) == 2)

  # FC matrix should have 2 columns (one per comparison)
  expect_equal(ncol(result$FC_matrix), 2)
})

test_that("plotDifferentialAbundanceOverview protein name cleaning works", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  set.seed(123)
  diff_result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  diff_result$diffExpr_result_dt[, comparison := "A_vs_B"]

  pdf(NULL)
  result <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = diff_result$diffExpr_result_dt,
    remove_from_protein_names = "_HUMAN|_YEAST|_ECOLI",
    label_prefix = tempfile("test_clean"),
    browsable_html = FALSE,
    heatmap = FALSE
  )
  dev.off()

  # Should have cleaned protein names
  expect_true("Protein.Names.Short" %in% names(result$diffExpr_result_dt))

  # Cleaned names should not contain the removed patterns
  has_human <- any(grepl("_HUMAN", result$diffExpr_result_dt$Protein.Names.Short))
  has_yeast <- any(grepl("_YEAST", result$diffExpr_result_dt$Protein.Names.Short))
  has_ecoli <- any(grepl("_ECOLI", result$diffExpr_result_dt$Protein.Names.Short))

  expect_false(has_human)
  expect_false(has_yeast)
  expect_false(has_ecoli)
})

test_that("plotDifferentialAbundanceOverview label_prefix works", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  set.seed(123)
  diff_result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  diff_result$diffExpr_result_dt[, comparison := "A_vs_B"]

  # Test with custom label prefix
  temp_label <- tempfile("custom_prefix")

  pdf(NULL)
  result <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = diff_result$diffExpr_result_dt,
    label_prefix = temp_label,
    browsable_html = FALSE,
    heatmap = FALSE
  )
  dev.off()

  # Volcano plot should have custom title
  expect_true(grepl(basename(temp_label), as.character(result$volcanos$labels$title)))
})

test_that("plotDifferentialAbundanceOverview heatmap parameter works", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  set.seed(123)
  diff_result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  diff_result$diffExpr_result_dt[, comparison := "A_vs_B"]

  # Test with heatmap = TRUE
  pdf(NULL)
  result_with_hm <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = diff_result$diffExpr_result_dt,
    label_prefix = tempfile("with_hm"),
    browsable_html = FALSE,
    heatmap = TRUE
  )
  dev.off()

  # Should have heatmap
  expect_false(is.null(result_with_hm$heatmap))
  expect_true(is.matrix(result_with_hm$FC_matrix))

  # Test with heatmap = FALSE
  pdf(NULL)
  result_no_hm <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = diff_result$diffExpr_result_dt,
    label_prefix = tempfile("no_hm"),
    browsable_html = FALSE,
    heatmap = FALSE
  )
  dev.off()

  # Both should still return structure, but heatmap might be NULL or empty
  expect_type(result_no_hm, "list")
})

test_that("plotDifferentialAbundanceOverview protein highlighting works", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  set.seed(123)
  diff_result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  diff_result$diffExpr_result_dt[, comparison := "A_vs_B"]

  # Test with custom highlight tag
  pdf(NULL)
  result <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = diff_result$diffExpr_result_dt,
    protein_highlight_tag = "HUMAN",
    label_prefix = tempfile("highlight"),
    browsable_html = FALSE,
    heatmap = FALSE
  )
  dev.off()

  # Should have target_prot column
  expect_true("target_prot" %in% names(result$diffExpr_result_dt))

  # target_prot should be logical
  expect_type(result$diffExpr_result_dt$target_prot, "logical")

  # Some proteins should be highlighted if pattern exists in data
  if (any(grepl("HUMAN", result$diffExpr_result_dt$Protein.Group))) {
    expect_true(any(result$diffExpr_result_dt$target_prot))
  }
})

test_that("plotDifferentialAbundanceOverview FC matrix has correct dimensions", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  set.seed(123)
  diff_result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  # Create two comparisons
  comp1 <- copy(diff_result$diffExpr_result_dt)
  comp1[, comparison := "comp1"]

  comp2 <- copy(diff_result$diffExpr_result_dt)
  comp2[, comparison := "comp2"]

  multi_comp <- rbind(comp1, comp2)

  pdf(NULL)
  result <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = multi_comp,
    significance_threshold_p = 0.05,
    significance_threshold_fc = 2,
    label_prefix = tempfile("fc_matrix"),
    browsable_html = FALSE,
    heatmap = TRUE
  )
  dev.off()

  # FC matrix should have correct dimensions
  expect_true(is.matrix(result$FC_matrix))

  # Number of columns should match number of comparisons
  expect_equal(ncol(result$FC_matrix), 2)

  # Rownames should be protein names
  expect_true(length(rownames(result$FC_matrix)) > 0)

  # All values should be numeric
  expect_type(result$FC_matrix, "double")
})

test_that("plotDifferentialAbundanceOverview handles n_precursors column", {
  # Load example data
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  set.seed(123)
  diff_result <- testDifferentialAbundance(
    input_dt = DiffTestR_example_data_wide,
    study_design = DiffTestR_example_study_design,
    condition_1 = "A",
    condition_2 = "B",
    min_n_obs = 4,
    normalize_data = FALSE,
    plot_pdf = FALSE,
    write_tsv_tables = FALSE
  )

  diff_result$diffExpr_result_dt[, comparison := "A_vs_B"]

  pdf(NULL)
  result <- plotDifferentialAbundanceOverview(
    diffExpr_result_dt = diff_result$diffExpr_result_dt,
    label_prefix = tempfile("n_precursors"),
    browsable_html = FALSE,
    heatmap = FALSE
  )
  dev.off()

  # Should have n_proteins column calculated
  expect_true("n_proteins" %in% names(result$diffExpr_result_dt))

  # n_proteins should be numeric and positive
  expect_type(result$diffExpr_result_dt$n_proteins, "integer")
  expect_true(all(result$diffExpr_result_dt$n_proteins > 0))
})
