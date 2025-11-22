# Unit tests for example datasets

test_that("DiffTestR_example_data_long loads correctly", {
  data(DiffTestR_example_data_long, envir = environment())

  # Test structure
  expect_s3_class(DiffTestR_example_data_long, "data.table")

  # Test key columns
  expected_cols <- c("Run", "Protein.Group", "Protein.Names",
                     "Precursor.Id", "Precursor.Quantity")
  expect_true(all(expected_cols %in% names(DiffTestR_example_data_long)))

  # Test dimensions
  expect_true(nrow(DiffTestR_example_data_long) > 0)
  expect_equal(ncol(DiffTestR_example_data_long), 5)
})

test_that("DiffTestR_example_data_wide loads correctly", {
  data(DiffTestR_example_data_wide, envir = environment())

  # Test structure
  expect_s3_class(DiffTestR_example_data_wide, "data.table")

  # Test key columns
  expect_true("Protein.Group" %in% names(DiffTestR_example_data_wide))
  expect_true("Protein.Names" %in% names(DiffTestR_example_data_wide))
  expect_true("Precursor.Id" %in% names(DiffTestR_example_data_wide))

  # Should have 3 annotation columns + 8 sample columns = 11 total
  expect_equal(ncol(DiffTestR_example_data_wide), 11)

  # Test dimensions
  expect_true(nrow(DiffTestR_example_data_wide) > 0)
})

test_that("DiffTestR_example_study_design loads correctly", {
  data(DiffTestR_example_study_design, envir = environment())

  # Test structure
  expect_s3_class(DiffTestR_example_study_design, "data.table")

  # Test columns
  expect_named(DiffTestR_example_study_design,
               c("filename", "condition", "replicate"))

  # Test content
  expect_equal(nrow(DiffTestR_example_study_design), 8)
  expect_true(all(DiffTestR_example_study_design$condition %in% c("A", "B")))
  expect_true(all(DiffTestR_example_study_design$replicate %in% 1:4))

  # Should have 4 replicates per condition
  expect_equal(sum(DiffTestR_example_study_design$condition == "A"), 4)
  expect_equal(sum(DiffTestR_example_study_design$condition == "B"), 4)
})

test_that("Example datasets are consistent with each other", {
  data(DiffTestR_example_data_long, envir = environment())
  data(DiffTestR_example_data_wide, envir = environment())
  data(DiffTestR_example_study_design, envir = environment())

  # Same number of unique precursors
  expect_equal(length(unique(DiffTestR_example_data_long$Precursor.Id)),
               nrow(DiffTestR_example_data_wide))

  # Same number of unique proteins
  expect_equal(length(unique(DiffTestR_example_data_long$Protein.Group)),
               length(unique(DiffTestR_example_data_wide$Protein.Group)))

  # Filenames in study design should match runs in long data
  expect_true(all(DiffTestR_example_study_design$filename %in%
                  unique(DiffTestR_example_data_long$Run)))

  # Column names in wide data (excluding first 3) should match study design filenames
  wide_sample_cols <- names(DiffTestR_example_data_wide)[4:11]
  expect_true(all(wide_sample_cols %in% DiffTestR_example_study_design$filename))
})

test_that("Example data contains expected species", {
  data(DiffTestR_example_data_long, envir = environment())

  # Should contain proteins from all three species
  has_human <- any(grepl("_HUMAN", DiffTestR_example_data_long$Protein.Names))
  has_yeast <- any(grepl("_YEAST", DiffTestR_example_data_long$Protein.Names))
  has_ecoli <- any(grepl("_ECOLI", DiffTestR_example_data_long$Protein.Names))

  expect_true(has_human)
  expect_true(has_yeast)
  expect_true(has_ecoli)

  # Count proteins per species
  proteins_per_species <- DiffTestR_example_data_long[
    !grepl(";", Protein.Group),
    .(n_proteins = length(unique(Protein.Group))),
    by = .(species = ifelse(grepl("_HUMAN", Protein.Names), "Human",
           ifelse(grepl("_YEAST", Protein.Names), "Yeast",
           ifelse(grepl("_ECOLI", Protein.Names), "Ecoli", "Other"))))
  ]

  # Should have roughly equal representation (100 each)
  expect_true(all(proteins_per_species$n_proteins >= 90))
})

test_that("Example data has appropriate data quality", {
  data(DiffTestR_example_data_long, envir = environment())

  # Quantification values should be positive
  expect_true(all(DiffTestR_example_data_long$Precursor.Quantity >= 0,
                  na.rm = TRUE))

  # Should have some missing values (realistic data)
  data(DiffTestR_example_data_wide, envir = environment())
  has_missing <- any(is.na(DiffTestR_example_data_wide[, 4:11]))
  expect_true(has_missing)

  # But not all missing
  all_missing <- all(is.na(DiffTestR_example_data_wide[, 4:11]))
  expect_false(all_missing)
})
