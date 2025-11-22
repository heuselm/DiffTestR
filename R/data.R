#' Example proteomics dataset (long format)
#'
#' A subset of 300 proteins from the Exploris480 HYE (Human-Yeast-E.coli)
#' species mix benchmark dataset. Contains precursor-level quantification data
#' for differential abundance testing between Mix A and Mix B (4 replicates each).
#'
#' @format A data.table with 18,822 rows and columns:
#' \describe{
#'   \item{Run}{MS run filename}
#'   \item{Protein.Group}{Protein group identifier}
#'   \item{Protein.Names}{Protein names with species suffix (_HUMAN, _YEAST, _ECOLI)}
#'   \item{Precursor.Id}{Unique precursor identifier}
#'   \item{Precursor.Quantity}{Precursor-level quantification value}
#' }
#'
#' @details
#' This dataset contains:
#' \itemize{
#'   \item 300 protein groups (100 from each species: H. sapiens, S. cerevisiae, E. coli)
#'   \item 2,835 unique precursors
#'   \item 8 MS runs (4 x Mix A, 4 x Mix B)
#'   \item Acquired on Thermo Exploris 480 Orbitrap
#' }
#'
#' The dataset is designed for testing differential abundance analysis functions
#' and demonstrating package workflows.
#'
#' @source Subset from HYE species mix benchmark (2023)
#' @seealso \code{\link{DiffTestR_example_data_wide}}, \code{\link{DiffTestR_example_study_design}}
"DiffTestR_example_data_long"

#' Example proteomics dataset (wide format)
#'
#' Wide format version of the example dataset with precursor quantities
#' arranged in columns by MS run.
#'
#' @format A data.table with 2,835 rows (precursors) and 11 columns:
#' \describe{
#'   \item{Protein.Group}{Protein group identifier}
#'   \item{Protein.Names}{Protein names with species suffix}
#'   \item{Precursor.Id}{Unique precursor identifier}
#'   \item{20230314_EXPL1_SA_Protoc_OT2_MixA_200ng_TR1_01}{Mix A replicate 1}
#'   \item{20230314_EXPL1_SA_Protoc_OT2_MixA_200ng_TR2_01}{Mix A replicate 2}
#'   \item{20230314_EXPL1_SA_Protoc_OT2_MixA_200ng_TR3_01}{Mix A replicate 3}
#'   \item{20230314_EXPL1_SA_Protoc_OT2_MixA_200ng_TR4_01}{Mix A replicate 4}
#'   \item{20230314_EXPL1_SA_Protoc_OT2_MixB_200ng_TR1_01}{Mix B replicate 1}
#'   \item{20230314_EXPL1_SA_Protoc_OT2_MixB_200ng_TR2_01}{Mix B replicate 2}
#'   \item{20230314_EXPL1_SA_Protoc_OT2_MixB_200ng_TR3_01}{Mix B replicate 3}
#'   \item{20230314_EXPL1_SA_Protoc_OT2_MixB_200ng_TR4_01}{Mix B replicate 4}
#' }
#'
#' @details
#' This is the wide format representation suitable for direct input to
#' \code{\link{testDifferentialAbundance}}.
#'
#' @source Subset from HYE species mix benchmark (2023)
#' @seealso \code{\link{DiffTestR_example_data_long}}, \code{\link{DiffTestR_example_study_design}}
"DiffTestR_example_data_wide"

#' Example study design
#'
#' Study design table describing the experimental design for the example dataset.
#'
#' @format A data.table with 8 rows and 3 columns:
#' \describe{
#'   \item{filename}{MS run filename (without .raw extension)}
#'   \item{condition}{Experimental condition ("A" or "B")}
#'   \item{replicate}{Biological replicate number (1-4)}
#' }
#'
#' @details
#' This table links each MS run to its experimental condition and replicate number.
#' Required for differential abundance testing with \code{\link{testDifferentialAbundance}}.
#'
#' @examples
#' data(DiffTestR_example_study_design)
#' print(DiffTestR_example_study_design)
#'
#' @seealso \code{\link{DiffTestR_example_data_wide}}, \code{\link{DiffTestR_example_data_long}}
"DiffTestR_example_study_design"
