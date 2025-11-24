#' ddata: Differential Expression Data S3 Class
#'
#' @description
#' The `ddata` class is an S3 class for storing and managing differential
#' expression analysis results in DiffTestR. It provides a standardized structure
#' for differential abundance testing results along with the processed data matrices
#' and statistical summaries.
#'
#' @section Components:
#' A ddata object is a list with the following components:
#'
#' \describe{
#'   \item{data_source}{Character string with path to the source data file}
#'   \item{comparison}{Character string describing the comparison (e.g., "A/B")}
#'   \item{conditions}{Named list with two elements:
#'     \itemize{
#'       \item condition_1: Name of the first condition
#'       \item condition_2: Name of the second condition
#'     }
#'   }
#'   \item{diffExpr_result_dt}{Data.table containing the complete differential
#'     expression results with columns:
#'     \itemize{
#'       \item Protein.Group: Protein group identifier
#'       \item Precursor.Id: Precursor identifier
#'       \item Protein.Names: Protein names/descriptions
#'       \item log2_fold_change: Precursor-level log2 fold change
#'       \item p_value: Precursor-level p-value
#'       \item p_value_BHadj: Precursor-level BH-adjusted p-value
#'       \item log2_fold_change_protein: Protein-level log2 fold change
#'       \item p_value_protein: Protein-level p-value
#'       \item p_value_BHadj_protein: Protein-level BH-adjusted p-value
#'       \item n_precursors: Number of precursors per protein
#'       \item Plus raw and processed intensity columns
#'     }
#'   }
#'   \item{mat_quant_log2_qnorm_imp_minObs}{Numeric matrix of processed
#'     intensities (log2-transformed, quantile-normalized, imputed, and filtered
#'     by minimum observations)}
#'   \item{mat_quant_log2_qnorm}{Numeric matrix of log2-transformed and
#'     quantile-normalized intensities (before imputation and filtering)}
#'   \item{mat_quant_log2}{Numeric matrix of log2-transformed raw intensities}
#'   \item{mat_quant}{Numeric matrix of raw intensities}
#'   \item{study_design}{Data.table containing the study design with columns:
#'     filename, condition, replicate}
#'   \item{summary_stats}{List of summary statistics including:
#'     \itemize{
#'       \item n_precursors_input: Number of input precursors
#'       \item n_precursors_tested: Number of tested precursors
#'       \item n_proteins_tested: Number of tested proteins
#'       \item n_significant_precursors: Number of significant precursors
#'       \item n_significant_proteins: Number of significant proteins
#'       \item n_upregulated: Number of up-regulated proteins
#'       \item n_downregulated: Number of down-regulated proteins
#'       \item median_log2fc_precursor: Median precursor-level log2FC
#'       \item median_log2fc_protein: Median protein-level log2FC
#'       \item median_cv_condition1: Median CV for condition 1
#'       \item median_cv_condition2: Median CV for condition 2
#'     }
#'   }
#'   \item{candidates_condition1}{Data.table of proteins significantly higher
#'     in condition 1 (p <= 0.01 and log2FC >= 1)}
#'   \item{candidates_condition2}{Data.table of proteins significantly higher
#'     in condition 2 (p <= 0.01 and log2FC <= -1)}
#' }
#'
#' @section Creating ddata objects:
#' ddata objects are created by the differential abundance testing function:
#' \itemize{
#'   \item \code{\link{testDifferentialAbundance}} - Main differential testing function
#' }
#'
#' The typical workflow is:
#' \enumerate{
#'   \item Import data using \code{\link{importFromDIANN}} or \code{\link{importFromTables}}
#'     to create a qdata object
#'   \item Pass the qdata object to \code{\link{testDifferentialAbundance}}
#'     to perform statistical testing and create a ddata object
#'   \item Use the ddata S3 methods (print, summary, plot) to explore results
#' }
#'
#' @section S3 Methods:
#' \itemize{
#'   \item \code{print.ddata} - Print concise summary information
#'   \item \code{summary.ddata} - Detailed summary statistics and top candidates
#'   \item \code{plot.ddata} - Create volcano plots or MA plots
#'     \itemize{
#'       \item Static ggplot2 plots (default)
#'       \item Interactive plotly plots (with \code{interactive = TRUE})
#'       \item Volcano plot (default, \code{which = "volcano"})
#'       \item MA plot (\code{which = "ma"})
#'     }
#'   \item \code{validate_ddata} - Validate object structure
#' }
#'
#' @section Design Rationale:
#' The ddata class provides several benefits:
#' \itemize{
#'   \item \strong{Consistency}: Standardized structure for all differential
#'     expression results
#'   \item \strong{Self-contained}: All relevant data and results in one object
#'   \item \strong{Visualization}: Built-in plotting methods for volcano and MA plots
#'   \item \strong{Interactivity}: Easy conversion to interactive plotly visualizations
#'   \item \strong{Extensibility}: Can add new methods and components as needed
#'   \item \strong{Reproducibility}: Stores all processing parameters and
#'     intermediate matrices
#' }
#'
#' @section Comparison with qdata:
#' While \code{qdata} objects store raw quantification data, \code{ddata} objects
#' contain the results of statistical testing:
#' \itemize{
#'   \item \code{qdata}: Input data (precursor intensities, annotations, study design)
#'   \item \code{ddata}: Analysis results (statistical tests, fold changes, p-values)
#' }
#'
#' The two classes work together in the DiffTestR workflow:
#' \code{qdata -> testDifferentialAbundance() -> ddata}
#'
#' @examples
#' \dontrun{
#' # Complete workflow
#' qdata <- importFromTables(
#'   path_to_tsv = "data.tsv",
#'   study_design = "design.txt"
#' )
#'
#' ddata <- testDifferentialAbundance(
#'   data = qdata,
#'   condition_1 = "Control",
#'   condition_2 = "Treatment"
#' )
#'
#' # Print concise summary
#' print(ddata)
#'
#' # Detailed summary with top candidates
#' summary(ddata)
#'
#' # Validate structure
#' validate_ddata(ddata)
#'
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
#'
#' # Without labels
#' plot(ddata, label_top = 0)
#'
#' # More labels
#' plot(ddata, label_top = 20)
#'
#' # Access components
#' head(ddata$diffExpr_result_dt)
#' ddata$summary_stats
#' ddata$candidates_condition1
#'
#' # Filter for specific proteins
#' my_proteins <- ddata$diffExpr_result_dt[
#'   grepl("MYC", Protein.Names, ignore.case = TRUE)
#' ]
#' }
#'
#' @name ddata
#' @aliases ddata-class
NULL
