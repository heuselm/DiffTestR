#' qdata: Proteomics Quantification Data S3 Class
#'
#' @description
#' The `qdata` class is an S3 class for storing and managing proteomics
#' quantification data in DiffTestR. It provides a standardized structure
#' for precursor-level quantification matrices along with their annotations.
#'
#' @section Components:
#' A qdata object is a list with the following components:
#'
#' \describe{
#'   \item{data_long}{Data.table in long format with columns: Protein.Group,
#'     Precursor.Id, File.Name, Precursor.Quantity, condition, replicate}
#'   \item{data_wide}{Data.table in wide format with Protein.Group and
#'     Precursor.Id as first columns, followed by sample columns}
#'   \item{matrix_raw}{Numeric matrix of raw intensities with precursors
#'     in rows and samples in columns}
#'   \item{ann_col}{Data.table of sample annotations with columns: filename,
#'     condition, replicate (and optional additional columns)}
#'   \item{ann_row}{Data.table of precursor annotations with columns:
#'     Protein.Group, Precursor.Id (and optional: Protein.Ids, Protein.Names,
#'     Genes, Stripped.Sequence, Modified.Sequence, Precursor.Charge, etc.)}
#'   \item{study_design}{Data.table containing the original study design}
#'   \item{source_file}{Character string with path to the source data file}
#' }
#'
#' @section Creating qdata objects:
#' qdata objects are typically created using import functions:
#' \itemize{
#'   \item \code{\link{importFromDIANN}} - Import from DIA-NN long format
#'   \item \code{\link{importFromTables}} - Import from wide format tables
#' }
#'
#' @section Working with qdata objects:
#' Many DiffTestR functions now accept qdata objects directly:
#' \itemize{
#'   \item \code{\link{assessDataQuality}} - Quality assessment
#'   \item More functions to be updated...
#' }
#'
#' @section S3 Methods:
#' \itemize{
#'   \item \code{print.qdata} - Print summary information
#'   \item \code{summary.qdata} - Detailed summary statistics
#'   \item \code{validate_qdata} - Validate object structure
#' }
#'
#' @section Design Rationale:
#' The qdata class provides several benefits:
#' \itemize{
#'   \item \strong{Type safety}: Functions can validate input using inherits()
#'   \item \strong{Convenience}: No need to pass matrix + study_design separately
#'   \item \strong{Consistency}: Standardized structure across all functions
#'   \item \strong{Extensibility}: Can add methods for plotting, subsetting, etc.
#'   \item \strong{Backwards compatibility}: Functions can still accept
#'     matrix + study_design for legacy code
#' }
#'
#' @examples
#' \dontrun{
#' # Import data to create qdata object
#' qdata <- importFromTables(
#'   path_to_tsv = "data.tsv",
#'   study_design = "design.txt"
#' )
#'
#' # Print summary
#' print(qdata)
#'
#' # Detailed summary
#' summary(qdata)
#'
#' # Validate structure
#' validate_qdata(qdata)
#'
#' # Use with downstream functions
#' qc <- assessDataQuality(qdata)
#'
#' # Access components
#' head(qdata$data_wide)
#' qdata$matrix_raw[1:5, 1:3]
#' qdata$ann_row[1:5, ]
#' }
#'
#' @name qdata
#' @aliases qdata-class
NULL
