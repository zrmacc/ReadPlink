# Path to bundled sample PLINK data

#' Path to bundled sample PLINK data
#'
#' Returns the path stem (no extension) for the sample BED/BIM/FAM set shipped
#' with the package (3 individuals, 4 SNPs). Use with \code{\link{ReadGeno}}
#' or \code{\link{readbed}}.
#'
#' @return Character path to the sample stem (e.g. \code{.../plink_data/sample}).
#'   Returns \code{""} if the package is not installed.
#' @export
#'
#' @examples
#' stem <- SamplePath()
#' if (nchar(stem) > 0L) {
#'   geno <- ReadGeno(stem)
#'   dim(geno)
#' }
SamplePath <- function() {
  bed_path <- system.file("plink_data", "sample.bed", package = "ReadPlink")
  if (bed_path == "") {
    return("")
  }
  return(sub("\\.bed$", "", bed_path))
}
