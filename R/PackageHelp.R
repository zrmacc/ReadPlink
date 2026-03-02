# Package documentation
# Updated: 2023-07-06

## usethis namespace: start
#' @useDynLib ReadPlink, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' ReadPlink: Read PLINK binary genotype files
#'
#' Read compressed PLINK binary genotype files (BED/BIM/FAM) into R.
#' Provides \code{\link{ReadGeno}} for reading genotypes into a matrix,
#' \code{\link{ReadPlinkData}} for genotypes plus BIM/FAM metadata,
#' \code{\link{readbed}} for low-level single-SNP access, and
#' \code{\link{SamplePath}} for the path to bundled sample data.
#'
#' @author Zachary R. McCaw
#' @keywords internal
"_PACKAGE"
