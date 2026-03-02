# Updated: 2026-03-02

#' Read Genotypes from Plink
#'
#' Reads compressed PLINK genotypes into R. Specify the \code{stem} of the input
#' files, excluding extensions. The user can choose to retain only certain loci
#' from the BIM file, or subjects from the FAM file. Genotypes are formatted
#' into a numeric matrix, with subjects as rows and loci as columns.
#'
#' @param stem Stem of the BED, BIM, FAM files (path without extension).
#' @param bim_rows Row numbers of desired SNPs in the BIM file. If omitted, all
#'   loci are returned.
#' @param fam_rows Row numbers of desired subjects in the FAM file. If omitted,
#'   all subjects are returned.
#' @return Numeric subject-by-SNP genotype matrix (subjects as rows, SNPs as
#'   columns).
#' @export
ReadGeno <- function(stem, bim_rows = NULL, fam_rows = NULL) {
  bed_file <- paste0(stem, ".bed")
  bim_file <- paste0(stem, ".bim")
  fam_file <- paste0(stem, ".fam")

  if (!file.exists(bed_file)) {
    stop("BED file not found.")
  }
  if (!file.exists(bim_file)) {
    stop("BIM file not found.")
  }
  if (!file.exists(fam_file)) {
    stop("FAM file not found.")
  }

  bim_data <- data.table::fread(
    file = bim_file,
    header = FALSE,
    select = 2,
    col.names = "snp"
  )
  n_snps <- nrow(bim_data)

  if (is.null(bim_rows)) {
    bim_rows <- seq_len(n_snps)
  } else {
    bim_rows <- as.integer(bim_rows)
    if (any(is.na(bim_rows)) || any(bim_rows < 1L) || any(bim_rows > n_snps)) {
      stop("bim_rows must be integers between 1 and the number of SNPs in the BIM file.")
    }
  }

  fam_data <- data.table::fread(
    file = fam_file,
    header = FALSE,
    select = 2,
    col.names = "id"
  )
  n_subj <- nrow(fam_data)

  if (is.null(fam_rows)) {
    fam_rows <- seq_len(n_subj)
  } else {
    fam_rows <- as.integer(fam_rows)
    if (any(is.na(fam_rows)) || any(fam_rows < 1L) || any(fam_rows > n_subj)) {
      stop("fam_rows must be integers between 1 and the number of subjects in the FAM file.")
    }
  }
  total_subj <- n_subj

  geno_matrix <- readbed_multi(
    bed = bed_file,
    obs = total_subj,
    snp_indices = bim_rows
  )
  geno_matrix <- geno_matrix[fam_rows, , drop = FALSE]
  dimnames(geno_matrix) <- list(fam_data$id[fam_rows], bim_data$snp[bim_rows])
  return(geno_matrix)
}

