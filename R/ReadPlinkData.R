# Convenience function returning geno matrix plus BIM and FAM metadata

#' Read PLINK data with genotype matrix and metadata
#'
#' Reads the BED/BIM/FAM set and returns the genotype matrix plus full BIM and
#' FAM tables (optionally subset by \code{bim_rows} and \code{fam_rows}). Useful
#' when you need SNP/sample metadata together with genotypes.
#'
#' @param stem Stem of the BED, BIM, FAM files (path without extension).
#' @param bim_rows Row numbers of desired SNPs in the BIM file. If \code{NULL},
#'   all SNPs are returned.
#' @param fam_rows Row numbers of desired subjects in the FAM file. If
#'   \code{NULL}, all subjects are returned.
#' @return A list with components \code{geno} (numeric subject-by-SNP matrix),
#'   \code{bim} (BIM data with columns chromosome, id, cm, pos, a1, a2), and
#'   \code{fam} (FAM data with columns fid, iid, pid, mid, sex, pheno).
#' @export
ReadPlinkData <- function(stem, bim_rows = NULL, fam_rows = NULL) {
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
    col.names = c("chromosome", "id", "cm", "pos", "a1", "a2")
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
    col.names = c("fid", "iid", "pid", "mid", "sex", "pheno")
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

  geno <- ReadGeno(stem = stem, bim_rows = bim_rows, fam_rows = fam_rows)
  bim_sub <- as.data.frame(bim_data[bim_rows, ])
  fam_sub <- as.data.frame(fam_data[fam_rows, ])
  return(list(geno = geno, bim = bim_sub, fam = fam_sub))
}
