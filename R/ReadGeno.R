# Updated: 2023-07-06

#' Read Genotypes from Plink
#' 
#' Reads compressed plink genotypes into R. Specify the \code{stem} of the input
#' files, excluding extensions. The user can choose to retain only certain loci 
#' from the BIM file, or subjects from the FAM file. Genotypes are formatted
#' into a numeric matrix, with subjects as rows and loci as columns.
#'
#' @param stem Stem of the BED, BIM, FAM files.
#' @param bim_rows Row numbers of desired SNPs in the BIM file. If omitted, all loci
#'   are returned.
#' @param fam_rows Row numbers of desired subjects in the FAM files. If omitted, all
#'   subjects are turned.
#' @return Numeric subject by SNP genotype matrix. 
#' @export
ReadGeno <- function(
  stem, 
  bim_rows = NULL, 
  fam_rows = NULL
){
  
  # PLINK file names.
  bed_file <- paste0(stem, ".bed")
  bim_file <- paste0(stem, ".bim")
  fam_file <- paste0(stem, ".fam")
  
  # Check PLINK files all exist. 
  if (!file.exists(bed_file)) {stop("BED file DNE.")}
  if (!file.exists(bim_file)) {stop("BIM file DNE.")}
  if (!file.exists(fam_file)) {stop("FAM file DNE.")}
  
  # Import BIM.
  bim_data <- data.table::fread(
    file = bim_file,
    header = FALSE,
    select = 2,
    col.names = "snp"
  )
  
  if (is.null(bim_rows)) {
    bim_rows <- seq_len(nrow(bim_data))
  } 
  n_snps <- length(bim_rows)
  
  # Import FAM file.
  fam_data <- data.table::fread(
    file = fam_file,
    header = FALSE,
    select = 2,
    col.names = "id"
  )
  
  if (is.null(fam_rows)) {
    fam_rows <- seq_len(nrow(fam_data))
  }
  total_subj <- nrow(fam_data)
  
  # Genotypes.
  geno_matrix <- lapply(bim_rows, function(r) {
    genotypes <- readbed(bed = bed_file, obs = total_subj , snp = r)
    genotypes <- genotypes[fam_rows, , drop = FALSE]
  })
  geno_matrix <- do.call(cbind, geno_matrix)
  
  # Output.
  dimnames(geno_matrix) <- list(fam_data$id[fam_rows], bim_data$snp[bim_rows])
  return(geno_matrix)
}

