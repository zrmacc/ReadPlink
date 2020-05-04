# Purpose: Function to import genotypes from plink 
# Updated: 20/05/03

#' Read Genotypes from Plink
#' 
#' Reads compressed plink genotypes into R. Specify the \code{stem} of the input
#' files, excluding extensions. The user can choose to retain only certain loci 
#' from the BIM file, or subjects from the FAM file. Genotypes are formatted
#' into a numeric matrix, with subjects as rows and loci as columns.
#'
#' @param stem Stem of the BED, BIM, FAM files.
#' @param snp_rows Row numbers of desired SNPs in the BIM file. If omitted, all loci
#'   are returned.
#' @param subj_rows Row numbers of desired subjects in the FAM files. If omitted, all
#'   subjects are turned.
#' @param label_geno_cols Label columns of genotype matrix with SNP names?
#' @param label_geno_rows Label rows of genotype matrix with subject IIDs?
#' 
#' @return Number subject by SNP genotype matrix. 
#'
#' @importFrom data.table fread
#' @export

ReadGeno <- function(stem, 
                     snp_rows = NULL, 
                     subj_rows = NULL, 
                     label_geno_cols = FALSE,
                     label_geno_rows = FALSE){
  
  # PLINK file names.
  BED <- paste0(stem, ".bed")
  BIM <- paste0(stem, ".bim")
  FAM <- paste0(stem, ".fam")
  
  # Check PLINK files all exist. 
  if(!file.exists(BED)){stop("BED file DNE.")}
  if(!file.exists(BIM)){stop("BIM file DNE.")}
  if(!file.exists(FAM)){stop("FAM file DNE.")}
  
  # Import BIM, if necessary.
  if(is.null(snp_rows)){
    bim <- fread(file = BIM, sep = "\t", header = FALSE, select = c(2))
    n_snps <- nrow(bim)
    snp_rows <- seq(1:n_snps)
  } else {
    snp_rows <- unique(snp_rows)
    n_snps <- length(snp_rows)
  }
  
  # Import FAM file, if necessary.
  if(is.null(subj_rows)){
    fam <- fread(file = FAM, sep = "\ ", header = FALSE, select = c(2))
    n_subj <- nrow(fam)
    subj_rows <- seq(1:n_subj)
  } else {
    subj_rows <- unique(subj_rows)
    # Get total number of subjects.
    sys_call <- paste0('wc -l ', FAM)
    sys_output <- trimws(system(sys_call, intern = TRUE))
    n_subj <- as.integer(gsub(pattern = "([0-9]+).*", replacement = "\\1", x = sys_output))
  }
  
  # Genotypes
  aux <- function(j){
    genotypes <- readbed(bed = BED, obs = n_subj , snp = snp_rows[j])
    # Select subjects
    genotypes <- genotypes[subj_rows, , drop = FALSE]
  }
  
  geno_matrix <- lapply(seq(1:n_snps), aux)
  geno_matrix <- do.call(cbind, geno_matrix)
  
  # Apply col names, if requested. 
  if(label_geno_cols){
    if(!exists("bim")){
      bim <- fread(file = BIM, sep = "\t", header = FALSE, select = c(2))
    }
    colnames(geno_matrix) <- bim[[1]][snp_rows]
  }
  
  # Apply row names, if requested.
  if(label_geno_rows){
    if(!exists("fam")){
      fam <- fread(file = FAM, sep = "\ ", header = FALSE, select = c(2))
    }
    rownames(geno_matrix) <- fam[[1]][subj_rows]
  }
  
  # Output
  return(geno_matrix)
}