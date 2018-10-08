# Purpose: Function to import genotypes from plink 
# Updated: 180918

#' Read Genotypes from Plink
#' 
#' Reads compressed plink genotypes into R. Specify the \code{stem} of the input
#' files, excluding extensions. The user can choose to retain only certain loci 
#' from the BIM file, or subjects from the FAM file. Genotypes are formatted
#' into a numeric matrix, with subjects as rows and loci as columns. If \code{keep=T},
#' the BIM and FAM files for selected loci and subjects are additionally returned. 
#' If \code{parallel=T}, the genotypes are imported in parallel. 
#'
#' @param stem Stem of the BED, BIM, FAM files.
#' @param loci Row numbers of desired loci in the BIM file. If omitted, all loci
#'   are returned.
#' @param subj Row numbers of desired subjects in the FAM files. If omitted, all
#'   subjects are turned.
#' @param keep Return BIM and FAM files in addition to genotypes?
#' @param parallel Run in parallel? Must register parallel backend first. 
#' 
#' @return Either a numeric genotype matrix, or a list containing the genotype
#'   matrix together with the imported BIM and FAM files.
#'
#' @importFrom data.table fread
#' @importFrom foreach foreach '%do%' '%dopar%' registerDoSEQ
#' @export

ReadGeno = function(stem,loci=NULL,subj=NULL,keep=F,parallel=F){
  # Plink files
  BED = paste0(stem,".bed");
  BIM = paste0(stem,".bim");
  FAM = paste0(stem,".fam");
  # Input check
  if(!file.exists(BED)){stop("BED file DNE.")};
  if(!file.exists(BIM)){stop("BIM file DNE.")};
  if(!file.exists(FAM)){stop("FAM file DNE.")};
  # Import BIM
  BIM = fread(file=BIM,sep="\t",data.table=F,header=F);
  colnames(BIM) = c("Chr","Var","Pos","Coord","Minor","Major");
  snps = nrow(BIM);
  # Loci to import
  if(is.null(loci)){loci = seq(1:snps)};
  L = length(loci);
  # Import FAM
  FAM = fread(file=FAM,sep="\ ",data.table=F,header=F);
  colnames(FAM) = c("FID","IID","Father","Mother","Sex","Phenotype");
  obs = nrow(FAM);
  # Subjects to retain
  if(is.null(subj)){subj=rep(T,obs)} else {subj=(seq(1:obs)%in%subj)};
  # Import genotypes
  if(!parallel){registerDoSEQ()};
  j = NULL;
  G = foreach(j=1:L,.combine=cbind,.inorder=T) %dopar% {
    # Import
    g = readbed(bed=BED,obs=obs,snp=loci[j]);
    # Select subjects
    g = g[subj,,drop=F];
    return(g);
  }
  # Formatting
  rownames(G) = FAM[subj,"IID"];
  colnames(G) = BIM[loci,"Var"];
  # Output
  if(!keep){
    Out = G;
  } else {
    Out = list("G"=G,"BIM"=BIM[loci,,drop=F],"FAM"=FAM[subj,,drop=F]);
  }
  return(Out);
}