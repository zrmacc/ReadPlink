// Purpose: Function to read .bed files
// Updated: 181007

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// For ifstream
#include <fstream>
// For ceil 
#include <math.h> 

//' Read BED Files
//' 
//' @param bed Name of the bed file.
//' @param obs Number of subjects.
//' @param snp Locus to select.  
//' 
//' @return Numeric vector of genotypes at the select locus. 
//' 
// [[Rcpp::export]]

SEXP readbed(const char* bed, const double obs, const int snp){
  // Open file for reading
  std::ifstream is(bed,std::ios::binary);
  // Bytes per snp
  const int bps = ceil(obs/4);
  // Get header
  char head[3];
  is.read((char *)head,3);
  // Check header
  if(head[0]!=0x6C || head[1]!=0x1B){
    throw std::range_error("Input in not a plink .bed file.");
  };
  // Check allele coding
  if(head[2]!=1){
    throw std::range_error("Genotypes not coded in snp major.");
  }
  // First byte in snp
  const int idx0 = 3+(snp-1)*bps;
  // Last byte in snp
  // const int idx1 = 3+snp*bps-1;
  // Read bytes in snp
  char bytes[bps];
  // arma::mat all_bits(8,bps);
  // Note: bytes is a C-type array
  is.seekg(idx0);
  is.read((char *)bytes,bps);
  // Current individual
  int current_n;
  current_n = 0;
  // Declare current byte
  char current_byte;
  // Declare current bit string
  arma::vec current_bits(8);
  // Declare working bits
  int b0;
  int b1;
  // Declare genotype vector
  arma::vec geno((int)obs);
  // Loop over bytes
  for(int i=0; i<bps; i++){
    // Current byte
    current_byte = bytes[i];
    // Convert byte to bits
    for(int j=0; j<8; j++){
      current_bits(j) = ((current_byte&(1<<j)) != 0);
    };
    // all_bits.col(i) = current_bits;
    // Loop over subjects in byte
    for(int j=0; j<7; j+=2){
      // Current subject
      current_n += 1;
      // Store genotype if not all subjects have been read
      if(current_n<=obs){
        // Working bits
        b0 = current_bits(j);
        b1 = current_bits(j+1);
        // Check for missing value
        if((b0==1)&&(b1==0)){
          geno(current_n-1) = NA_REAL;
        } else {
          geno(current_n-1) = 2-(b0+b1);
        }; // End missingness check 
      }; // End observation index check 
    }; // End loop over bits 
  }; // End loop over bytes
  // Output
  return Rcpp::wrap(geno);
}
