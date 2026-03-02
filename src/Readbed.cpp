// Purpose: Read PLINK .bed (binary genotype) files.
// Updated: 2026-03-02

#include <Rcpp.h>
#include <cmath>
#include <fstream>
#include <vector>

// Decode one byte into at most 4 genotypes; write into geno starting at idx.
// PLINK encoding: 2 bits per subject, (b0,b1): 00->2, 01->1, 10->NA, 11->0.
inline void decode_byte(unsigned char b, int n_obs, int idx, Rcpp::NumericVector& geno) {
  for (int k = 0; k < 4 && idx + k < n_obs; k++) {
    int b0 = (b >> (2 * k)) & 1;
    int b1 = (b >> (2 * k + 1)) & 1;
    if (b0 == 1 && b1 == 0) {
      geno[idx + k] = NA_REAL;
    } else {
      geno[idx + k] = 2.0 - (b0 + b1);
    }
  }
}

//' Read BED Files
//'
//' Low-level reader for a single SNP from a PLINK binary BED file (SNP-major
//' mode).
//'
//' @param bed Path to the BED file.
//' @param obs Total number of subjects (rows) in the FAM file.
//' @param snp 1-based index of the SNP to read.
//' @return Numeric vector of genotypes (0, 1, 2, or NA) at the selected SNP.
//' @export
// [[Rcpp::export]]
SEXP readbed(const char* bed, const double obs, const int snp) {
  std::ifstream is(bed, std::ios::binary);
  const int n_obs = static_cast<int>(obs);
  const int bps = static_cast<int>(std::ceil(obs / 4.0));

  char head[3];
  is.read(head, 3);
  if (head[0] != 0x6C || head[1] != 0x1B) {
    throw std::range_error("Input is not a PLINK .bed file.");
  }
  if (head[2] != 1) {
    throw std::range_error("Genotypes not coded in SNP-major order.");
  }

  const int idx0 = 3 + (snp - 1) * bps;
  std::vector<char> bytes(bps);
  is.seekg(idx0);
  is.read(bytes.data(), bps);
  is.close();

  Rcpp::NumericVector geno(n_obs);
  int idx = 0;
  for (int i = 0; i < bps && idx < n_obs; i++) {
    decode_byte(static_cast<unsigned char>(bytes[i]), n_obs, idx, geno);
    idx += 4;
  }
  return geno;
}

//' Read multiple SNPs from a BED file in one pass
//'
//' Reads requested SNPs from a PLINK binary BED file with a single file open.
//' Used internally by \code{\link{ReadGeno}} for efficiency.
//'
//' @param bed Path to the BED file.
//' @param obs Total number of subjects (rows) in the FAM file.
//' @param snp_indices 1-based indices of SNPs to read (integer vector).
//' @return Numeric matrix with \code{obs} rows and \code{length(snp_indices)}
//'   columns (one column per SNP).
//' @keywords internal
// [[Rcpp::export]]
SEXP readbed_multi(const char* bed, const double obs, const Rcpp::IntegerVector& snp_indices) {
  std::ifstream is(bed, std::ios::binary);
  const int n_obs = static_cast<int>(obs);
  const int bps = static_cast<int>(std::ceil(obs / 4.0));
  const int n_snps = snp_indices.size();

  char head[3];
  is.read(head, 3);
  if (head[0] != 0x6C || head[1] != 0x1B) {
    throw std::range_error("Input is not a PLINK .bed file.");
  }
  if (head[2] != 1) {
    throw std::range_error("Genotypes not coded in SNP-major order.");
  }

  Rcpp::NumericMatrix geno(n_obs, n_snps);
  std::vector<char> bytes(bps);

  for (int j = 0; j < n_snps; j++) {
    int snp = snp_indices[j];
    const int idx0 = 3 + (snp - 1) * bps;
    is.seekg(idx0);
    is.read(bytes.data(), bps);
    int row = 0;
    for (int i = 0; i < bps && row < n_obs; i++) {
      unsigned char b = static_cast<unsigned char>(bytes[i]);
      for (int k = 0; k < 4 && row < n_obs; k++) {
        int b0 = (b >> (2 * k)) & 1;
        int b1 = (b >> (2 * k + 1)) & 1;
        geno(row, j) = (b0 == 1 && b1 == 0) ? NA_REAL : (2.0 - (b0 + b1));
        row++;
      }
    }
  }
  is.close();
  return geno;
}
