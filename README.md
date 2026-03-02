# ReadPlink

Zachary R. McCaw <br> Updated: 2026-03-02

<!-- badges: start -->

[![R-CMD-check](https://github.com/zrmacc/ReadPlink/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zrmacc/ReadPlink/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Read PLINK binary genotype files (BED/BIM/FAM) into R.

## Installation

``` r
devtools::install_github("yourusername/ReadPlink")
```

## Usage

Use the shared file stem (path without `.bed`/`.bim`/`.fam`) to read
genotypes into a subject-by-SNP matrix. The package contains a small
sample in `inst/plink_data` (installed as `plink_data/`); use
`SamplePath()` to get its path:

``` r
library(ReadPlink)

# Load sample plink_data (3 individuals, 4 SNPs)
stem <- SamplePath()
if (nchar(stem) > 0L) {
  geno <- ReadGeno(stem)
}
print(geno)
```

    ##      snp1 snp2 snp3 snp4
    ## ind1    0    1    2    0
    ## ind2    1    1    0    0
    ## ind3    2    0    1    2

General pattern:

``` r
geno <- ReadGeno("/path/to/plink_stem")

# Optional: subset by BIM row indices (SNPs) and FAM row indices (subjects)
geno <- ReadGeno("/path/to/plink_stem", bim_rows = 1:100, fam_rows = 1:50)
```
