# Helper to create minimal valid PLINK BED/BIM/FAM for tests.
# BED: PLINK magic 0x6C 0x1B 0x01 (SNP-major), then 2 bits per subject per SNP.
# For n subjects, ceil(n/4) bytes per SNP. Encoding: 00=0, 01=1, 10=NA, 11=2.

create_plink_fixture <- function(dir = tempdir(), stem = "plink_test") {
  # 3 individuals, 2 SNPs
  n_subj <- 3L
  n_snp <- 2L
  bps <- ceiling(n_subj / 4)  # 1 byte per SNP

  # BED: magic + SNP1 byte (0,1,2) + SNP2 byte (1,1,0)
  # C++ pairs bits (0,1),(2,3),(4,5); geno=2-(b0+b1); (1,0)=missing. 0=(1,1), 1=(0,1), 2=(0,0).
  # SNP1: 0,1,2 -> 0x0B; SNP2: 1,1,0 -> 0x3A (pairs (b0,b1): 0=(1,1), 1=(0,1), 2=(0,0))
  bed_path <- file.path(dir, paste0(stem, ".bed"))
  writeBin(
    as.raw(c(0x6C, 0x1B, 0x01, 0x0B, 0x3A)),
    bed_path
  )

  # BIM: chr id cm pos a1 a2 (tab, no header)
  bim_path <- file.path(dir, paste0(stem, ".bim"))
  writeLines(
    c("1\tsnp1\t0\t1000\tA\tG", "1\tsnp2\t0\t2000\tC\tT"),
    bim_path
  )

  # FAM: fid iid pid mid sex pheno (tab, no header)
  fam_path <- file.path(dir, paste0(stem, ".fam"))
  writeLines(
    c("fam1\ts1\t0\t0\t1\t-9", "fam1\ts2\t0\t0\t2\t-9", "fam1\ts3\t0\t0\t1\t-9"),
    fam_path
  )

  return(list(
    stem = file.path(dir, stem),
    bed = bed_path,
    bim = bim_path,
    fam = fam_path,
    n_subj = n_subj,
    n_snp = n_snp
  ))
}
