test_that("ReadGeno errors when files are missing", {
  tmp <- tempfile()
  expect_error(ReadGeno(tmp), "not found")
})

test_that("ReadGeno errors when only some files exist", {
  tmp <- tempdir()
  stem <- file.path(tmp, "partial")
  writeLines("1\ts1\t0\t0\t1\t-9", paste0(stem, ".fam"))
  expect_error(ReadGeno(stem), "BED file not found")
  unlink(paste0(stem, ".fam"))
})

test_that("ReadGeno reads minimal fixture and returns correct structure", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  geno <- ReadGeno(fix$stem)

  expect_equal(dim(geno), c(fix$n_subj, fix$n_snp))
  expect_equal(rownames(geno), c("s1", "s2", "s3"))
  expect_equal(colnames(geno), c("snp1", "snp2"))
})

test_that("ReadGeno returns correct genotype values for fixture", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  geno <- ReadGeno(fix$stem)

  # SNP1: 0, 1, 2; SNP2: 1, 1, 0
  expect_true(geno["s1", "snp1"] == 0 && geno["s2", "snp1"] == 1 && geno["s3", "snp1"] == 2)
  expect_true(geno["s1", "snp2"] == 1 && geno["s2", "snp2"] == 1 && geno["s3", "snp2"] == 0)
})

test_that("ReadGeno respects bim_rows and fam_rows", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  geno_all <- ReadGeno(fix$stem)
  geno_sub <- ReadGeno(fix$stem, bim_rows = 2L, fam_rows = c(1L, 3L))

  expect_equal(dim(geno_sub), c(2L, 1L))
  expect_true(geno_sub["s1", "snp2"] == 1 && geno_sub["s3", "snp2"] == 0)
  expect_equal(geno_sub, geno_all[c("s1", "s3"), "snp2", drop = FALSE])
})

test_that("ReadGeno bim_rows returns subset of SNPs in requested order", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  geno_all <- ReadGeno(fix$stem)
  geno_snp2_only <- ReadGeno(fix$stem, bim_rows = 2L)
  expect_equal(geno_snp2_only, geno_all[, "snp2", drop = FALSE])

  geno_reversed <- ReadGeno(fix$stem, bim_rows = c(2L, 1L))
  expect_equal(colnames(geno_reversed), c("snp2", "snp1"))
  expect_equal(geno_reversed[, "snp1"], geno_all[, "snp1"])
  expect_equal(geno_reversed[, "snp2"], geno_all[, "snp2"])
})

test_that("ReadGeno fam_rows returns subset of subjects in requested order", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  geno_all <- ReadGeno(fix$stem)
  geno_subj_1_3 <- ReadGeno(fix$stem, fam_rows = c(1L, 3L))
  expect_equal(geno_subj_1_3, geno_all[c("s1", "s3"), ])

  geno_reorder <- ReadGeno(fix$stem, fam_rows = c(3L, 1L))
  expect_equal(rownames(geno_reorder), c("s3", "s1"))
  expect_equal(geno_reorder["s1", ], geno_all["s1", ])
  expect_equal(geno_reorder["s3", ], geno_all["s3", ])
})

test_that("ReadGeno rejects out-of-range bim_rows and fam_rows", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  expect_error(ReadGeno(fix$stem, bim_rows = 0L), "bim_rows")
  expect_error(ReadGeno(fix$stem, bim_rows = 10L), "bim_rows")
  expect_error(ReadGeno(fix$stem, fam_rows = 0L), "fam_rows")
  expect_error(ReadGeno(fix$stem, fam_rows = 10L), "fam_rows")
})
