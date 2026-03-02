test_that("ReadPlinkData returns list with geno, bim, fam", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  out <- ReadPlinkData(fix$stem)
  expect_named(out, c("geno", "bim", "fam"))
  expect_true(is.matrix(out$geno))
  expect_true(is.data.frame(out$bim))
  expect_true(is.data.frame(out$fam))
})

test_that("ReadPlinkData geno matches ReadGeno", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  out <- ReadPlinkData(fix$stem)
  geno_direct <- ReadGeno(fix$stem)
  expect_equal(out$geno, geno_direct)
})

test_that("ReadPlinkData bim has expected columns and rows", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  out <- ReadPlinkData(fix$stem)
  expect_equal(nrow(out$bim), fix$n_snp)
  expect_equal(
    names(out$bim),
    c("chromosome", "id", "cm", "pos", "a1", "a2")
  )
  expect_equal(out$bim$id, c("snp1", "snp2"))
})

test_that("ReadPlinkData fam has expected columns and rows", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  out <- ReadPlinkData(fix$stem)
  expect_equal(nrow(out$fam), fix$n_subj)
  expect_equal(
    names(out$fam),
    c("fid", "iid", "pid", "mid", "sex", "pheno")
  )
  expect_equal(out$fam$iid, c("s1", "s2", "s3"))
})

test_that("ReadPlinkData respects bim_rows and fam_rows", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  out <- ReadPlinkData(fix$stem, bim_rows = 2L, fam_rows = c(1L, 3L))
  expect_equal(dim(out$geno), c(2L, 1L))
  expect_equal(nrow(out$bim), 1L)
  expect_equal(out$bim$id, "snp2")
  expect_equal(nrow(out$fam), 2L)
  expect_equal(out$fam$iid, c("s1", "s3"))
})

test_that("ReadPlinkData errors when files missing", {
  tmp <- tempfile()
  expect_error(ReadPlinkData(tmp), "not found")
})
