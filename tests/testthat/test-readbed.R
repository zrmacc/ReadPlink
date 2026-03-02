test_that("readbed errors on non-BED file", {
  tmp <- tempfile()
  writeLines("not a bed file", tmp)
  on.exit(unlink(tmp))
  expect_error(readbed(tmp, obs = 3, snp = 1), "not a PLINK")
})

test_that("readbed errors on wrong allele order", {
  # Valid magic but individual-major (0x00) instead of SNP-major (0x01)
  tmp <- tempfile()
  writeBin(as.raw(c(0x6C, 0x1B, 0x00)), tmp)
  on.exit(unlink(tmp))
  expect_error(readbed(tmp, obs = 3, snp = 1), "SNP-major")
})

test_that("readbed returns correct length and values for fixture", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  v1 <- readbed(fix$bed, obs = fix$n_subj, snp = 1)
  v2 <- readbed(fix$bed, obs = fix$n_subj, snp = 2)

  expect_length(v1, fix$n_subj)
  expect_length(v2, fix$n_subj)
  expect_true(all(as.numeric(v1) == c(0, 1, 2)))
  expect_true(all(as.numeric(v2) == c(1, 1, 0)))
})

test_that("readbed SNP index is 1-based", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  v_first <- readbed(fix$bed, obs = fix$n_subj, snp = 1)
  v_second <- readbed(fix$bed, obs = fix$n_subj, snp = 2)
  expect_true(all(as.numeric(v_first) == c(0, 1, 2)))
  expect_true(all(as.numeric(v_second) == c(1, 1, 0)))
})

test_that("readbed_multi matches readbed results", {
  fix <- create_plink_fixture()
  on.exit({
    unlink(fix$bed)
    unlink(fix$bim)
    unlink(fix$fam)
  })

  by_multi <- ReadPlink:::readbed_multi(
    fix$bed,
    obs = fix$n_subj,
    snp_indices = c(1L, 2L)
  )
  by_single <- cbind(
    readbed(fix$bed, obs = fix$n_subj, snp = 1),
    readbed(fix$bed, obs = fix$n_subj, snp = 2)
  )
  expect_equal(dim(by_multi), c(fix$n_subj, 2L))
  expect_true(all(as.numeric(by_multi) == as.numeric(by_single)))
})
