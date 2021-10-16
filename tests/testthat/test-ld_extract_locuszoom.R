test_that("Incorrect function argument type returns errors", {
  expect_error(ld_extract_locuszoom(chrom = "1", pos = 53830055, ref = "C", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "GRCh37", population = "ALL", metric = "rsquare"))
  expect_error(ld_extract_locuszoom(chrom = 1, pos = "53830055", ref = "C", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "GRCh37", population = "ALL", metric = "rsquare"))
  expect_error(ld_extract_locuszoom(chrom = 1, pos = 53830055, ref = "X", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "GRCh37", population = "ALL", metric = "rsquare"))
  expect_error(ld_extract_locuszoom(chrom = 1, pos = 53830055, ref = "C", alt = "X", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "GRCh37", population = "ALL", metric = "rsquare"))
  expect_error(ld_extract_locuszoom(chrom = 1, pos = 53830055, ref = "C", alt = "G", start = "53830055 - 5e5", stop = 53830055 + 5e5, genome_build = "GRCh37", population = "ALL", metric = "rsquare"))
  expect_error(ld_extract_locuszoom(chrom = 1, pos = 53830055, ref = "C", alt = "G", start = 53830055 - 5e5, stop = "53830055 + 5e5", genome_build = "GRCh37", population = "ALL", metric = "rsquare"))
  expect_error(ld_extract_locuszoom(chrom = 1, pos = 53830055, ref = "C", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "X", population = "ALL", metric = "rsquare"))
  expect_error(ld_extract_locuszoom(chrom = 1, pos = 53830055, ref = "C", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "GRCh37", population = "X", metric = "rsquare"))
  expect_error(ld_extract_locuszoom(chrom = 1, pos = 53830055, ref = "C", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "GRCh37", population = "ALL", metric = "X"))
})

test_that("ld_extract_locuszoom returns a tibble", {
  ld_df <- ld_extract_locuszoom(chrom = 16, pos = 53830055, ref = "C", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "GRCh37", population = "ALL", metric = "rsquare")
  expect_s3_class(ld_df, "tbl_df")
})

