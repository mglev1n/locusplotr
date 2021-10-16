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

test_that("ld_extract_locuszoom returns a tibble, or error if variant is missing", {
  # Appropriate ref/alt
  ld_df <- ld_extract_locuszoom(chrom = 16, pos = 53830055, ref = "C", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "GRCh37", population = "ALL", metric = "rsquare")
  expect_s3_class(ld_df, "tbl_df")
  # Flipped variants
  ld_df <- ld_extract_locuszoom(chrom = 16, pos = 53830055, ref = "G", alt = "C", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "GRCh37", population = "ALL", metric = "rsquare")
  expect_s3_class(ld_df, "tbl_df")
  # Variant missing from LD panel gives error
  expect_error(ld_extract_locuszoom(chrom = 16, pos = 54065268, ref = "C", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "GRCh37", population = "ALL", metric = "rsquare"), regexp = "Reference variant is missing from the specified LD panel")
})
