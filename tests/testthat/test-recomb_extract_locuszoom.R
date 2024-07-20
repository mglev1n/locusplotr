test_that("recomb_locuszoom returns a tibble", {
  recomb_df <- recomb_extract_locuszoom(chrom = 1, start = 1000, end = 150000, genome_build = "GRCh38")
  expect_s3_class(recomb_df, "tbl_df")

  recomb_df <- recomb_extract_locuszoom(chrom = 1, start = 1000, end = 150000)
  expect_s3_class(recomb_df, "tbl_df")
})
