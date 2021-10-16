test_that("Incorrect function argument type returns error", {
  expect_error(gg_geneplot("1", 170054349 - 1e6, 170054349 + 1e6, "GRCh37"))
  expect_error(gg_geneplot(1, "170054349 - 1e6", 170054349 + 1e6, "GRCh37"))
  expect_error(gg_geneplot(1, 170054349 - 1e6, "170054349 + 1e6", "GRCh37"))
  expect_error(gg_geneplot(1, 170054349 - 1e6, 170054349 + 1e6, "X"))
})

test_that("gg_geneplot returns a ggplot object", {
  # GRCh37
  geneplot_res <- gg_geneplot(1, 170054349 - 1e6, 170054349 + 1e6, "GRCh37")
  expect_s3_class(geneplot_res, "ggplot")

  # GRCh38
  geneplot_res <- gg_geneplot(1, 170054349 - 1e6, 170054349 + 1e6, "GRCh38")
  expect_s3_class(geneplot_res, "ggplot")
})
