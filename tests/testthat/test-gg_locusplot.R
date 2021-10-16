test_that("Incorrect function argument type returns error", {
  expect_error(gg_locusplot(df = "fto_locus_df", lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_distance = 1e6, path = NULL))
  expect_error(gg_locusplot(df = fto_locus_df, lead_snp = 1, rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL))
  expect_error(gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_title = NULL, plot_subtitle = NULL, plot_distance = "X", path = NULL))
  expect_error(gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL, plot_genes = "X"))
  })

test_that("Default gg_locusplot returns a ggplot object", {
  locusplot_res <- gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res, "ggplot")
})

test_that("gg_locusplot with plot_genes returns a patchwork object", {
  locusplot_res <- gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = TRUE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res, "patchwork")
})
