test_that("Incorrect function argument type returns error", {
  expect_error(gg_locusplot(df = "fto_locus_df", lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_distance = 1e6, path = NULL))
  # expect_error(gg_locusplot(df = fto_locus_df, lead_snp = 1, rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL))
  expect_error(gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_title = NULL, plot_subtitle = NULL, plot_distance = "X", path = NULL))
  expect_error(gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL, plot_genes = "X"))
})

test_that("Default gg_locusplot returns a ggplot object", {
  locusplot_res <- gg_locusplot(df = fto_locus_df %>% mutate(trait = "BMI"), lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL, trait = trait)
  expect_s3_class(locusplot_res, "ggplot")
  # ggplot object returnd when lead_snp = NULL
  locusplot_res <- gg_locusplot(df = fto_locus_df, lead_snp = NULL, rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res, "ggplot")
})

test_that("Multitrait gg_locusplot returns a ggplot object", {
  locusplot_res <- gg_locusplot(df = fto_locus_df %>% mutate(trait = "BMI") %>% bind_rows(fto_locus_df %>% mutate(trait = "HbA1c")) %>% unique, rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL, trait = trait)
  expect_s3_class(locusplot_res, "ggplot")
  # ggplot object returned when lead_snp = NULL
  locusplot_res <- gg_locusplot(df = fto_locus_df, lead_snp = NULL, rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res, "ggplot")
})

test_that("Default gg_locusplot with plot_genes returns a patchwork object", {
  locusplot_res <- gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = TRUE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res, "patchwork")
})

test_that("gg_locusplot with missing LD returns error and ggplot object", {
  expect_message(gg_locusplot(df = fto_locus_df, lead_snp = "rs142090714", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL), regexp = "No linkage disequilibrium information found")
  locusplot_res <- gg_locusplot(df = fto_locus_df, lead_snp = "rs142090714", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res, "ggplot")
})

test_that("gg_locusplot with missing lead SNP returns error and ggplot object", {
  expect_message(gg_locusplot(df = fto_locus_df, lead_snp = "missing_snp", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL), regexp = "Lead snp not present in supplied locus data")
  locusplot_res <- gg_locusplot(df = fto_locus_df, lead_snp = "missing_snp", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res, "ggplot")
})

test_that("gg_locusplot saves to file", {
  .dir <- tempdir()
  expect_equal(file.size(paste0(.dir, "/rs62033413.pdf")), NA_real_) # ensure no plot file exists at baseline
  suppressWarnings(gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = paste0(.dir, "/")))
  expect_true(file.size(paste0(.dir, "/rs62033413.pdf")) > 0) # ensure gg_locusplot wrote a file
  unlink(paste0(.dir, "/rs62033413.pdf")) # remove temporary file
})
