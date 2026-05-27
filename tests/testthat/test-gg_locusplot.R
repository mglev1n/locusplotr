test_that("Incorrect function argument type returns error", {
  expect_error(gg_locusplot(df = "fto_locus_df", lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_distance = 1e6, path = NULL))
  # expect_error(gg_locusplot(df = fto_locus_df, lead_snp = 1, rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL))
  expect_error(gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_title = NULL, plot_subtitle = NULL, plot_distance = "X", path = NULL))
  expect_error(gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL, plot_genes = "X"))
})

test_that("Default gg_locusplot returns a ggplot object", {
  locusplot_res <- gg_locusplot(df = fto_locus_df %>% mutate(trait = "BMI"), lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_recombination = TRUE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL, trait = trait)
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


# Generating sim'd LD data ( No checks to the validity of the data itself (i.e, Same SNP not having 1 or NA with itself))
mock_ld <- data.frame(
  rsid1 = rep("rs62033413", 2),
  rsid2 = c("rs62033413", "rs142090714"),
  r = c(1.0, 0.8)
)

mock_ld_caps <- data.frame(
  RSID1 = rep("rs62033413", 2),
  RSID2 = c("rs62033413", "rs142090714"),
  R = c(1.0, 0.8)
)

mock_ld_invalid <- data.frame(
  snp1 = rep("rs62033413", 2),
  snp2 = c("rs62033413", "rs142090714"),
  correlation = c(NA, 0.8)
)

mock_ld_negative <- data.frame(
  rsid1 = rep("rs62033413", 2),
  rsid2 = c("rs62033413", "rs142090714"),
  r = c(1.0, -0.8) # Negative r value to test the abs() conversion
)

mock_ld_scrambled_headers <- data.frame(
  correlation_value = c(1.0, 0.8),
  second_snp = c("rs62033413", "rs142090714"),
  first_snp = rep("rs62033413", 2)
)

test_that("Custom ld_df returns a ggplot object", {
  locusplot_res <- gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", ld_df = mock_ld, rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res, "ggplot")
})

test_that("Custom ld_df with case-insensitive columns returns a ggplot object", {
  locusplot_res <- gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", ld_df = mock_ld_caps, rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res, "ggplot")
})

test_that("Custom ld_df with missing columns returns error", {
  expect_error(
    gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", ld_df = mock_ld_invalid, rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL),
    regexp = "ld_df must contain 'RSID1', 'RSID2', and 'r' columns"
  )
})


test_that("Genome build variations translate correctly without error", {
  locusplot_res_hg19 <- gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", genome_build = "hg19", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res_hg19, "ggplot")
  locusplot_res_hg38 <- gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", genome_build = "HG38", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res_hg38, "ggplot")
})

test_that("Invalid genome build returns error", {
  expect_error(
    gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", genome_build = "hg39", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL),
    regexp = "Assertion on 'toupper\\(genome_build\\)' failed"
  )
})

test_that("Custom ld_df handles negative correlation values", {
  locusplot_res <- gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", ld_df = mock_ld_negative, rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
  expect_s3_class(locusplot_res, "ggplot")
})

test_that("Custom ld_df throws error if not a data frame", {
  expect_error(
    gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", ld_df = "not_a_dataframe", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL),
    regexp = "Assertion on 'ld_df' failed"
  )
})

test_that("Unsupported GRC genome builds return an error", {
  expect_error(
    gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", genome_build = "GRCh36", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL),
    regexp = "Assertion on 'toupper\\(genome_build\\)' failed"
  )
})

test_that("Custom ld_df with completely wrong and out-of-order headers returns error", {
  expect_error(
    gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", ld_df = mock_ld_scrambled_headers, rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL),
    regexp = "ld_df must contain 'RSID1', 'RSID2', and 'r' columns"
  )
})

test_that("Custom ld_df alerts the user about genome build alignment when background tracks are active", {
  expect_message(
    gg_locusplot(
      df = fto_locus_df, lead_snp = "rs62033413",ld_df = mock_ld,  genome_build = "GRCh37",
      plot_recombination = TRUE, # Triggers the API track needing the build
      rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL
    ),
    regexp = "Ensure that your supplied 'genome_build'"
  )
})
