# Function to plot regional association with LD
#' gg_locusplot
#'
#' Returns a ggplot object containing a regional association plot (-log10(p-value) as a function of chromosomal position, with variants colored by linkage disequilibrium to reference variant).
#' This function allows the user to integrate genome wide association study (GWAS) summary statistics for a locus of interest with linkage disequilibrium information (obtained using the University of Michigan LocusZoom API <https://portaldev.sph.umich.edu/>) for that locus to create a regional association plot.
#'
#' @param df Dataframe containing columns with rsid, chromosome, position, reference/effect allele, alternate/non-effect allele, and p-value for all variants within the range of interest
#' @param lead_snps A character vector containing the lead variant of interest
#' @param rsid Rsid column
#' @param chromosome Chromosome column
#' @param position Position column
#' @param ref Reference/effect allele column
#' @param alt Alternate/non-effect allele column
#' @param p_value P-value column
#' @param plot_pvalue_threshold Threshold for plotting p-value on regional association plot (default = 0.1) - reducing the number of points decreases file size and improves performance
#' @param plot_distance Integer corresponding to the size of the locus that should be plotted
#' @param genome_build Character - one of "GRCh37" or "GRCh38"
#' @param population Character - one of "ALL", "AFR", "AMR", "EAS", "EUR", "SAS" referring to the reference population of interest for obtaining linkage disequilibrium information
#' @param plot_genes Logical - Include a plot of genes/transcripts within the region of interest beneath the regional association plot (default = FALSE)
#' @param plot_title A character string corresponding to plot title (default = NULL)
#' @param plot_subtitle A character string corresponding to plot subtitle (default = NULL)
#' @param path Character string (default = NULL) - if a path is supplied a .pdf of the plot will be saved
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' library(locusplotr)
#' library(tidyverse)
#' test <- arrow::read_parquet("../Varicose-Veins-GWAMA/Data/HyPrColoc/vv_risk_loci_transancestry_20210324_FUMA_20210928.parquet") %>%
#'   filter(lead_rsid == "rs6025")

#' test %>%
#'   gg_regional_assoc(lead_snps = "rs6025", , rsid = rsid, chromosome = CHR, position = POSITION, ref = Allele1, alt = Allele2, p_value = `P-value`, plot_genes = TRUE, plot_title = NULL, plot_subtitle = NULL, plot_distance = 1e6, path = NULL)
#' }


gg_locusplot <- function(df, lead_snps, rsid = rsid, chromosome = chromosome, position = position, ref = ref, alt = alt, p_value = p_value, plot_pvalue_threshold = 0.1, plot_distance = 500000, genome_build = "GRCh37", population = "ALL", plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, path = NULL) {
  df <- df %>%
    select(rsid = {{ rsid }}, chromosome = {{ chromosome }}, position = {{ position }}, ref = {{ ref }}, alt = {{ alt }}, p_value = {{ p_value }}) %>%
    mutate_if(is.factor, as.character) %>%
    mutate(ref = stringr::str_to_upper(ref), alt = stringr::str_to_upper(alt))

  indep_snps <- df %>%
    select(lead_rsid = rsid, lead_chromosome = chromosome, lead_position = position, lead_ref = ref, lead_alt = alt) %>%
    filter(lead_rsid %in% lead_snps)

  # return(indep_snps)

  suppressMessages(locus_snps <- df %>%
    filter(rsid %in% indep_snps$lead_rsid) %>%
    select(chromosome, position, lead_rsid = rsid) %>%
    pmap_dfr(function(chromosome_filter = first, position_filter = second, lead_rsid = third) {
      df %>%
        filter(chromosome == chromosome_filter & between(position, position_filter - plot_distance / 2, position_filter + plot_distance / 2)) %>%
        dplyr::mutate(lead_rsid = lead_rsid) %>%
        left_join(indep_snps)
    }))

  # return(locus_snps)

  # Extract LD and format colors
  # consider adding error handling if we can't extract LD - eg. just plot without colors
  possibly_ld_extract_locuszoom <- purrr::possibly(ld_extract_locuszoom, otherwise = NULL)

  ld_extracted <- possibly_ld_extract_locuszoom(chrom = indep_snps$lead_chromosome, pos = indep_snps$lead_position, ref = indep_snps$lead_ref, alt = indep_snps$lead_alt, start = min(locus_snps$position), stop = max(locus_snps$position), build = genome_build, population = population)

  if ((dim(ld_extracted)[1] != 0)) {
    locus_snps_ld <- ld_extracted %>%
      select(chromosome = chromosome2, position = position2, variant2, correlation) %>%
      mutate(chromosome = as.numeric(chromosome), position = as.numeric(position)) %>%
      separate(variant2, into = c("chr_pos", "ref_alt"), sep = "_") %>%
      separate(ref_alt, into = c("ref", "alt"), sep = "/") %>%
      right_join(locus_snps, by = c("chromosome" = "chromosome", "position" = "position")) %>%
      filter((ref.x == ref.y & alt.x == alt.y) | (ref.x == alt.y & alt.x == ref.y)) %>%
      select(-ends_with(".y"), -chr_pos) %>%
      rename_with(~ str_replace(.x, ".x", ""), .cols = ends_with(".x")) %>%
      mutate(color_code = as.character(cut(correlation, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("blue4", "skyblue", "darkgreen", "orange", "red"), include.lowest = TRUE))) %>%
      mutate(legend_label = as.character(cut(correlation, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1"), include.lowest = TRUE))) %>%
      mutate(lead = rsid == lead_rsid) %>%
      mutate(label = case_when(
        rsid == lead_rsid ~ lead_rsid,
        TRUE ~ NA_character_
      )) %>%
      mutate(color_code = case_when(
        rsid == lead_rsid ~ "purple",
        TRUE ~ color_code
      )) %>%
      mutate(color_code = fct_relevel(color_code, "purple", "red", "orange", "darkgreen", "skyblue", "blue4")) %>%
      mutate(legend_label = case_when(
        rsid == lead_rsid ~ "Reference",
        TRUE ~ legend_label
      )) %>%
      mutate(legend_label = fct_relevel(legend_label, "Reference", "0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2"))
  } else {
    locus_snps_ld <- locus_snps %>%
      mutate(correlation = NA_integer_) %>%
      mutate(lead = rsid == lead_rsid) %>%
      mutate(label = case_when(
        rsid == lead_rsid ~ lead_rsid,
        TRUE ~ NA_character_
      )) %>%
      mutate(color_code = case_when(
        rsid == lead_rsid ~ "purple",
        TRUE ~ color_code
      )) %>%
      mutate(color_code = fct_relevel(color_code, "purple", NA_character_)) %>%
      mutate(legend_label = case_when(
        rsid == lead_rsid ~ "Reference",
        TRUE ~ legend_label
      )) %>%
      mutate(legend_label = fct_relevel(legend_label, "Reference", NA_character_))
  }


  # return(locus_snps_ld)

  # Make Plot
  suppressMessages(regional_assoc_plot <- locus_snps_ld %>%
    filter(p_value < plot_pvalue_threshold) %>%
    # improve overplotting
    arrange(desc(color_code)) %>%
    ggplot(aes(position, -log10(p_value), fill = factor(color_code), size = lead, alpha = lead, shape = lead)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = label), size = 4, color = "black", fontface = "bold", min.segment.length = 0, box.padding = 1, alpha = 1, nudge_x = -0.1 * max(locus_snps_ld$position)) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
    scale_fill_identity(parse(text = "r^2"), guide = "legend", labels = levels(locus_snps_ld$legend_label), na.translate = FALSE) +
    scale_size_manual(values = c(3, 8), guide = "none") +
    scale_shape_manual(values = c(21, 23), guide = "none") +
    scale_alpha_manual(values = c(0.8, 1), guide = "none") +
    scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1/1e6)) +
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 6))) +
    # facet_grid(trait ~ lead_rsid, scales = "free") +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = glue::glue("Position on Chromosome {unique(indep_snps$lead_chromosome)} (Mb)"),
      y = bquote(-log[10]("P-value"))
    ) +
    theme_light(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title.align = 0.5,
      legend.key = element_rect(size = 3, fill = NA, colour = NA),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      # legend.position = "top",
      legend.margin = margin(1,1,1,1),
      legend.justification = c("right", "top"),
      legend.position = c(0.99, 0.99),
      legend.spacing = unit(0, "pt"),
      # legend.box.background = element_rect(color = "grey30"),
      # strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(color = "black", face = "bold"),
      strip.text.x = element_blank()
    ))

  if(plot_genes) {
    cli::cli_alert_info("Extracting genes for the region {indep_snps$lead_chromosome}:{indep_snps$lead_position - plot_distance/2}-{indep_snps$lead_position + plot_distance/2}")
    gene_plot <- callr::r(function(chr, start, end, build) {locusplotr::gg_gene_plot(chr, start, end, build)}, args = list(chr = indep_snps$lead_chromosome, start = indep_snps$lead_position - plot_distance/2, end = indep_snps$lead_position + plot_distance/2, build = genome_build)) +
      labs(x = glue::glue("Position on Chromosome {indep_snps$lead_chromosome} (Mb)")) +
      # scale_fill_brewer(palette = "Set3", guide = "none") +
      scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1/1e6), limits = c(indep_snps$lead_position - plot_distance/2, indep_snps$lead_position + plot_distance/2)) +
      theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

    suppressWarnings(suppressMessages(regional_assoc_plot <- patchwork::wrap_plots(list(regional_assoc_plot +
                                                        labs(x = "") +
                                                        xlim(indep_snps$lead_position - plot_distance/2, indep_snps$lead_position + plot_distance/2) +
                                                        ggplot2::theme(axis.text.x = element_blank(),
                                                                       axis.ticks.x = element_blank(),
                                                                       axis.title.x = element_blank(),
                                                                       plot.margin = margin(5.5, 5.5, 0, 5.5)),
                                                      gene_plot), nrow = 2, heights = c(2, 1))))
  }

  if (is.null(path)) {
    return(suppressWarnings(suppressMessages(regional_assoc_plot)))
  } else {
    suppressWarnings(suppressMessages(ggsave(regional_assoc_plot, filename = paste0(path, unique(indep_snps$lead_rsid), ".pdf"), units = "in", height = 8.5, width = 11, device = "pdf")))
    return(suppressWarningssuppressMessages(regional_assoc_plot))
  }
}
