# Function to plot regional association with LD
#' Create a regional association plot
#'
#' Returns a ggplot object containing a regional association plot (-log10(p-value) as a function of chromosomal position, with variants colored by linkage disequilibrium to reference variant).
#' This function allows the user to integrate genome wide association study (GWAS) summary statistics for a locus of interest with linkage disequilibrium information (obtained using the University of Michigan LocusZoom API <https://portaldev.sph.umich.edu/>) for that locus to create a regional association plot.
#'
#' @param df Dataframe containing columns with rsid, chromosome, position, reference/effect allele, alternate/non-effect allele, and p-value for all variants within the range of interest
#' @param lead_snp A character vector containing a lead variant of interest. When NULL (default), the variant with the lowest p-value will be selected as the lead variant.
#' @param rsid Rsid column
#' @param chrom Chromosome column
#' @param pos Position column
#' @param ref Reference/effect allele column
#' @param alt Alternate/non-effect allele column
#' @param p_value P-value column
#' @param plot_pvalue_threshold Threshold for plotting p-value on regional association plot (default = 0.1) - reducing the number of points decreases file size and improves performance
#' @param plot_subsample_prop Proportion of points above p-value threshold to plot (default = 0.1; range = 0-1) - reducing the number of points decreases file size and improves performance
#' @param plot_distance Integer corresponding to the size of the locus that should be plotted
#' @param genome_build Character - one of "GRCh37" or "GRCh38"
#' @param population Character - one of "ALL", "AFR", "AMR", "EAS", "EUR", "SAS" referring to the reference population of interest for obtaining linkage disequilibrium information (default = "ALL")
#' @param plot_genes Logical - Include a plot of genes/transcripts within the region of interest beneath the regional association plot (default = FALSE)
#' @param plot_title A character string corresponding to plot title (default = NULL)
#' @param plot_subtitle A character string corresponding to plot subtitle (default = NULL)
#' @param path Character string (default = NULL) - if a path is supplied a .pdf of the plot will be saved
#' @param trait (optional) Column containing the name of the trait
#'
#' @return A ggplot object containing a regional association plot for the locus of interest
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic regional association plot
#' gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value)
#'
#' # Use "plot_genes = TRUE" to add a plot of genes within the region beneath the regional association plot
#' gg_locusplot(df = fto_locus_df, lead_snp = "rs62033413", rsid = rsid, chrom = chromosome, pos = position, ref = effect_allele, alt = other_allele, p_value = p_value, plot_genes = TRUE)
#' }
#'
gg_locusplot <- function(df, lead_snp = NULL, rsid = rsid, chrom = chrom, pos = pos, ref = ref, alt = alt, p_value = p_value, trait = NULL, plot_pvalue_threshold = 0.1, plot_subsample_prop = 0.1, plot_distance = 500000, genome_build = "GRCh37", population = "ALL", plot_genes = FALSE, plot_title = NULL, plot_subtitle = NULL, path = NULL) {
  # Check input arguments to ensure they are of the correct type and within reasonable ranges
  checkmate::assert_data_frame(df)
  # checkmate::assert_string(lead_snp)
  checkmate::assert_numeric(plot_pvalue_threshold, upper = 1)
  checkmate::assert_numeric(plot_subsample_prop, lower = 0, upper = 1)
  checkmate::assert_numeric(plot_distance, lower = 0)
  checkmate::assert_logical(plot_genes)

  # trait <- rlang::enquo(trait)

  if (rlang::quo_is_null(rlang::enquo(trait))) {
    df <- df %>%
      select(rsid = {{ rsid }}, chromosome = {{ chrom }}, position = {{ pos }}, ref = {{ ref }}, alt = {{ alt }}, p_value = {{ p_value }}) %>%
      mutate_if(is.factor, as.character) %>%
      mutate(ref = stringr::str_to_upper(ref), alt = stringr::str_to_upper(alt)) %>%
      group_by(rsid) %>%
      slice_min(p_value) %>%
      ungroup() %>%
      tidyr::drop_na()
  } else {
    df <- df %>%
      select(rsid = {{ rsid }}, chromosome = {{ chrom }}, position = {{ pos }}, ref = {{ ref }}, alt = {{ alt }}, p_value = {{ p_value }}, trait = {{ trait }}) %>%
      mutate_if(is.factor, as.character) %>%
      mutate(ref = stringr::str_to_upper(ref), alt = stringr::str_to_upper(alt)) %>%
      arrange(p_value) %>%
      group_by(trait, rsid) %>%
      slice_min(p_value) %>%
      ungroup() %>%
      tidyr::drop_na()
  }


  # Create df containing information about lead SNP (by default, select SNP with lowest p-value, otherwise take user-supplied value)
  if (is.null(lead_snp)) {
    indep_snps <- df %>%
      slice_min(p_value, with_ties = FALSE, n = 1) %>%
      select(lead_rsid = rsid, lead_chromosome = chromosome, lead_position = position, lead_ref = ref, lead_alt = alt)

    cli::cli_alert_info("No lead_snp supplied. Defaulting to {indep_snps$lead_rsid} - {indep_snps$lead_chromosome}:{indep_snps$lead_position}, which has the lowest p-value in the region")
  } else if (!(lead_snp %in% df$rsid)) {
    # ensure lead_snp is in the supplied data; if not, use minimum p-value at locus
    indep_snps <- df %>%
      slice_min(p_value, with_ties = FALSE, n = 1) %>%
      select(lead_rsid = rsid, lead_chromosome = chromosome, lead_position = position, lead_ref = ref, lead_alt = alt)

    cli::cli_alert_info("Lead snp not present in supplied locus data. Defaulting to {indep_snps$lead_rsid} - {indep_snps$lead_chromosome}:{indep_snps$lead_position}, which has the lowest p-value in the region")
  } else {
    # Ensure lead_snp supplied by user is a string
    # checkmate::assert_string(lead_snp)

    indep_snps <- df %>%
      select(lead_rsid = rsid, lead_chromosome = chromosome, lead_position = position, lead_ref = ref, lead_alt = alt) %>%
      filter(lead_rsid == lead_snp) %>%
      distinct(lead_rsid, .keep_all = TRUE)
  }

  # Create dataframe of variants within the region size specified by the user
  suppressMessages(locus_snps <- df %>%
                     filter(rsid %in% indep_snps$lead_rsid) %>%
                     select(chromosome, position, lead_rsid = rsid) %>%
                     purrr::pmap_dfr(function(chromosome_filter = first, position_filter = second, lead_rsid = third) {
                       df %>%
                         filter(chromosome == chromosome_filter & between(position, position_filter - plot_distance / 2, position_filter + plot_distance / 2)) %>%
                         mutate(lead_rsid = lead_rsid) %>%
                         left_join(indep_snps)
                     }))

  # Extract LD and format colors
  possibly_ld_extract_locuszoom <- purrr::possibly(locusplotr::ld_extract_locuszoom, otherwise = NULL)

  ld_extracted <- possibly_ld_extract_locuszoom(chrom = indep_snps$lead_chromosome, pos = indep_snps$lead_position, ref = indep_snps$lead_ref, alt = indep_snps$lead_alt, start = min(locus_snps$position), stop = max(locus_snps$position), genome_build = genome_build, population = population)

  # Create dataframe with variants at locus, LD information, color codes, and labels in preparation for plotting
  if (!(is.null(ld_extracted))) {
    # Join GWAS locus df with LD information
    locus_snps_ld <- ld_extracted %>%
      select(chromosome = chromosome2, position = position2, variant2, correlation) %>%
      mutate(chromosome = as.numeric(chromosome), position = as.numeric(position)) %>%
      tidyr::separate(variant2, into = c("chr_pos", "ref_alt"), sep = "_") %>%
      tidyr::separate(ref_alt, into = c("ref", "alt"), sep = "/") %>%
      right_join(locus_snps, by = c("chromosome" = "chromosome", "position" = "position"), relationship = "many-to-many") %>%
      filter((ref.x == ref.y & alt.x == alt.y) | (ref.x == alt.y & alt.x == ref.y)) %>%
      select(-ends_with(".y"), -chr_pos) %>%
      rename_with(~ stringr::str_replace(.x, ".x", ""), .cols = ends_with(".x"))

    # Create color codes and labels
    locus_snps_ld <- locus_snps_ld %>%
      mutate(color_code = as.character(cut(as.numeric(correlation), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("blue4", "skyblue", "darkgreen", "orange", "red"), include.lowest = TRUE))) %>%
      mutate(legend_label = as.character(cut(as.numeric(correlation), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1"), include.lowest = TRUE))) %>%
      mutate(lead = rsid == lead_rsid) %>%
      mutate(label = case_when(
        rsid == lead_rsid ~ lead_rsid,
        TRUE ~ NA_character_
      )) %>%
      mutate(color_code = case_when(
        rsid == lead_rsid ~ "purple",
        TRUE ~ color_code
      )) %>%
      mutate(color_code = forcats::fct_relevel(color_code, "purple", "red", "orange", "darkgreen", "skyblue", "blue4")) %>%
      mutate(legend_label = case_when(
        rsid == lead_rsid ~ "Ref",
        TRUE ~ legend_label
      )) %>%
      mutate(legend_label = forcats::fct_relevel(legend_label, "Ref", "0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2"))
  } else {
    # Deal with scenario where lead variant is not present in LD database
    cli::cli_alert_info("No linkage disequilibrium information found")
    locus_snps_ld <- locus_snps %>%
      mutate(correlation = NA_integer_) %>%
      mutate(lead = rsid == lead_rsid) %>%
      mutate(label = case_when(
        rsid == lead_rsid ~ lead_rsid,
        TRUE ~ NA_character_
      )) %>%
      mutate(color_code = case_when(
        rsid == lead_rsid ~ "purple",
        TRUE ~ "grey50"
      )) %>%
      mutate(legend_label = case_when(
        rsid == lead_rsid ~ "Ref",
        TRUE ~ "Other"
      ))
  }

  # group locus by trait if necessary
  if (!rlang::quo_is_null(rlang::enquo(trait))) {
    locus_snps_ld <- locus_snps_ld %>%
      group_by(.data = ., trait)
    }

  # Make plot (sample non-significant p-values to reduce overplotting)
  suppressMessages(regional_assoc_plot <- locus_snps_ld %>%
                     distinct(rsid, .keep_all = TRUE) %>%
                     filter(p_value < plot_pvalue_threshold | correlation > 0.2 | legend_label == "Ref") %>% # improve overplotting
                     bind_rows(locus_snps_ld %>%
                                 filter(p_value >= plot_pvalue_threshold & correlation < 0.2 & legend_label != "Ref") %>%
                                 slice_sample(prop = plot_subsample_prop)) %>%
                     arrange(desc(color_code)) %>%
                     ggplot(aes(position, -log10(p_value), fill = factor(color_code), size = lead, alpha = lead, shape = lead)) +
                     geom_point() +
                     ggrepel::geom_label_repel(aes(label = label),
                                               size = 4,
                                               color = "black",
                                               fontface = "bold",
                                               fill = "white",
                                               min.segment.length = 0,
                                               box.padding = 1,
                                               alpha = 1,
                                               # nudge_x = -0.05 * max(locus_snps_ld$position),
                                               # nudge_y = 0.25 * -log10(min(locus_snps_ld$p_value))
                                               nudge_y = 5
                     ) +
                     geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
                     scale_fill_identity(parse(text = "r^2"), guide = "legend", labels = levels(forcats::fct_drop(locus_snps_ld$legend_label)), na.translate = FALSE) +
                     scale_size_manual(values = c(3, 5), guide = "none") +
                     scale_shape_manual(values = c(21, 23), guide = "none") +
                     scale_alpha_manual(values = c(0.8, 1), guide = "none") +
                     scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1 / 1e6)) +
                     scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
                     guides(fill = guide_legend(override.aes = list(shape = 22, size = 6))) +
                     labs(
                       title = plot_title,
                       subtitle = plot_subtitle,
                       x = glue::glue("Position on Chromosome {unique(indep_snps$lead_chromosome)} (Mb)"),
                       y = "-log<sub>10</sub>(P-value)"
                     ) +
                     theme_bw(base_size = 16) +
                     theme(
                       plot.title = element_text(face = "bold"),
                       legend.title.align = 0.5,
                       legend.text = element_text(size = 10),
                       legend.title = element_text(size = 10),
                       # legend.margin = margin(1, 1, 1, 1),
                       legend.justification = c("right", "top"),
                       legend.position = c(0.99, 0.99),
                       # legend.spacing = unit(0, "pt"),
                       strip.text = element_text(color = "black"),
                       strip.text.x = element_blank(),
                       axis.title.y = ggtext::element_markdown(),
                       legend.spacing.y = unit(0, "pt")
                     ))

  if (!rlang::quo_is_null(enquo(trait))) {
    regional_assoc_plot <- regional_assoc_plot +
      facet_grid(rows = vars(trait))
  }

  # Add plot of genes if reuested by user
  if (plot_genes) {
    cli::cli_alert_info("Extracting genes for the region {indep_snps$lead_chromosome}:{indep_snps$lead_position - plot_distance/2}-{indep_snps$lead_position + plot_distance/2}")
    geneplot <- callr::r(function(chr, start, end, genome_build) {
      locusplotr::gg_geneplot(chr, start, end, genome_build) # nocov
    }, args = list(chr = indep_snps$lead_chromosome, start = indep_snps$lead_position - plot_distance / 2, end = indep_snps$lead_position + plot_distance / 2, genome_build = genome_build)) +
      labs(x = glue::glue("Position on Chromosome {indep_snps$lead_chromosome} (Mb)")) +
      # scale_fill_brewer(palette = "Set3", guide = "none") +
      scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1 / 1e6), limits = c(indep_snps$lead_position - plot_distance / 2, indep_snps$lead_position + plot_distance / 2)) +
      theme(plot.margin = margin(0, 5.5, 5.5, 5.5))

    suppressWarnings(suppressMessages(regional_assoc_plot <- patchwork::wrap_plots(list(
      regional_assoc_plot +
        labs(x = "") +
        xlim(indep_snps$lead_position - plot_distance / 2, indep_snps$lead_position + plot_distance / 2) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(5.5, 5.5, 0, 5.5)
        ),
      geneplot
    ), nrow = 2, heights = c(2, 1))))
  }

  # Return +/- save ggplot object
  if (!is.null(path)) {
    ggsave(regional_assoc_plot, filename = paste0(path, stringr::str_replace_all(unique(indep_snps$lead_rsid), "[^[:alnum:]]", "_"), ".pdf"), units = "in", height = 8.5, width = 11, device = "pdf")
  }
  # } else {
  #   ggsave(regional_assoc_plot, filename = paste0(path, stringr::str_replace_all(unique(indep_snps$lead_rsid), "[^[:alnum:]]", "_"), ".pdf"), units = "in", height = 8.5, width = 11, device = "pdf")
  #   return(regional_assoc_plot)
  # }

  return(regional_assoc_plot)
}
