#' gg_gene_plot
#'
#' Returns a ggplot containing the genes within a specified genomic region. The function uses database connections to EnsDb.Hsapiens.v75 (hg19/GRCh37) or EnsDb.Hsapiens.v86 (hg38/GRCh38) to identify genes within the specified region, and uses the ggbio package to create the plot.
#'
#' @param chr Integer - chromosome
#' @param start Integer - starting position for region of interest
#' @param end Integer - ending position for region of interest
#' @param build Character - one of "GRCh37" or "GRCh38"
#'
#' @return
#' @export
#'
#' @examples
gg_gene_plot <- function(chr, start, end, build) {

  if(build == "GRCh37") {
    ensdb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
  } else if(build == "GRCh38") {
    ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  } else {
    cli::cli_abort("Please provide a valid genome build")
  }

  gr <- GenomicRanges::GRanges(seqnames = chr, IRanges::IRanges(start, end), strand = "*")

  ensembldb::genes(ensdb) %>%
    plyranges::join_overlap_inner(gr) %>%
    plyranges::filter(gene_biotype == "protein_coding") %>%
    as_tibble() %>%
    mutate(direction = case_when(
      strand == "*" ~ "+",
      TRUE ~ as.character(strand)
    )) %>%
    dplyr::select(gene = gene_name, seq_id, start, end, strand) %>%
    gggenomes::gggenomes(infer_start = 0) +
    # gggenomes::geom_seq() +
    gggenomes::geom_gene(position = gggenomes::position_pile(offset = 0.2), shape = 4, fill = "gray80") +
    gggenomes::geom_gene_tag(aes(label = gene), size = 3, position = gggenomes::position_pile(offset = 0.2), angle = 15) +
    scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1/1e6)) +
    labs(x = glue::glue("Position on Chromosome {chr} (Mb)")) +
    ggplot2::theme_light(base_size = 16) +
    ggplot2::theme(panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.text.y = element_blank())

  # ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75) %>%
  #   plyranges::join_overlap_inner(GenomicRanges::GRanges(1, ranges = IRanges::IRanges(169519049 - 500000, 169519049 + 500000),
  #                                                        strand = "*")) %>%
  #   plyranges::filter(gene_biotype == "protein_coding") %>%
  #   as_tibble() %>%
  #   mutate(direction = case_when(
  #     strand == "*" ~ "+",
  #     TRUE ~ as.character(strand)
  #   )) %>%
  #   dplyr::select(gene = gene_name, seq_id, start, end, strand) %>%
  #   gggenomes::gggenomes(infer_start = 0) +
  #   # gggenomes::geom_seq() +
  #   gggenomes::geom_gene(position = gggenomes::position_pile(offset = 0.2), shape = 4, fill = "gray80") +
  #   gggenomes::geom_gene_tag(aes(label = gene), size = 3, position = gggenomes::position_pile(offset = 0.2), angle = 15) +
  #   scale_x_continuous(breaks = scales::extended_breaks(n = 5), labels = scales::label_number(scale = 1/1e6)) +
  #   ggplot2::theme_light(base_size = 16) +
  #   ggplot2::theme(panel.grid.major.y = element_blank(),
  #                  panel.grid.minor.y = element_blank(),
  #                  axis.title.y = element_blank(),
  #                  axis.ticks.y = element_blank(),
  #                  axis.text.y = element_blank())

}



