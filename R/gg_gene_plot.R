#' gg_gene_plot
#'
#' Returns a ggplot containing the genes within a specified genomic region. The function uses database connections to EnsDb.Hsapiens.v75 (hg19/GRCh37) or EnsDb.Hsapiens.v86 (hg38/GRCh38) to identify genes within the specified region, and uses the ggbio package to create the plot.
#'
#' @param chr Integer - chromosome
#' @param start Integer - starting position for region of interest
#' @param end Integer - ending position for region of interest
#' @param genome_build Character - genome build - one of "GRCh37" or "GRCh38"
#'
#' @return A ggplot object containing a plot of genes within the region of interest
#' @export
#'
#' @examples
#' \dontrun{
#' gg_gene_plot(1, 170054349 - 1e6, 170054349 + 1e6, "GRCh37")
#' }
#'
gg_gene_plot <- function(chr, start, end, genome_build) {

  # ensure function arguments are of the correct type
  checkmate::assert_numeric(chr)
  checkmate::assert_numeric(start)
  checkmate::assert_numeric(end)
  checkmate::assert_choice(genome_build, choices = c("GRCh37", "GRCh38"))

  # select appropriate genome database for identifying genes within the specified region
  if (genome_build == "GRCh37") {
    txb <- AnnotationDbi::loadDb(system.file("extdata", "txb_hg19.sqlite", package = "locusplotr"))
    # txb <- txb_hg19
  } else if (genome_build == "GRCh38") {
    txb <- AnnotationDbi::loadDb(system.file("extdata", "txb_hg19.sqlite", package = "locusplotr"))
    # txb <- txb_hg38
  }

  # Identify genes in region, create model with exons and gene symbols as labels
  gr <- GenomicRanges::GRanges(paste0("chr", chr), IRanges::IRanges(start, end), strand = "*")
  gr.txdb <- biovizBase::crunch(txb, which = gr)
  colnames(values(gr.txdb))[4] <- "model"
  grl <- split(gr.txdb, gr.txdb$gene_id)
  suppressMessages(symbols <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = names(grl), columns = "SYMBOL", keytype = "ENTREZID"))

  names(grl) <- symbols[match(symbols$ENTREZID, names(grl), nomatch = 0), "SYMBOL"]

  grl <- GenomicRanges::GRangesList(grl, compress = TRUE)

  # override angle for plotting gene labels
  update_geom_defaults("text", list(angle = 30, hjust = 0))

  # plot genes
  suppressMessages(plot_res <- ggbio::autoplot(grl, aes(type = model)) +
    theme_light(base_size = 16) +
    scale_y_discrete(expand = expansion(mult = c(0.15, 0.25))))

  # return ggplot object with results
  suppressMessages(plot_res@ggplot +
    theme_light(base_size = 16))
}
