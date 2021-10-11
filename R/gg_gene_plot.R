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
#' gg_gene_plot(1, 170054349 - 1e6, 170054349 + 1e6, "GRCh37")

gg_gene_plot <- function(chr, start, end, build) {

  if(build == "GRCh37") {
    txb <- AnnotationDbi::loadDb("inst/extdata/txb_hg19.sqlite")
    # txb <- txb_hg19
  } else if(build == "GRCh38") {
    txb <- AnnotationDbi::loadDb("inst/extdata/txb_hg38.sqlite")
    # txb <- txb_hg38
  } else {
    cli::cli_abort("Please provide a valid genome build")
  }

  # txb <- do.call(GenomicFeatures::makeTxDb, as.list(txb))

  gr <- GenomicRanges::GRanges(paste0("chr",chr), IRanges::IRanges(start, end), strand = "*")
  gr.txdb <- biovizBase::crunch(txb, which = gr)
  colnames(values(gr.txdb))[4] <- "model"
  grl <- split(gr.txdb, gr.txdb$gene_id)
  symbols <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=names(grl), columns="SYMBOL", keytype="ENTREZID")

  names(grl) <- symbols[match(symbols$ENTREZID, names(grl), nomatch=0), "SYMBOL"]

  # return(grl)

  ggplot2::update_geom_defaults("text", list(angle = 30, hjust = 0))

  plot_res <- ggbio::autoplot(grl, ggplot2::aes(type = model)) +
    ggplot2::theme_light(base_size = 16) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult=c(0.15,0.25)))

  plot_res@ggplot +
    ggplot2::theme_light(base_size = 16)
}
