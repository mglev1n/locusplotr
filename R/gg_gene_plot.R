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

  ggbio::autoplot(ensdb, AnnotationFilter::GRangesFilter(gr), names.expr="gene_name") +
    ggplot2::theme_minimal()

}
