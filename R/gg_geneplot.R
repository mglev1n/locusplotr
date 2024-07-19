#' Plot genes located within a genomic region of interest
#'
#' Returns a ggplot containing the genes within a specified genomic region. The function uses database connections to EnsDb.Hsapiens.v75 (hg19/GRCh37) or EnsDb.Hsapiens.v86 (hg38/GRCh38) to identify genes within the specified region, and uses the ggbio package to create the plot.
#'
#' @param chr Integer - chromosome
#' @param start Integer - starting position for region of interest
#' @param end Integer - ending position for region of interest
#' @param genome_build Character - genome build - one of "GRCh37" or "GRCh38"
#' @param max_levels Integer - maximum number of levels for gene tracks
#'
#' @return A ggplot object containing a plot of genes within the region of interest
#' @export
#'
#' @examples
#' \dontrun{
#' gg_geneplot(1, 170054349 - 1e6, 170054349 + 1e6, "GRCh37")
#' }
#'

gg_geneplot <- function(chr, start, end, genome_build = "GRCh38", max_levels = 5) {
  checkmate::assert_numeric(chr)
  checkmate::assert_numeric(start)
  checkmate::assert_numeric(end)
  checkmate::assert_choice(genome_build, choices = c("GRCh37", "GRCh38"))
  # Select the appropriate gene table based on the genome version
  if (genome_build == "GRCh38") {
    gene_table <- snpsettest::gene.curated.GRCh38
  } else if (genome_build == "GRCh37") {
    gene_table <- snpsettest::gene.curated.GRCh37
  } else {
    stop("Invalid genome version. Use 'GRCh37' or 'GRCh38'.")
  }
  chromosome <- chr
  filter_start <- start
  filter_end <- end
  # Filter genes within the specified region
  genes <- gene_table %>%
    filter(chr == chromosome,
           start <= filter_end,
           end >= filter_start) %>%
    select(gene = gene.name, start, end, strand)

  # Trim genes that extend beyond the specified region
  genes <- genes %>%
    mutate(
      start = pmax(start, !!start),
      end = pmin(end, !!end)
    )

  # Check if any genes were found
  if (nrow(genes) == 0) {
    warning("No genes found in the specified region.")
    return(NULL)
  }

  # Assign y-levels to genes
  genes <- assign_y_levels(genes, max_levels, min_center_distance = 200000, min_end_distance = 20000)

  # Create the plot
  p <- ggplot(genes, aes(xmin = start, xmax = end, y = y_level)) +
    geom_segment(aes(x = start, xend = end, yend = y_level), linewidth = 2, color = "darkblue") +
    geom_text(aes(x = (start + end) / 2, label = gene), vjust = -0.5, size = 3) +
    scale_x_continuous(breaks = scales::extended_breaks(n = 5),
                       labels = scales::label_number(scale = 1 / 1e6),
                       limits = c(start, end)) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    labs(x = glue::glue("Position on Chromosome {chr} (Mb)"),
         y = "") +
    theme_bw(base_size = 16) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())

  return(p)
}

#' Assign y-levels to genes for even distribution in a plot
#'
#' This function takes a dataframe of genes and assigns y-levels to them,
#' ensuring even distribution and preventing overlap in the resulting plot.
#'
#' @param genes A dataframe containing gene information. Must have columns 'start' and 'end'.
#' @param max_levels Integer. The initial maximum number of y-levels to use. Default is 5.
#' @param min_center_distance Integer. The minimum distance between the centers of two genes on the same level, in base pairs. Default is 200000 (200kb).
#' @param min_end_distance Integer. The minimum distance between the end of one gene and the start of the next on the same level, in base pairs. Default is 20000 (20kb).
#'
#' @return A dataframe similar to the input, with an additional column 'y_level' indicating the assigned level for each gene.
#'
#' @details
#' The function sorts genes by their start position and then assigns them to levels.
#' It ensures that genes on the same level are sufficiently spaced apart, both in terms of
#' their center positions and their end-to-start distances. If a gene cannot be placed on
#' any existing level, a new level is created.
#'
#' @note
#' The function may create more levels than the initial `max_levels` if necessary to
#' accommodate all genes while maintaining the specified distances.
#'
#' @examples
#' genes_df <- data.frame(
#'   gene = c("Gene1", "Gene2", "Gene3"),
#'   start = c(1000, 5000, 10000),
#'   end = c(2000, 7000, 12000)
#' )
#' result <- assign_y_levels(genes_df, max_levels = 3, min_center_distance = 5000, min_end_distance = 1000)
#'
#' @seealso \code{\link{gg_geneplot}} for the main function that uses this to create gene plots.
#' @noRd
assign_y_levels <- function(genes, max_levels = 5, min_center_distance = 200000, min_end_distance = 20000) {
  genes <- genes[order(genes$start), ]
  n_genes <- nrow(genes)

  # Initialize levels
  levels <- vector("list", max_levels)

  for (i in 1:n_genes) {
    placed <- FALSE
    for (j in 1:max_levels) {
      if (length(levels[[j]]) == 0) {
        levels[[j]] <- c(levels[[j]], i)
        placed <- TRUE
        break
      } else {
        last_gene_index <- levels[[j]][length(levels[[j]])]
        last_gene_center <- (genes$start[last_gene_index] + genes$end[last_gene_index]) / 2
        current_gene_center <- (genes$start[i] + genes$end[i]) / 2
        center_distance <- current_gene_center - last_gene_center
        end_distance <- genes$start[i] - genes$end[last_gene_index]

        if (center_distance >= min_center_distance && end_distance >= min_end_distance) {
          levels[[j]] <- c(levels[[j]], i)
          placed <- TRUE
          break
        }
      }
    }
    if (!placed) {
      # If we couldn't place the gene, create a new level
      max_levels <- max_levels + 1
      levels[[max_levels]] <- i
    }
  }

  # Assign y-coordinates
  y_coords <- rep(0, n_genes)
  for (i in 1:length(levels)) {
    y_coords[levels[[i]]] <- i
  }

  genes$y_level <- y_coords
  return(genes)
}
