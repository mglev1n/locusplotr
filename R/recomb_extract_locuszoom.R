#' Query LocusZoom API for Recombination Data
#'
#' This function queries the LocusZoom API to retrieve recombination data
#' for a specified genomic region and returns the result as a tibble.
#'
#' @param chrom A numeric value specifying the chromosome (e.g., 1, 2, ..., 22, 23 for X, 24 for Y)
#' @param start An integer specifying the start position of the region of interest
#' @param end An integer specifying the end position of the region of interest
#' @param genome_build A character string specifying the genome build (default: "GRCh37")
#'
#' @return A tibble containing the parsed recombination data from the LocusZoom API
#'
#' @examples
#' \dontrun{
#' result <- recomb_locuszoom(chrom = 1, start = 1000, end = 150000)
#' print(result)
#'
#' # Using a different genome build
#' result_grch38 <- recomb_locuszoom(chrom = 1, start = 1000, end = 150000, genome_build = "GRCh38")
#' }
#'
#' @export
recomb_extract_locuszoom <- function(chrom, start, end, genome_build = "GRCh37") {
  checkmate::assert_numeric(chrom)
  checkmate::assert_numeric(start)
  checkmate::assert_numeric(end)
  checkmate::assert_choice(genome_build, choices = c("GRCh37", "GRCh38"))

  # Construct the filter string
  filter <- glue::glue("chromosome eq '{chrom}' and position ge {start} and position le {end}")

  # Make the API request
  res <- httr::GET(
    url = "https://portaldev.sph.umich.edu/api/v1/annotation/recomb/results/",
    query = list(
      filter = filter,
      build = genome_build
    )
  )

  # Check if the request was successful
  if (httr::status_code(res) == 200) {
    # Parse the content
    content <- httr::content(res, "parsed")

    # Extract the data and convert to a tibble
    data_tibble <- tibble::tibble(
      chromosome = purrr::map_chr(content$data$chromosome, ~.x[[1]]),
      id = purrr::map_chr(content$data$id, ~as.character(.x[[1]])),
      pos_cm = purrr::map_dbl(content$data$pos_cm, ~.x[[1]]),
      position = purrr::map_int(content$data$position, ~.x[[1]]),
      recomb_rate = purrr::map_dbl(content$data$recomb_rate, ~.x[[1]])
    )

    return(data_tibble)
  } else {
    # If the request failed, return an error message
    stop(glue::glue("API request failed with status code {httr::status_code(res)}"))
  }
}

# Example usage:
# result <- recomb_locuszoom(chrom = 1, start = 1000, end = 150000)
# print(result)
