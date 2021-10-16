# Function to extract LD
#' ld_extract_locuszoom
#'
#' This function allows the user to extract linkage disequilibrium statistics between a variant of interest and surrounding variants within a contiguous genomic region. This function uses the University of Michigan LocusZoom API (<https://portaldev.sph.umich.edu/>) to obtain LD information, and allows the user to specify genome-build and ancestry of interest.
#'
#' @param chrom Integer - chromosome of reference variant
#' @param pos Integer - position of reference variant
#' @param ref Character - reference allele (or effect allele) for reference variant
#' @param alt Character - alternate allele (or non-effect allele) for reference variant
#' @param start Integer - starting position of range of interest
#' @param stop Integer - ending position of range of interest
#' @param genome_build Character - Genome build - one of "GRCh37" or "GRCh38"
#' @param population Character - 1000G genetic superpopulation - one of "ALL", "AFR", "AMR", "EAS", "EUR", "SAS"
#' @param metric Character - one of "r", "rsquare", or "cov", referring to the correlation statistic of interest
#'
#' @return A tibble containing each variant within the supplied range surrounding the variant of interest, with the requested linkage disequilibrium information with respect to the variant of interest
#' @export
#'
#' @examples
#' \dontrun{
#' ld_extract_locuszoom(chrom = 16, pos = 53830055, ref = "C", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5, build = "GRCh37", population = "ALL", metric = "rsquare")
#' }
#'
ld_extract_locuszoom <- function(chrom, pos, ref, alt, start, stop, genome_build = "GRCh37", population = "ALL", metric = "rsquare") {

  # Check function arguments
  checkmate::assert_numeric(chrom)
  checkmate::assert_numeric(pos)
  checkmate::assert_numeric(start)
  checkmate::assert_numeric(stop)
  checkmate::assert_true(stringr::str_detect(ref, stringr::regex("A|C|G|T", ignore_case = TRUE)), .var.name = "ref must contain A/C/G/T")
  checkmate::assert_true(stringr::str_detect(alt, stringr::regex("A|C|G|T", ignore_case = TRUE)), .var.name = "alt must contain A/C/G/T")
  checkmate::assert_choice(genome_build, choices = c("GRCh37", "GRCh38"))
  checkmate::assert_choice(population, choices = c("ALL", "AMR", "AFR", "EAS", "EUR", "SAS"))
  checkmate::assert_choice(metric, choices = c("r", "rsquare", "cov"))

  # Message
  cli::cli_alert_info("Extracting LD for {chrom}:{pos}_{stringr::str_to_upper(ref)}/{stringr::str_to_upper(alt)} for the region {chrom}:{start}-{stop}")
  # Format and run initial query
  .json_res <- httr::GET(
    url = glue::glue("https://portaldev.sph.umich.edu/ld/genome_builds/{genome_build}/references/1000G/populations/{population}/variants"),
    query = list(
      correlation = metric,
      variant = glue::glue("{chrom}:{pos}_{stringr::str_to_upper(ref)}/{stringr::str_to_upper(alt)}"),
      chrom = chrom,
      start = start,
      stop = stop
    )
  )

  # Parse query
  .json_res_parsed <- .json_res %>% httr::content(as = "parsed", type = "application/json")

  # Convert to tibble
  .json_res_parsed_df <- .json_res_parsed$data %>%
    tibble::as_tibble() %>%
    tidyr::unnest(cols = everything())

  # If tibble is empty, reverse ref/alt
  if (dim(.json_res_parsed_df)[1] == 0) {
    # Message
    cli::cli_alert_info("Initial query failed - trying again with flipped Ref/Alt")
    # Format and run initial query
    .json_res <- httr::GET(
      url = glue::glue("https://portaldev.sph.umich.edu/ld/genome_builds/{genome_build}/references/1000G/populations/{population}/variants"),
      query = list(
        correlation = metric,
        variant = glue::glue("{chrom}:{pos}_{stringr::str_to_upper(alt)}/{stringr::str_to_upper(ref)}"),
        chrom = chrom,
        start = start,
        stop = stop
      )
    )

    # Parse query
    .json_res_parsed <- .json_res %>% httr::content(as = "parsed", type = "application/json")

    # Convert to tibble
    .json_res_parsed_df <- .json_res_parsed$data %>%
      tibble::as_tibble() %>%
      tidyr::unnest(cols = everything())
  }

  if (dim(.json_res_parsed_df)[1] == 0) { # If tibble still empty, reference variant is missing from LD panel
    cli::cli_abort("Reference variant is missing from the specified LD panel")
  } else {
    return(suppressMessages(.json_res_parsed_df %>% readr::type_convert()))
  }
}
