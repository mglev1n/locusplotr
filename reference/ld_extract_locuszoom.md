# Fetch linkage disequilibrium statistics for a locus of interest

This function allows the user to extract linkage disequilibrium
statistics between a variant of interest and surrounding variants within
a contiguous genomic region. This function uses the University of
Michigan LocusZoom API (<https://portaldev.sph.umich.edu/>) to obtain LD
information, and allows the user to specify genome-build and ancestry of
interest.

## Usage

``` r
ld_extract_locuszoom(
  chrom,
  pos,
  ref,
  alt,
  start,
  stop,
  genome_build = "GRCh37",
  population = "ALL",
  metric = "rsquare"
)
```

## Arguments

- chrom:

  Integer - chromosome of reference variant

- pos:

  Integer - position of reference variant

- ref:

  Character - reference allele (or effect allele) for reference variant

- alt:

  Character - alternate allele (or non-effect allele) for reference
  variant

- start:

  Integer - starting position of range of interest

- stop:

  Integer - ending position of range of interest

- genome_build:

  Character - Genome build - one of "GRCh37" or "GRCh38"

- population:

  Character - 1000G genetic superpopulation - one of "ALL", "AFR",
  "AMR", "EAS", "EUR", "SAS"

- metric:

  Character - one of "r", "rsquare", or "cov", referring to the
  correlation statistic of interest

## Value

A tibble containing each variant within the supplied range surrounding
the variant of interest, with the requested linkage disequilibrium
information with respect to the variant of interest

## Examples

``` r
if (FALSE) { # \dontrun{
ld_extract_locuszoom(chrom = 16, pos = 53830055, ref = "C", alt = "G", start = 53830055 - 5e5, stop = 53830055 + 5e5, genome_build = "GRCh37", population = "ALL", metric = "rsquare")
} # }
```
