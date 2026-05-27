# Query LocusZoom API for Recombination Data

This function queries the LocusZoom API to retrieve recombination data
for a specified genomic region and returns the result as a tibble.

## Usage

``` r
recomb_extract_locuszoom(chrom, start, end, genome_build = "GRCh37")
```

## Arguments

- chrom:

  A numeric value specifying the chromosome (e.g., 1, 2, ..., 22, 23 for
  X, 24 for Y)

- start:

  An integer specifying the start position of the region of interest

- end:

  An integer specifying the end position of the region of interest

- genome_build:

  A character string specifying the genome build (default: "GRCh37")

## Value

A tibble containing the parsed recombination data from the LocusZoom API

## Examples

``` r
if (FALSE) { # \dontrun{
result <- recomb_locuszoom(chrom = 1, start = 1000, end = 150000)
print(result)

# Using a different genome build
result_grch38 <- recomb_locuszoom(chrom = 1, start = 1000, end = 150000, genome_build = "GRCh38")
} # }
```
