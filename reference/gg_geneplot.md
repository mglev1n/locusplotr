# Plot genes located within a genomic region of interest

Returns a ggplot containing the genes within a specified genomic region.
The function uses database connections to EnsDb.Hsapiens.v75
(hg19/GRCh37) or EnsDb.Hsapiens.v86 (hg38/GRCh38) to identify genes
within the specified region, and uses the ggbio package to create the
plot.

## Usage

``` r
gg_geneplot(chr, start, end, genome_build = "GRCh38", max_levels = 5)
```

## Arguments

- chr:

  Integer - chromosome

- start:

  Integer - starting position for region of interest

- end:

  Integer - ending position for region of interest

- genome_build:

  Character - genome build - one of "GRCh37" or "GRCh38"

- max_levels:

  Integer - maximum number of levels for gene tracks

## Value

A ggplot object containing a plot of genes within the region of interest

## Examples

``` r
if (FALSE) { # \dontrun{
gg_geneplot(1, 170054349 - 1e6, 170054349 + 1e6, "GRCh37")
} # }
```
