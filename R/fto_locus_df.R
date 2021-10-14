#' Example genome wide association study results at the FTO locus
#'
#' A dataset containing genome wide association study stummary results at the FTO locus on chromosome 16
#'
#' @format A data frame with 19119 rows and 9 variables:
#' \describe{
#'   \item{chromosome}{chromosome}
#'   \item{position}{position on chromosome}
#'   \item{rsid}{rsid number from dbSNP}
#'   \item{effect_allele}{effect allele (or reference allele)}
#'   \item{other_allele}{non-effect allele (or alternate allele)}
#'   \item{eaf}{effect allele frequenct (or minor allele frequency)}
#'   \item{effect}{effect of allele on outcome}
#'   \item{std_err}{standard error of the effect estimate}
#'   \item{p_value}{p-value of effect estimate}
#' }
"fto_locus_df"
