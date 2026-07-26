
#### Import the pipe operator from magrittr ####
#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#### Import functions from packages ####
#' @importFrom dplyr n n_distinct left_join mutate rename select arrange bind_rows count desc filter last
#' @importFrom ggplot2 aes geom_abline geom_col geom_point ggplot xlab ylab
#' @importFrom purrr map2
#' @importFrom readr write_rds write_tsv
#' @importFrom stats quantile rmultinom runif sd setNames
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_wider separate unnest
#' @importFrom utils read.table unzip write.table




#### Declare names of columns  to keep CRAN checks from barking  ####
# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if(getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "AbsoluteIndex",
      "AlleIdx",
      "AlleLine",
      "Allele",
      "Chrom",
      "D1a",
      "D1a_1",
      "D1a_2",
      "D2a",
      "D2a_1",
      "D2a_2",
      "FNR",
      "FPR",
      "Freq",
      "Intercept",
      "Kid",
      "Lambda_star",
      "LocIdx",
      "LocLine",
      "LocScrub",
      "Locus",
      "Ma",
      "Manew",
      "Pa",
      "Panew",
      "Pos",
      "Qijs_linked",
      "Qijs_unlinked",
      "Rep",
      "Sex",
      "Sexnew",
      "TwinStatus",
      "a1",
      "a2",
      "alle1",
      "alle2",
      "alleidx",
      "data",
      "desc",
      "ender",
      "estimate",
      "g",
      "genostr",
      "gp",
      "h",
      "hp",
      "id",
      "impwt",
      "indiv_1",
      "indiv_2",
      "kappas",
      "lambda",
      "list_name",
      "locidx",
      "log10_FPR",
      "name",
      "newfreq",
      "numnonmissingloci",
      "pedname",
      "mc",
      "num_non_miss",
      "num_non_missing_loci",
      "se",
      "simplm",
      "term",
      "tidy",
      "xxx",
      "gnames",
      "1",
      "2",
      "D1_indiv",
      "D2_indiv",
      "GenoIdx",
      "Indiv",
      "NumA",
      "chrom",
      "gene_copy",
      "idx",
      "ind",
      "ind1",
      "ind2",
      "loc_name",
      "num_bases",
      "num_markers",
      "num_loc",
      "num_mismatch",
      "pos",
      "scaled_length",
      "value"
    )
  )
}
