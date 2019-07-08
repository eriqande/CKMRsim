
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
#' @importFrom dplyr n n_distinct rename
#' @importFrom stats quantile runif sd setNames
#' @importFrom tibble tibble
#' @importFrom tidyr separate
#' @importFrom utils read.table write.table




#### Declare names of columns  to keep CRAN checks from barking  ####
# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if(getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "AlleIdx",
      "AlleLine",
      "Allele",
      "Chrom",
      "FNR",
      "Freq",
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
      "Rep",
      "Sex",
      "Sexnew",
      "TwinStatus",
      "a1",
      "a2",
      "alle1",
      "alle2",
      "alleidx",
      "desc",
      "ender",
      "g",
      "genostr",
      "gp",
      "h",
      "hp",
      "id",
      "impwt",
      "kappas",
      "lambda",
      "list_name",
      "locidx",
      "name",
      "newfreq",
      "pedname",
      "xxx",
      "gnames"
    )
  )
}
