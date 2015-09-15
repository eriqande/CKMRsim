



#' A list of pedigrees specifying different pairwise relationships for MORGAN's genedrop program.
#'
#' A list of pedigrees, each specified as a data frame with five columns.
#'
#'
#' @format A named list of pedigrees, each one specified in a data frame with five columns in the following order:
#' \describe{
#'   \item{Kid}{alphanumeric character specifier of the kid's name.  Must be <15 characters.}
#'   \item{Pa}{alphanumeric character specifier of the dad's name.  Must be <15 characters.}
#'   \item{Ma}{alphanumeric character specifier of the mom's name.  Must be <15 characters.}
#'   \item{Sex}{Integer specifier of sex of the kid: 0 - unknown, 1 - male, 2 - female.}
#'   \item{Observed}{Integer specifying if the kid is observed.  Typically there will be two observed
#'                   individuals: i and j. }
#' }
"pedigrees"




#' An example of how to specify markers mapped onto chromosomes
#'
#' This is a named list.  Each component is named with the number of the chromosome.
#' In this case, those numbers are: "4", "9", "14", "19", and "21". Each  component
#' is itself a named list that holds the position in centiMorgans and the allele frequencies
#' of the alleles at the marker.
#'
#'
#' @format A named list of chromosomes upon which mapped markers are available.  The component
#' corresponding to a chromome is itself a list which may or may not be named.  If it is named
#' then its names are the names of the loci. The loci should be in sorted order along the chromosome.
#' Each locus within a chromosome is represented as a list with the following components.
#' \describe{
#'   \item{pos}{Position in centiMorgans along the chromosome. This is a little nonstandard, but you
#'   can think of the first locus as having a certain rate of recombination with the very beginning of
#'   the chromosome, and then pos is the cumulative sum of the inter-locus recombination fractions. }
#'   \item{freqs}{The allele frequencies at the locus}
#' }
"markers_on_map"


