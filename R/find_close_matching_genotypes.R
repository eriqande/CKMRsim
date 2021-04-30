#' Return every pair of individuals that mismatch at no more than max_miss loci
#'
#' This is used for identifying duplicate individuals/genotypes in large
#' data sets. I've specified this in terms of the max number of missing loci because
#' I think everyone should already have tossed out individuals with a lot of
#' missing data, and then it makes it easy to toss out pairs without even
#' looking at all the loci, so it is faster for all the comparisons.
#'
#' @param LG a long genotypes data frame
#' @param CK a ckmr object created from the allele frequencies computed from LG.
#' @param max_mismatch maximum allowable number of mismatching genotypes betwen the pairs.
#' @return a data frame with columns:
#' \describe{
#'   \item{indiv_1}{the id (from the rownames in S) of the firt member of the pair}
#'   \item{indiv_2}{the id (from the rownames in S) of the second individual of the pair}
#'   \item{num_mismatch}{the number of loci at which the pair have mismatching genotypes}
#'   \item{num_loc}{the total number of loci missing in neither individual}
#' }
#' @export
find_close_matching_genotypes <- function(LG, CK, max_mismatch) {

  S <- create_integer_genotype_matrix(LG = LG,
                                      AF = CK$orig_data)

  matchers <- pairwise_geno_id(S, max_miss = max_mismatch) %>%
    dplyr::arrange(num_mismatch) %>%
    dplyr::mutate(indiv_1 = rownames(S)[ind1],
           indiv_2 = rownames(S)[ind2])

  matchers
}
