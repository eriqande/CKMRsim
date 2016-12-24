

#' compute the all-locus genotype pair probs for pairs simulated given specified relationships
#'
#' After you have created a ckmr object with \code{\link{create_ckmr}}, then
#' it is time to simulate multilocus genotype pairs under all the
#' relationships that you want to simulate from, and compute the likelihood of those
#' relationships under different relationship hypotheses.
#' This function does that.
#' @param C the ckmr object upon which to base the simulations.
#' @param froms a vector of names of the relationship IDs to simulate from. For each relationship ID in
#' froms, genotype values will get simulated from the Y_l_true values in C.
#' @param tos a vector of names of the relationship IDs to calculate the genotype log probs of the simulated
#' genotypes from.  Genotype log probs are calculated using the Y_l matrices.
#' @param reps An integer. The number of pairs to simulate for each relationship in \code{froms}.
#' Default is 10,000.
#' @param unlinked A logical indicating whether to simulate the markers as unlinked.  By default this is
#' TRUE.  If FALSE, then genotypes at linked markers will be simulated using the program MENDEL, genotyping
#' errors will be applied to them, and the Q_ij values themselves will still be computed under the assumption
#' of no linkage. However, they will be simulated under the no-linkage model for relationships "U", "PO", and
#' "MZ", because, in the absence of LD, related pairs under those relationships are not affected by
#' physical linkage.
#' @param pedigree_list  If you specify unlinked == FALSE, then you have to supply a pedigree_list.
#' @param forceLinkagePO If you really want to force simulation to be done under physical linkage for the PO case
#' (perhaps to verify that you get the same result as with unlinked.  Pass TRUE to this while \code{unlinked} is
#' FALSE.)
#' @param miss_mask_mat  A logical matrix with length(YL) columns and reps rows.  The (r,c)-th is TRUE if
#' the c-th locus should be considered missing in the r-th simulated sample.  This type of specification
#' lets the user simulate either a specific pattern of missingness, if desired, or to simulate patterns of
#' missing data given missing data rates, etc. (Need to make another function to do that.)
#' @param rando_miss_wts weights to be given to different loci that influence whether they
#' will be one of the rando_miss_n missing loci in any iteration.  These will be recycled
#' (or truncated) to have length equal to the number of loci, and they will be normalized
#' to sum to one as appropriate (so you can provide
#' them in unnormalized form.) The idea of this is to be able to use observed rates of
#' missingness amongst loci to mask some loci as missing.  Given as a comma-delimited
#' string in column "rando_miss_wts" in the output.
#' @param rando_miss_n a single number less than the number of loci.  Each iteration,
#' rando_miss_n loci will be considered missing, according to the rando_miss_wts.  This
#' let's you get a sense for how well you will do, on average, with a certain number of
#' missing loci.

simulate_Qij <- function(C, froms, tos, reps = 10^4, unlinked = TRUE, forceLinkagePO = FALSE, pedigree_list = NULL,
                         miss_mask_mat = NULL,
                         rando_miss_wts = NULL,
                         rando_miss_n = 0) {
  if(unlinked == FALSE && is.null(pedigree_list)) stop("If you choose unlinked == TRUE, you must supply a pedigree_list")
  if(unlinked == TRUE) {
    ret <- simulate_and_calc_Q(C$loci, reps = reps, froms = froms, tos = tos, miss_mask_mat = miss_mask_mat,
                               rando_miss_wts = rando_miss_wts, rando_miss_n = rando_miss_n)
  } else {
    ret <- simulate_and_calc_Q(YL = C$loci, reps = reps, froms = froms, tos = tos, df = C$orig_data, pedigrees = pedigree_list,
                               forceLinkagePO = forceLinkagePO, miss_mask_mat = miss_mask_mat,
                               rando_miss_wts = rando_miss_wts, rando_miss_n = rando_miss_n)
  }


  Qij_class(ret, unlinked, forceLinkagePO, miss_mask_mat, rando_miss_wts, rando_miss_n)
}
