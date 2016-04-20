# functions for simulating genotypes and calculating the Q_{ij}'s and the Lambdas




#' compute all-locus genotype pair probs for pairs simulated given specified relationships
#'
#' Once the Y_l and Y_l_true matrices have been inserted into the list of locus computations, then
#' it is time to simulate multilocus genotype pairs under a variety of relationships,
#' and then compute the log probability of them assuming various relationships.
#' @param YL the list that holds the Y_l and Y_l_true matrices
#' @param reps number of pairs to simulate for each relationship in froms.  This recycles
#' to be the same length as froms, but you can specify different numbers of reps for each
#' of the relationships.
#' @param froms a vector of names of the relationship IDs to simulate from. If this is NULL then
#' the function will just simulate from all the values that have had Y_l_true matrices computed for them
#' in the first component of YL.  Genotype values get simulated from the Y_l_true values
#' @param tos a vector of names of the relationship IDs to calculate the genotype log probs of the simulated
#' genotypes from.  Genotype log probs are calculated using the Y_l matrices. If this is NULL then
#' the function will just compute probs for all the relationships that have had Y_l matrices computed for them
#' in YL.
#' @param rando_miss_n NOT IMPLEMENTED YET how many loci to mask at random on each
#' iteration.  If this is a vector then it does the whole analysis for each
#' one of the values.  Corresponds to column "rando_miss_n" in the output.
#' Default is 0. If any value is greater than or equal to the total number of
#' loci available, it is removed.
#' @param rando_miss_wts weights NOT IMPLEMENTED YET to be given to different loci that influence whether they
#' will be one of the rando_miss_n missing loci in any iteration.  These will be recycled
#' (or truncated) to have length equal to the number of loci (or to the \code{first_n_loci} if that
#' is in effect), and they will be normalized to sum to one as appropriate (so you can provide
#' them in unnormalized form.)  The idea of this is to be able to use observed rates of
#' missingness amongst loci to mask some loci as missing.  Given as a comma-delimited
#' string in column "rando_miss_wts" in the output.

#' @return This returns a list with components that are the relationships that were simulated
#' from.  Inside each of those components is a list with components referring to the relationships
#' that had genotype log probabilities calculated (the "tos"). The contents of each is a vector
#' of length reps which are the log probability of the multilocus genotypes of that each of reps
#' simulated pairs.
#' @examples
#' example(insert_Y_l_matrices)
#' mh_sacQ_example <- simulate_and_calc_Q(mh_yl_example, 10^4, froms = c("U", "PO"), tos = c("U", "PO"))
#'
#' \dontrun{
#' # and we could plot the distribution of those if we wanted to
#' m <- mh_sacQ_example
#' D <- data.frame(U = m$U$PO - m$U$U, PO = m$PO$PO - m$PO$U)
#' DD <- tidyr::gather(D, Relat, logL)
#' library(ggplot2)
#' ggplot(DD, aes(x = logL, fill = Relat)) + geom_density(alpha = 0.6)
#' }
simulate_and_calc_Q <- function(YL, reps = 10^4, froms = NULL, tos = NULL) {
  if(is.null(froms)) {
    froms <- names(YL[[1]]$Y_l_true)
  }
  if(is.null(tos)) {
    tos <- names(YL[[1]]$Y_l)
  }

  # deal with recycling the reps as necessary
  reps <- rep(reps, length.out = length(froms))
  names(reps) <- froms

  # name these so that if we lapply over them, we get the names in the result
  names(froms) <- froms
  names(tos) <- tos

  # cycle over the relationships to simulate from
  lapply(froms, function(r) {
    # then simulate the indexes of the genotype pairs that are simulated at each locus
    gp <- lapply(YL, function(y) { # this is cycling over loci
      sample.int(n = length(y$Y_l_true[[r]]), size = reps[r], replace = TRUE, prob = y$Y_l_true[[r]])
    })

    # now with those in hand we cycle over the tos relationships and compute the log probs
    # for them and then cbind them and rowSum them...
    lapply(tos, function(t) {
      sapply(names(YL), function(y) {
        log(YL[[y]]$Y_l[[t]])[gp[[y]]]   # right about here is where I think I can add randomly missing loci, etc.
      }) %>%
        rowSums
    })
  })
}
