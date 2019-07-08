

#' a simple "length-aware" genotyping error model for microsatellites
#'
#' In this genotyping error model, mis-calls are more likely to nearby
#' allele lengths than distant ones, and larger alleles are more likely
#' to drop out than small or average-sized ones.
#'
#' @param L an element of the list created by \code{\link{long_markers_to_X_l_list}}. Such
#' an element basically holds the information at a single locus.  The idea here
#' is that every ge_mod_* function takes in an object like L, and then can
#' use any piece of information in it about alleles or genotypes to
#' configure a genotyping error model. In the case of microsatellites,
#' the names of the alleles must be the allele lengths, like "64" or "112".
#' @param miscall_rate the rate at which alleles are miscalled.
#' @param miscall_decay the parameter (or a geometric distribution) that determines
#' how quickly the probability of the miscall being to a certain allele length
#' decreases as you get further and further away from the true allele length.
#' The larger this values, the more quickly the probability decays.  For every
#' microsatellite allele step size further from the true allele size, the probability
#' increases by a factor of 1 - miscall_decay.  Thus, if miscall_decay = 1
#' then all errors are a step-size of one away, and if it is 0, then
#' the miscall is equally likely to occur to any other allele length.
#' @param dropout_rate the base rate at which alleles of "typical" size dropout.
#' @param dropout_scale_factor the rate of drop out increases for alleles
#' larger than the (frequency-weighted) average size of allele by a factor of
#' dropout_scale_factor times the z-score of the allele length.
#' @export
#'
ge_model_microsat1 <- function(L,
                               miscall_rate = 0.005,
                               miscall_decay = 0.8,
                               dropout_rate = 0.01,
                               dropout_scale_factor = 0.25
) {

  # first, compute the step-size as the minimum distance between any two alleles
  lengths <- as.integer(names(L$freqs))
  diffs <- outer(lengths, lengths, function(x,y) abs(x - y))
  diag(diffs) <- NA
  min_diff = min(diffs, na.rm = TRUE)

  step_size = min_diff

  # now, we compute the decay factors between all pairs of alleles,
  # and then we multiply that by the base-miscall rate to get the
  # matrix W of allele-to-allele miscall rates.
  steps <- diffs/step_size

  preW <- miscall_decay * (1 - miscall_decay) ^ (steps - 1)

  # and those need to be normalized so that the off-diagonals sum to miscall_rate
  W <- t(apply(preW, 1, function(x) x * miscall_rate / sum(x, na.rm = TRUE)))
  # and we put 1 - miscall_rate on the diagonal
  diag(W) <- 1 - miscall_rate
  # and we MUST put the allele lengths on there as the names
  rownames(W) <- names(L$freqs)
  colnames(W) <- names(L$freqs)



  # Now, let's make the vector D of dropout rates.
  # first get the mean allele length
  mean_len <- sum(lengths * L$freqs)

  # and then we get the variance and the standard deviation.
  # Note that we have the freqs of each lengths, but not the sample
  # sizes, so we can't do the n-1 sort of correction for variance.
  var_len <- sum( L$freqs * (lengths - mean_len) ^ 2)
  sd_len <- sqrt(var_len)

  # now compute the z-scores
  zscores <- (lengths - mean_len) / sd_len

  # now we can compute the vector D of dropout probabilities
  D <- rep(dropout_rate, length(L$freqs))
  D[zscores > 0] <- dropout_rate * (1 + zscores[zscores > 0] * dropout_scale_factor)
  names(D) <- names(L$freqs)

  # And now we create the matrix C by pushing all these into the function
  # that combines allelic mis-calls and dropouts. Note that by setting Ws
  # to W here we assume that the occurrence of a dropout does not affect the
  # rate at which mis-calls occur.
  ret <- combine_miscalls_and_dropouts(geno_names = names(L$geno_freqs),
                                D = D,
                                W = W,
                                Ws = W)

  ret

}


