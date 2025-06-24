#' implements a simple true-genotype-independent genotyping error model
#'
#' Note that this is sort of silly in the context of parentage inference
#' because for non-parentals, it can lead to probability of zero for
#' Mendelian incompatibilities.  Probably better to use the ge_model_microhap1().
#'
#' @param L required locus specific information
#' @param epsilon the rate at which genotypes are incorrectly observed.
#' @export
#' @examples
#' L <- list(
#'  freqs = c(A = 0.6, B = 0.4),
#'  geno_freqs = c(`A / A` = 0.36, `A / B` = 0.48, `B / B` = 0.16)
#' )
#' ge_model_TGIE(L, m = 0.185)
ge_model_pure_het_miscall <- function(L, m) {
  # first make a matrix G x G matrix of the genotype frequencies
  gmat <- matrix(rep(L$geno_freqs, length(L$geno_freqs)), byrow = TRUE, nrow = length(L$geno_freqs))

  # since this is a pure het miscall model, we don't have any
  # errors in allelic types from the reads.

  # put the genotype names on for the row and column names since that can be
  # handy down the road.
  rownames(gmat) <- names(L$geno_freqs)
  colnames(gmat) <- names(L$geno_freqs)

  # Now we make the W matrix, which is just diagonal ones because there
  # is no allelic miscall in this model
  W <- matrix(0, nrow = length(L$freqs), ncol = length(L$freqs))
  diag(W) <- 1
  rownames(W) <- names(L$freqs)
  colnames(W) <- names(L$freqs)

  # now we can compute the vector D of dropout probabilities.  These are m/2
  # because then we have a total het miscall rate of m
  D <- rep(m/2, length(L$freqs))
  names(D) <- names(L$freqs)

  # And now we create the matrix C by pushing all these into the function
  # that combines allelic mis-calls and dropouts. Note that by setting Ws
  # to W here we assume that the occurrence of a dropout does not affect the
  # rate at which mis-calls occur.
  ret <- combine_miscalls_and_dropouts(
    geno_names = names(L$geno_freqs),
    D = D,
    W = W,
    Ws = W
  )

  ret
}

