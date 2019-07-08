

#' implements a simple true-genotype-independent genotyping error model
#'
#' @param L required locus specific information
#' @param epsilon the rate at which genotypes are incorrectly observed.
#' @export
#' @examples
#' L <- list(
#'  freqs = c(A = 0.6, B = 0.4),
#'  geno_freqs = c(`A-A` = 0.36, `A-B` = 0.48, `B-B` = 0.16)
#' )
#' ge_model_TGIE(L, epsilon = 0.02)
ge_model_TGIE <- function(L, epsilon = 0.01) {
  # first make a matrix G x G matrix of the genotype frequencies
  gmat <- matrix(rep(L$geno_freqs, length(L$geno_freqs)), byrow = TRUE, nrow = length(L$geno_freqs))

  # Now, we make the diagonals 1 - epsilon exactly and rescale
  # everything else to sum to epsilon.
  diag(gmat) <- 0
  gmat <- t(apply(gmat, 1, function(x) x * (epsilon/sum(x))))
  diag(gmat) <- 1 - epsilon

  # put the genotype names on for the row and column names since that can be
  # handy down the road.
  rownames(gmat) <- names(L$geno_freqs)
  colnames(gmat) <- names(L$geno_freqs)

  # return gmat, which is the matrix C.
  gmat
}

