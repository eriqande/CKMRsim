

# these are functions that have to do with flattening a ckmr object into a single
# vector (and associated genotype counts, and other information)
# that lets us index into it from the genotypes in observed pairs of individuals,
# so that we can quickly and easily compute the observed Logls from within Rcpp.

#' flatten a from/to pair from a ckmr object (see \code{\link{ckmr_class}}into a long vector of probabilities
#'
#' blah
#' @param CK an object of class ckmr
#' @param rel The string name of the relationship you want to extract probs for.
#' @param useYl_true set to TRUE if you want to use the Y_l_true field instead of Y_l.
#' @param unname if TRUE (the default) this strips the names off the results.
#' @export
flatten_ckmr <- function(CK, rel, useYl_true = FALSE, unname = TRUE) {
  if(useYl_true == TRUE) {
    mats <- lapply(CK$loci, function(x) x$Y_l_true[[rel]])
  } else {
    mats <- lapply(CK$loci, function(x) x$Y_l[[rel]])
  }

  ret <- list()
  ret$nGenos <- unlist(lapply(mats, nrow))
  nGenoPairs <- unlist(lapply(mats, length))  # the number of genotypes pairs at each locus
  tmp <- cumsum(nGenoPairs)
  tmp <- c(0, tmp)
  names(tmp) <- names(tmp)[-1]
  ret$base0_locus_starts <- tmp[-length(tmp)]
  ret$probs <- unlist(lapply(mats, function(x) as.numeric(x)))

  if(unname == TRUE) {
    ret$nGenos <- unname(ret$nGenos)
    ret$base0_locus_starts <- unname(ret$base0_locus_starts)
    ret$probs <- unname(ret$probs)
  }

  # now, with that information I should be able to pick out the prob given the index of the locus and the of
  # the genotype at indiv1 and indiv2.
  ret
}
