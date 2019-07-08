#' implements the generalized miscall + dropout model
#'
#' This is the version into which you pass in your own matrices of allele
#' call probabilities W (corresponding to omega) and Ws (omega-star)
#' and vectors of allele-specific dropout rates. Notation follows the paper
#' for the most part. So see it for details.
#'
#' The old version (general_allele_based_geno_err_model) required that
#' the variables g, gp, h, and gp all
#' be input by the user.  This version just uses the names attribute of the
#' geno_freqs variable to take care of that.  In other words, you pass
#' it a character vector of genotype names (constructed by CKMRsim using the
#' standard allele separate, " / ", and this function creates the g, gp, h, and hp
#' vectors that are needed.)
#'
#' Note that this is vectorized
#' over g, gp, h, and hp.  Note that g, gp, h, and hp can be character vectors
#' giving the allele names, so long as the vector D has names which are the allele
#' names and W and Ws have rownames and colnames which are the allele names.
#' @param geno_names  the names of the genotypes in the order that CKMRsim has
#' stored them.
#' @param D vector of allele-specific dropout rates (delta) in the paper. Must be named
#' with the names of the alleles, and should be in the correct order (that which CKMRsim uses).
#' @param W the matrix omega of allelic call probabilities. Must have rownames and colnames that
#' are the allele names, and have to be in the correct order.
#' @param Ws the matrix omega* of allelic call probabilities.  Gotta have rownames and colnames.
#' @return This function returns matrix C, which is the genotype-call probability matrix.  The element
#' in the s-th row and the t-th column gives the probability that a true genotype
#' of index s is called as a genotype of index t.  This matrix is returned with rownames and colnames
#' which are the genotype names using the standard CKMRsim alelle separator of " / ".
#' @export
#' @examples
#' # simple example with 2-allele locus.
#' # first, look at the relevant information you would have about a locus,
#' # namely the allele frequencies and the genotype frequencies
#' example_L_biallelic
#'
#' # Now, we create a matrix of allele-to-allele mis-call rates.
#' # Let's say that with probability 0.01 you miscall a true A as a B
#' # and with probability 0.03 you miscall a true B as an A.  We set
#' # these as different just because it provides a better illustration here.
#' W <- matrix(c(0.99, 0.01, 0.03, 0.97), ncol = 2, byrow = TRUE)
#' rownames(W) <- c("A", "B")
#' colnames(W) <- c("A", "B")
#'
#' # Have a look at that:
#' W
#'
#' # Ws is a matrix of similar form that specifies the miscall rates given that
#' # a dropout occurred.  We will set it to be the same as W for this example:
#' Ws <- W
#'
#' # Now we need a vector of allele-specific dropout rates. Let's say that
#' # A is very likely to drop out (0.05) and B, not so much (0.02)
#' D <- c(A = 0.05, B = 0.02)
#'
#' # See what the resulting genotype-call probability matrix looks like:
#' combine_miscalls_and_dropouts(names(example_L_biallelic$geno_freqs), D, W, Ws)
#'
#' # now, let's see what it looks like if it is only the dropouts, no miscalls
#' W2 <- matrix(0, ncol = 2, nrow = 2)
#' diag(W2) <- 1
#' dimnames(W2) <- dimnames(W)  # recall, you must have row and column names on these
#' Ws2 <- W2
#'
#' combine_miscalls_and_dropouts(names(example_L_biallelic$geno_freqs), D, W2, Ws2)
#'
combine_miscalls_and_dropouts <- function(geno_names, D, W, Ws) {

  # we will first make a tibble of g gp h and hp.  This just has one row
  # for every possible combination of alleles, true and observed, in different
  # genotypes, and it is ordered with the true genotypes cycling fastest. This
  # is the order you get if you flatten the matrix C.
  allele_separator = " / "  # define a variable for it in case we decide to change it later
  tmp <- tibble(gnames = geno_names) %>%
    separate(gnames, into = c("g", "gp"), sep = allele_separator)

  g <- rep(tmp$g, nrow(tmp))
  gp <- rep(tmp$gp, nrow(tmp))
  h <- rep(tmp$g, each = nrow(tmp))
  hp <- rep(tmp$gp, each = nrow(tmp))

  # now, we make matrices of those for subscripting (see the paper for details of why/how this works)
  gh = cbind(g, h)
  gph = cbind(gp, h)
  ghp = cbind(g, hp)
  gphp = cbind(gp, hp)


  ret_vec <- (h == hp) * (D[g] * Ws[gph] + D[gp] * Ws[gh]) +
    (1 - D[g] - D[gp]) * (W[gh] * W[gphp] + (h != hp) * W[gph] * W[ghp])


  ret <- matrix(ret_vec, nrow = length(geno_names))
  colnames(ret) <- geno_names
  rownames(ret) <- geno_names

  ret
}

