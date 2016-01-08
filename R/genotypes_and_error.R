

#' a function to compute the genotype index from alleles a and b
#'
#' This implements the function to index genotypes so that we can put
#' them into matrices.  The paper mentions this indexing, in the section
#' "Implementation."
#' @param The first allele in the genotype
#' @param The second allele in the genotype
#' @param A the total number of alleles at the locus in the population.
#' @details If a is not an allele with index lower than b, then this function
#' automatically reorients the alleles so that it is. This is vectorized and
#' so should work over long vectors of alleles.
#' @examples
#' # Create some data to test/demonstrate this:
#' A <- 7
#' all_genos <- expand.grid(1:A, 1:A)[,c(2,1)] %>%
#'  setNames(c("a", "b"))
#'
#' # when a < b
#' right_ord <- all_genos %>%
#'  filter(a <= b) %>%
#'  mutate(Idx = 1:length(a),
#'       Indab = index_ab(a, b, A))
#'
#' # when a > b
#' wrong_ord <- all_genos %>%
#'  filter(a <= b) %>%
#'  mutate(Idx = 1:length(a),
#'       Indab = index_ab(b, a, A))
#'
#' # then check to make sure it all checks out
#' all(right_ord$Idx == right_ord$Indab)
#' all(wrong_ord$Idx == wrong_ord$Indab)
index_ab <- function(a, b, A) {
  newa <- ifelse(a > b, b, a)
  newb <- ifelse(a > b, a, b)

  a <- newa
  b <- newb

  2 + (a - 1) * (A + 2) - ( a * (a + 1) / 2) + (b - a)
}




#' convert a  long data frame of marker allele frequencies to a list of X_l matrices
#'
#' Given a data frame D in the format of long_markers, this function creates the matrix_X_l for
#' each locus, given the value of kappa. The rows and columns refer to genotypes that are in
#' the order given in the
#' function \code{\link{index_ab}}.
#' @param  D a data frame in the format of \code{\link{long_markers}}.  Monomorphic loci should have
#' been removed and the whole thing run through \code{\link{reindex_markers}} before passing it
#' to this function.
#' @param kappa A 3-vector giving the cotterman coefficients for a given relationship (like
#' (0.25, 0.5, 0.25) for full siblings, or (0,1,0) for parent-offspring.)
#' @return This returns a named list.  The names of the components are the Chrom.Locus of the
#' marker in D. The contents of each list component is a matrix with nA * (nA+1) / 2 rows and
#' columns
#' @examples
#' long_markers_to_X_l_list(long_markers[1:18,], c(0.25, 0.5, 0.25))
long_markers_to_X_l_list <- function(D, kappa) {

  D2 <- D %>%
    mutate(gname = paste(Chrom, Locus, sep = "."))

  split(D2, D2$gname) %>%
    lapply(., function(x) make_matrix_X_l(p = x$Freq, kappa = kappa))
}




#' a function to compute the X_l matrices for all loci at each desired kappa value
#'
#' still working on this
#' @inheritParams long_markers_to_X_l_list
#' @param kappa_matrix  A matrix like that in the supplied data \code{\link{kappas}}
#' @return this returns a named list of length equal to the number of rows in kappa_matrix.
#' The names are the rownames of the kappa_matrix.  Each component is itself a list named
#' by locus with each component holding the X_l matrix.
#' @examples
#' data(kappas)
#' compute_all_X_l_matrices(long_markers[1:18,], kappas)
compute_all_X_l_matrices <- function(D, kappa_matrix) {
  if(!is.matrix(kappa_matrix)) stop("kappa_matrix is not a matrix")
  if(is.null(rownames(kappa_matrix))) stop("kappa_matrix does not have rownames")
  if(ncol(kappa_matrix) != 3) stop("kappa_matrix should have 3 columns")

  ret <- lapply(rownames(kappa_matrix), function(x) {
    k <- kappa_matrix[x, ]
    long_markers_to_X_l_list(D, k)
  })
  names(ret) <- rownames(kappa_matrix)
  ret
}
