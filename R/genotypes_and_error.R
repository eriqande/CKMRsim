

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




#' convert a  long data frame of marker allele frequencies to a list of genotype frequency vectors
#'
#' Given a data frame D in the format of long_markers, this function creates the vectors G_l for
#' each locus. These are vectors of genotype frequencies with the genotypes ordered as in the
#' function \code{\link{index_ab}}.
#' @param  D a data frame in the format of \code{\link{long_markers}}.  Monomorphic loci should have
#' been removed and the whole thing run through \code{\link{reindex_markers}} before passing it
#' to this function.
#' @return This returns a named list.  The names of the components are the Chrom.Locus of the
#' marker in D.
#' @examples
#' long_markers_to_G_list(long_markers[1:18, ])
long_markers_to_G_list <- function(D) {
  # here is a quick function to compute the genotype freqs in the correct order
  # given a vector of allele freqs
  gfreqs <- function(d) {
    a <- d$Freq
    b <- outer(a, a, function(x,y) 2 * x * y)
    diag(b) <- diag(b) / 2
    b[lower.tri(b, diag = TRUE)]
  }

  D2 <- D %>%
    mutate(gname = paste(Chrom, Locus, sep = "."))

  split(D2, D2$gname) %>%
    lapply(., gfreqs)

}
