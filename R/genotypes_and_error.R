

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
#' Given a data frame D in the format of long_markers, this function creates the matrix X_l for
#' each locus, given all the values of kappa supplied as rows in a matrix like that in
#' \code{\link{kappas}}. The rows and columns refer to genotypes that are in
#' the order given in the
#' function \code{\link{index_ab}}.
#' @param  D a data frame in the format of \code{\link{long_markers}}.  Monomorphic loci should have
#' been removed and the whole thing run through \code{\link{reindex_markers}} before passing it
#' to this function.
#' @param kappa_matrix  A matrix like that in the supplied data \code{\link{kappas}}.  It should have
#' three columns, and rownames giving the abbreviated name of the pairwise relationship.  Each row
#' is a three-vector which are the Cotterman coefficients for the relationship. The first element is the
#' probability that the pair shares 0 gene copies IBD, the second is the prob that they share 1 gene
#' copy IBD, and the third is the prob that they share 2 gene copies IBD, all assuming no inbreeding.
#' @return This returns a named list.  The names of the components are the Chrom.Locus.Pos of the
#' marker in D. The contents of each list component is a list that includes the allele frequencies
#' (as a vector named with the Allele names), and also another list of matrices with nA * (nA+1) / 2 rows and
#' columns.  The rows and columns of this matrix are named by the genotypes.
#' @examples
#' data(kappas)
#' lm_example <- long_markers_to_X_l_list(long_markers[1:18,], kappa_matrix = kappas)
#' mh_example <- long_markers_to_X_l_list(microhaps, kappa_matrix = kappas)
long_markers_to_X_l_list <- function(D, kappa_matrix) {

  D2 <- D %>%
    mutate(gname = paste(Chrom, Locus, Pos, sep = "."))

  # this long expression splits D on the Chrom.Locus.Pos and lapplies over
  # the pieces and computes the X_l matrices for each locus, and then just
  # returns this
  split(D2, factor(D2$gname, levels = unique(D2$gname))) %>%
    lapply(., function(x) {

      # a list for output
      loc_ret <- list()

      # for the dimnames of the X_l matrices
      geno_names <- expand.grid(gp = x$Allele, g = x$Allele) %>%
        filter(as.integer(factor(g, levels = x$Allele)) <= as.integer(factor(gp, levels = x$Allele))) %>%
        mutate(name = paste(g, gp, sep = "-")) %>%
        select(name) %>%
        unlist %>%
        unname

      # store the allele frequencies with the allele names
      loc_ret$freqs <- x$Freq
      names(loc_ret$freqs) <- x$Allele

      loc_ret$X_l <- lapply(1:nrow(kappa_matrix), function(K) {
        tmp <- make_matrix_X_l(p = x$Freq, kappa_matrix[K, ])
        dimnames(tmp) <- list(geno_indiv_1 = geno_names, geno_indiv_2 = geno_names)
        tmp
      })

      names(loc_ret$X_l) <- rownames(kappa_matrix)
      loc_ret
    })
}




#' insert the matrix C_l into each element of the list of X_l matrices
#'
#' After you have gotten the X_l matrices using \code{\link{long_markers_to_X_l_list}}
#' you can add to the component for each locus the matrix C_l which has, as its (s,t)-th
#' entry, the probability that the genotype t is observed given that the true genotype
#' is s. I am going to a add a lot more parameters to this as soon as I figure out how to
#' pass in genotyping error models and their associated parameters. I would like to do that in
#' a way that easily let's people define their own functions.  But for now it just
#' does the microhaplotype error model with default parameters.
#' @param XL a list of the loci like that created using \code{\link{long_markers_to_X_l_list}}.
#' The key thing that each list component needs is the named vector freqs of the allele frequencies.
#' The functions that compute genotyping error use the names of the allele to compute the
#' probabilities of the observed genotype given the true genotype.
#' @param ... extra arguments (after haps) to be passed to microhaplotype_geno_err_matrix
#' #' @examples
#' example(long_markers_to_X_l_list)
#' mh_cl_example <- insert_C_l_matrices(mh_example)
insert_C_l_matrices <- function(XL, ...) {
  lapply(XL, function(x) {
    if(is.null(names(x$freqs))) stop("No names attribute on the allele frequencies")
    x$C_l <- microhaplotype_geno_err_matrix(haps = names(x$freqs), ...)
    x
  })
}




#' compute the Y_l matrices for each locus in a list
#'
#' Once you have the X_l matrices for each kappa, and the C_l matrices
#' have been inserted into the list, then this function cycles over the loci
#' and the kappa values and does the matrix multiplication to give you the
#' matrix Y_l.
#' @param L the list of X_l and C_l matrices
#' @examples
#' example(insert_C_l_matrices)
#' mh_yl_example <- insert_Y_l_matrices(mh_cl_example)

insert_Y_l_matrices <- function(L) {
  lapply(L, function(a) {
    a$Y_l <- lapply(a$X_l, function(b) {
      t(b %*% a$C_l) %*% t(a$C_l)
    })
    a
  })
}
