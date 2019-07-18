

#' a function to compute the genotype index from alleles a and b
#'
#' This implements the function to index genotypes so that we can put
#' them into matrices.  The paper mentions this indexing, in the section
#' "Implementation."
#' @param a The first allele in the genotype
#' @param b The second allele in the genotype
#' @param A the total number of alleles at the locus in the population.
#' @details If a is not an allele with index lower than b, then this function
#' automatically reorients the alleles so that it is. This is vectorized and
#' so should work over long vectors of alleles.
#' @export
#' @examples
#' # Create some data to test/demonstrate this:
#' A <- 7
#' all_genos <- expand.grid(1:A, 1:A)[,c(2,1)] %>%
#'  setNames(c("a", "b"))
#'
#' # when a < b
#' right_ord <- all_genos %>%
#'  dplyr::filter(a <= b) %>%
#' dplyr::mutate(Idx = 1:length(a),
#'       Indab = index_ab(a, b, A))
#'
#' # when a > b
#' wrong_ord <- all_genos %>%
#'  dplyr::filter(a >= b) %>%
#' dplyr::mutate(Idx = 1:length(a),
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
#' @param allele_separator This string is used to separate allele names in the names of genotypes.
#' By default it is " / ".  This should not get confused by any elements in actual allele names, as
#' long as the allele names don't have any spaces in them.  This gets used to parse the genotype
#' names down the road.
#' @return This returns a named list.  The names of the components are the Chrom.Locus.Pos of the
#' marker in D. The contents of each list component is a list that includes the allele frequencies
#' (as a vector named with the Allele names), and also another list of matrices with nA * (nA+1) / 2 rows and
#' columns.  The rows and columns of this matrix are named by the genotypes.
#' @export
#' @examples
#' data(kappas)
#' lm_example <- long_markers_to_X_l_list(long_markers[1:18,], kappa_matrix = kappas)
#' mh_example <- long_markers_to_X_l_list(microhaps, kappa_matrix = kappas)
long_markers_to_X_l_list <- function(D, kappa_matrix, allele_separator = " / ") {

  D2 <- D %>%
   dplyr::mutate(gname = paste(Chrom, Locus, Pos, sep = "."))

  # this long expression splits D on the Chrom.Locus.Pos and lapplies over
  # the pieces and computes the X_l matrices for each locus, and then just
  # returns this
  split(D2, factor(D2$gname, levels = unique(D2$gname))) %>%
    lapply(., function(x) {

      # a list for output
      loc_ret <- list()

      # for the dimnames of the X_l matrices
      geno_names <- expand.grid(gp = x$Allele, g = x$Allele) %>%
        dplyr::filter(as.integer(factor(g, levels = x$Allele)) <= as.integer(factor(gp, levels = x$Allele))) %>%
       dplyr::mutate(name = paste(g, gp, sep = allele_separator)) %>%
        dplyr::select(name) %>%
        unlist %>%
        unname

      # store the allele frequencies with the allele names
      loc_ret$freqs <- x$Freq
      names(loc_ret$freqs) <- x$Allele

      # while we are at it, let's also just get the genotype freqs, straight up,
      # under the assumption of HWE
      f <- x$Freq
      o <- outer(f, f) # this gets the freq of all ordered genotypes
      o[lower.tri(o)] <- 2 * o[lower.tri(o)] # factor of 2 for heterozygotes
      geno_freqs <- o[lower.tri(o, diag = TRUE)]
      names(geno_freqs) <- geno_names
      loc_ret$geno_freqs <- geno_freqs

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
#' does the microhaplotype error model with default parameters.  This creates and inserts the
#' C_l_true matrix and also the C_l matrix.  C_l_true is the "true" genotyping error model results
#' and C_l is what gets applied in the likelihood calculations.
#' @param XL a list of the loci like that created using \code{\link{long_markers_to_X_l_list}}.
#' The key thing that each list component needs is the named vector freqs of the allele frequencies.
#' The functions that compute genotyping error use the names of the allele to compute the
#' probabilities of the observed genotype given the true genotype.
#' @param ge_mod_assumed The genotyping error model assumed for the analysis.
#' @param ge_mod_true The actual, "true" genotyping error model for the simulation.
#' @param ge_mod_assumed_pars_list a named list of extra arguments (besides L)
#' to be passed to the ge_mod_assumed function.  Set it to NULL to use the default
#' values for the genotyping error model.
#' @param ge_mod_true_pars_list a named list of extra arguments (besides L)
#' to be passed to the ge_mod_true function. Set to NULL to use defaults.
#' @export
#' @examples
#' example(long_markers_to_X_l_list, package = "CKMRsim")
#' mh_cl_example <- insert_C_l_matrices(
#'     mh_example,
#'     ge_mod_assumed = ge_model_microhap1,
#'     ge_mod_true = ge_model_microhap1
#' )
insert_C_l_matrices <- function(
  XL,
  ge_mod_assumed,
  ge_mod_true,
  ge_mod_assumed_pars_list = NULL,
  ge_mod_true_pars_list = NULL
  ) {



  lapply(XL, function(x) {

    # take care of the extra parameters for the assumed genotyping model
    if(is.null(ge_mod_assumed_pars_list)) {
      ass_list <- list(L = x)
    } else {
      ass_list <- c(list(L = x), ge_mod_assumed_pars_list)
    }

    # take care of the extra parameters for the true genotyping model
    if(is.null(ge_mod_true_pars_list)) {
      true_list <- list(L = x)
    } else {
      true_list <- c(list(L = x), ge_mod_true_pars_list)
    }



    x$C_l_true <- do.call(ge_mod_true, true_list)
    x$C_l <- do.call(ge_mod_assumed, ass_list)
    x
  })
}




#' compute the Y_l matrices for each locus in a list
#'
#' Once you have the X_l matrices for each kappa, and the C_l matrices
#' have been inserted into the list, then this function cycles over the loci
#' and the kappa values and does the matrix multiplication to give you the
#' matrix Y_l.  This actually gives you the matrices Y_l_true and Y_l,
#' but for now they are the same.  I just want to put that into the code at this
#' point so I can work with them appropriately.
#' @param L the list of X_l and C_l matrices
#' @export
#' @examples
#' example(insert_C_l_matrices)
#' mh_yl_example <- insert_Y_l_matrices(mh_cl_example)
insert_Y_l_matrices <- function(L) {
  lapply(L, function(a) {
    a$Y_l_true <- lapply(a$X_l, function(b) {
      tmp <- t(b %*% a$C_l_true) %*% t(a$C_l_true)
      names(dimnames(tmp)) <- c("obs_geno_indiv_1", "obs_geno_indiv_2")
      tmp
    })
    a$Y_l <- lapply(a$X_l, function(b) {
      tmp <- t(b %*% a$C_l) %*% t(a$C_l)
      names(dimnames(tmp)) <- c("obs_geno_indiv_1", "obs_geno_indiv_2")
      tmp
    })
    a
  })
}
