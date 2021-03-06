#' create matrix C for probability of observed genotypes from microhaplotype data
#'
#' This is intended for the case where the genotypes in question are composed of alleles
#' that are actually the multi-SNP haplotypes obtained from next generation sequence data.
#' In other words, all the SNPs occur on a single read and the phase is known because they
#' are all together on the read.  It allows for a per-locus or a per-SNP-specific sequencing error rate.  The
#' haplotypes must be named as strings of A, C, G, or T, (though they could be strings of
#' any characters---the function isn't going to check that!) and for now  we assume that if the
#' SNPs are multiallelic then genotyping errors to either of the alternate alleles are
#' equally likely.  Currently assumes that genotyping errors are equally likely in either
#' direction at a SNP, too.
#' @param L an element of the list created by \code{\link{long_markers_to_X_l_list}}. Such
#' an element basically holds the information at a single locus.  The idea here
#' is that every ge_mod_* function takes in an object like L, and then can
#' use any piece of information in it about alleles or genotypes to
#' configure a genotyping error model.
#' @param miscall_rate The rate at which microhaplotype alleles are mis-called.  If
#' option \code{per_snp_rate} is TRUE, then this is the rate at which each SNPs might
#' get miscalled, such that the pverall miscall rate for microhaplotypes with more SNPs
#' is higher than for microhaplotypes with few SNPs.
#' @param dropout_rate Rate of allelic dropout.
#' @param per_snp_rate Logical.  If true, then the overall mis-call rate for the
#' microhaplotype locus is \code{miscall_rate} times the number of SNPs in the microhapotype.
#' The default is FALSE, and it this option is not really recommended.
#' @export
#' @examples
#' # here is some example information about a microhap:
#' example_L_microhap
#'
#' # now we can feed it in to the function with default parameter values.
ge_model_microhap1 <- function(L, miscall_rate = 0.005, dropout_rate = 0.005, per_snp_rate = FALSE) {

  if(is.null(names(L$freqs))) stop("No names attribute on the allele frequencies")

  # A word about this "haps" variable (below).  It is a character vector
  # of strings that denote the haplotypes at the locus.  For example
  # "CCAG", "CTAG", "GCAG", etc.   Each element of haps must be a
  # string of the same number of characters.  haps cannot be a factor.
  haps <- names(L$freqs)
  hap_length <- unique(nchar(haps))
  if(length(hap_length) > 1) stop("haps has variable number of SNPs: ", paste(hap_length, collapse = ", "))
  if(length(haps) != length(unique(haps))) stop("duplicate haplotype in haps")

  # first, get a matrix of the number of mismatches between each of the
  # microhap alleles:
  mism <- outer(haps, haps, function(x, y) {
    lapply(1:length(x), function(i) {
      sum(strsplit(x[i], "")[[1]] != strsplit(y[i], "")[[1]])
    }) %>%
      unlist()
  })
  colnames(mism) <- haps
  rownames(mism) <- haps

  # then make the matrix W
  if(per_snp_rate == TRUE) {
    miscall_rate = miscall_rate * nchar(haps[1])
  }

  W <- miscall_rate * t(apply(mism, 1, function(x) x / sum(x)))
  diag(W) <- 1 - miscall_rate

  # then create the dropout rates vector
  D <- rep(dropout_rate, length.out = length(haps))
  names(D) <- haps  # NOTE! It is very important that these be named with the allelic types

  # and then combine allelic miscalls and dropouts:
  ret <- combine_miscalls_and_dropouts(geno_names = names(L$geno_freqs),
                                D = D,
                                W = W,
                                Ws = W)

  ret
}

