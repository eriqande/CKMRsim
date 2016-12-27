


#' implements the generalized Allele-Based genotyping error model
#'
#' This is the version into which you pass in your own matrices of allele
#' call probabilities W (corresponding to omega) and Ws (omega-star)
#' and vectors of allele-specific dropout rates. Notation follows the paper
#' for the most part. So see it for details.  Note that this is vectorized
#' over g, gp, h, and hp.  Note that g, gp, h, and hp can be character vectors
#' giving the allele names, so long as the vector D has names which are the allele
#' names and W and Ws have rownames and colnames which are the allele names.
#' @param g allelic type of gene copy g in the true genotype.
#' @param gp allelic type of gene copy g' in the true genotype
#' @param h allelic type of one of the gene copies of the observed genotyped
#' @param hp allelic type of the other gene copy of the observed genotyped
#' @param D vector of allele-specific dropout rates (delta) in the paper. Can be named
#' with the names of the alleles.
#' @param W the matrix omega of allelic call probabilities. Can have rownames and colnames that
#' are the allele names.
#' @param Ws the matrix omega* of allelic call probabilities
#' @return This function returns a vector that is in the order of the genotypes
#' as they are listed according to  g, gp, h, and hp.
#' @examples
#' # let's fabricate a locus with 5 alleles named "a1" through "a5"
#' # and make W random, but with a much higher chance of observing the
#' # correct allele.
#' set.seed(5)
#' alles <- paste("a", 1:5, sep = "")
#'
#' # here is a data frame of all possible genotypes
#' pg <- expand.grid(g=alles, gp=alles, stringsAsFactors = FALSE) %>%
#'   dplyr::filter(g<=gp)
#'
#' # and then we replicate that to get all possible combinations of true
#' # genotypes to observed genotypes
#' input_df <- cbind(pg[rep(1:nrow(pg), nrow(pg)), ],
#'                   pg[rep(1:nrow(pg), each = nrow(pg)), ])
#' names(input_df)[3:4] <- c("h", "hp")
#' input_df <- dplyr::tbl_df(input_df)
#'
#'
#' # here we make the matrix W (just some random numbers)
#' W <- runif(length(alles)^2, 0.001, 0.02) %>%
#'   matrix(nrow = length(alles))
#' dimnames(W) <- list(alles, alles)
#' diag(W) <- diag(W) + 1  # make it most likely to call the right allele
#' W <- t(apply(W, 1, function(x) x/sum(x)))  # normalize them to be probablities
#'
#' # make Wstar to be similar but less likely to miscall
#' Ws <- W
#' diag(Ws) <- 2
#' Ws <- t(apply(Ws, 1, function(x) x/sum(x)))
#'
#' # finally make a random D
#' D <- runif(length(alles), min = 0.001, max = 0.05)
#' names(D) <- alles
#'
#' # then compute the probs
#' gprobs <- input_df %>%
#'  dplyr::mutate(prob = general_allele_based_geno_err_model(g, gp, h, hp, D, W, Ws))
#'
#' # now, confirm that, for each true genotype, when we sum over the
#' # probs of the observed genotypes, we get 1.
#' gprobs %>%
#'  dplyr::mutate(ggp = paste(g, gp, sep = "-")) %>%
#'   dplyr::group_by(ggp) %>%
#'   summarise(probSum = sum(prob))
#'
#' # Eureka! it all checks out!
general_allele_based_geno_err_model <- function(g, gp, h, hp, D, W, Ws) {
  gh = cbind(g, h)
  gph = cbind(gp, h)
  ghp = cbind(g, hp)
  gphp = cbind(gp, hp)


  (h == hp) * (D[g] * Ws[gph] + D[gp] * Ws[gh]) +
    (1 - D[g] - D[gp]) * (W[gh] * W[gphp] + (h != hp) * W[gph] * W[ghp])
}



#' a simple helper function
#'
#' This will ultimately not be exported. It gives the probability that
#' H2 is observed from a haplotype that is truly H1 given a super
#' simple model of independent sequencing errors. This is not
#' intended to be vectorized, at the moment.
#' @param H1 the DNA bases at the SNPs in the first haplotype
#' @param H2 the DNA bases at the SNPs in the second haplotype
#' @param e per snp error rates
#'
hap_obs_prob <- function(H1, H2, e) {
  h1s <- strsplit(H1, "")
  h2s <- strsplit(H2, "")
  unlist(lapply(seq_along(h1s), function(x) prod(e[h1s[[x]] != h2s[[x]]])))
}

#' create matrix C for probability of observed genotypes from microhaplotype data
#'
#' This is intended for the case where the genotypes in question are composed of alleles
#' that are actually the multi-SNP haplotypes obtained from next generation sequence data.
#' In other words, all the SNPs occur on a single read and the phase is known because they
#' are all together on the read.  It allows for a SNP-specific sequencing error rate.  The
#' haplotypes must be named as strings of A, C, G, or T, (though they could be strings of
#' any characters---the function isn't going to check that!) and for now  we assume that if the
#' SNPs are multiallelic then genotyping errors to either of the alternate alleles are
#' equally likely.  Currently assumes that genotyping errors are equally likely in either
#' direction at a SNP, too.
#' @param haps character vector of strings that denote the haplotypes at the locus.  For example
#' "CCAG", "CTAG", "GCAG", etc.  Note that these should be in the same order as they are given
#' in the allele frequency definitions (so that the ordering of genotypes made from them will
#' be correct).  Each element of haps must be a string of the same number
#' of characters.  haps cannot be a factor.
#' @param snp_err_rates Vector of rates at which sequencing errors are expected at each of the SNPs
#' that are in the haplotype.  This recycles if its length is less than the number of SNPs in the
#' haplotypes.
#' @param dropout_rates Haplotype-specific rates of allelic dropout.  Recycles if need be.
#' @examples
#' # five haplotypes in alphabetical order
#' haps <- c("AACC", "GACC", "GATA", "GTCC", "GTTC")
#'
#' # make the matrix C
#' C_mat <- microhaplotype_geno_err_matrix(haps)
#'
#' # look at the first part of it
#' C_mat[1:5, 1:5]
microhaplotype_geno_err_matrix <- function(haps, snp_err_rates = 0.005, dropout_rates = 0.005, scale_by_num_snps = FALSE) {
  hap_length <- unique(nchar(haps))
  if(length(hap_length) > 1) stop("haps has variable number of SNPs: ", paste(hap_length, collapse = ", "))
  if(length(haps) != length(unique(haps))) stop("duplicate haplotype in haps")

  ser <- rep(snp_err_rates, length.out = hap_length)
  dr <- rep(dropout_rates, length.out = hap_length)

  if(scale_by_num_snps == TRUE) {
    ser <- ser / length(haps)  # this let's us get roughly the same geno error rates for things with many
                               # snps as for few snps
  }

  names(haps) <- haps

  # first, create the matrix W and Ws
  W <- outer(haps, haps, hap_obs_prob, e = ser)
  tmp<-W
  diag(tmp) <- 0
  diag(W) <- 1 - rowSums(tmp)
  Ws <- W

  # then create the dropout rates vector
  D <- rep(dropout_rates, length(haps))

  # now make a data frame of allelic states in genotype order that we
  # would need to fill the matrix C.  That is, cycling over true
  # genotypes fastest...
  tmp <- expand.grid(gp = haps, g = haps) %>%
    dplyr::filter(as.integer(factor(g, levels = haps)) <= as.integer(factor(gp, levels = haps))) %>%
    dplyr::select(g, gp) %>%
   dplyr::mutate(h = g, hp = gp)

  gdf <- cbind(tmp[rep(1:nrow(tmp), nrow(tmp)), c("g", "gp")],
               tmp[rep(1:nrow(tmp), each = nrow(tmp)), c("h", "hp")])

  # now compute the obs probs given true genotypes of those
  probs <- gdf %>%
   dplyr::mutate(cprob = general_allele_based_geno_err_model(g, gp, h, hp, D, W, Ws))

  C_mat <- matrix(probs$cprob, nrow = sqrt(nrow(probs)))
  dimnames(C_mat) <- list(true_geno = paste(tmp$g, tmp$gp, sep = "-"),
                          obs_geno = paste(tmp$g, tmp$gp, sep = "-"))

  # OK, return that dude
  C_mat
}

