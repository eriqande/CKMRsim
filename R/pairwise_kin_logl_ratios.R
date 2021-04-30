

#' compute pairwise log likelihood ratios between two sets of individuals
#'
#' This computes log likelihoods for all pairs formed by comparing every individual
#' in D1 to every individual in D2.  If D1 and D2 are different objects, it is assumed
#' that they include disjoint groups of individuals, and all pairs are returned. If
#' D1 and D2 are the same object, then the output is filtered so that each
#' pairwise comparison appears only once, and so that the comparisons between the
#' same individual are dropped.
#' @param D1 long format genotype data frame of the first set of individuals.  If doing
#' parent-offspring, this should be the parents.
#' @param D2 long format genotype data frame of the second set of individuals. If doing
#' parent-offspring, this should be the offspring.
#' @param CK the ckmr object created previously from the data that D1 and D2 came from.
#' @param numer the relationship that should be the numerator of the logl ratio.
#' This relationship must appear in the ckmr object CK.
#' @param denom the relationship that should be the denominator of the logl ratio.
#' This relationship must appear in the ckmr object CK. Leave empty to only calculate the
#' log likelihood of of the numerator.
#' @param keep_top If this is given as an integer, then for each individual in D2, only the keep_top
#' individuals with the highest logl ratios are retained for the result.  The rest
#' are discarded.
#' @param num_cores The number of cores to use with mclapply.  This is ignored
#' if using windows (in which case it is always set to one.)  Otherwise, by
#' default it uses all the cores available.
#' @export
pairwise_kin_logl_ratios <- function(D1, D2, CK, numer, denom = NULL, keep_top = NULL,
                                     num_cores = parallel::detectCores()) {

  # first, make two genotype matrices
  M1 <- create_integer_genotype_matrix(D1, CK$orig_data)
  M2 <- create_integer_genotype_matrix(D2, CK$orig_data)

  # now, flatten the ckmr object into a useful flat form
  numer_flat <- flatten_ckmr(CK, numer)

  # then compute the logl ratios
  logl_flat <- numer_flat  # this just copies over some information

  if(! is.null(denom)) {
    denom_flat <- flatten_ckmr(CK, denom)
    logl_flat$probs <- log(numer_flat$probs / denom_flat$probs)
    rm(denom_flat)

  } else {
    logl_flat$probs <- log(numer_flat$probs)
  }

  # now, logl_flat is what we want to use
  rm(numer_flat)

  # now figure out how many cores to run this on
  if(.Platform$OS.type == "windows") {
    num_cores <- 1
  }

  idx <- 1:nrow(M2)
  names(idx) <- idx
  comps <- parallel::mclapply(idx,
                              mc.cores = num_cores,
                              FUN = function(i) {  # note, this can be parallelized with mclapply
    tmp <- comp_ind_pairwise(S = M1, T = M2, t = i, values = logl_flat$probs, nGenos = logl_flat$nGenos, Starts = logl_flat$base0_locus_starts)

    if(!is.null(keep_top)) {
      ret <- tmp[rev(top_index(tmp$value, min(keep_top, nrow(tmp)))), ]  # just take the top 5 from each
    } else {
      ret <- tmp
    }
    ret
  }) %>%
    dplyr::bind_rows(.id = "D2_indiv") %>%
    tibble::as_tibble() %>%
    dplyr::rename(D1_indiv = ind) %>%
    dplyr::mutate(D2_indiv = as.integer(D2_indiv))

  # now, if D1 and D2 are the same, give a message and
  # keep only one instance of each pair, and drop the
  # self-comparisons.
  if(identical(D1, D2)) {
    message("D1 and D2 are identical: dropping self comparisons and keeping only first instance of each pair")
    comps <- comps %>%
      dplyr::filter(D2_indiv < D1_indiv)
  }

  # now, put the IDs back in there and return the thing
  comps %>%
    dplyr::mutate(D2_indiv = rownames(M2)[D2_indiv],
           D1_indiv = rownames(M1)[D1_indiv]) %>%
    dplyr::rename(logl_ratio = value)
}
