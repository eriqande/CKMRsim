
# a variety of functions that all deal with sampling from the big lists of
# quantities that were simulated by, for example, simulate_and_calc_Q

# these functions are where the actual Monte Carlo and Importance Sampling
# takes place.




#' sample Q values to get and analyze a sample of Lambdas
#'
#' Once you have gotten an output from simulate_and_calc_Q you can pass that
#' to this function along with instructions on what quantities to compute
#' This version assumes that the denominator of Lambda is a simple single
#' relationship (typically, and by default, "U"), rather than a mixture of
#' possible relationships. For the latter see \code{\link{monte_carlo_mixed}}
#'
#' The output is a long format data frame.
#' @param Q the list that is the output of simulate_and_calc_Q.
#' @param nu the name of the relationship that is in the \strong{nu}merator
#' of the likelihood ratio (Lambda) whose distribution you wish to learn about.
#' It is a string, for example "FS", or "PO", or "U".  The Q values for that
#' relationship must be included in parameter Q.  If this is a vector, then
#' all different values are used. Corresponds to column "numerator" in the output
#' @param de the relationship that appears in the \strong{de}nominator of Lambda.
#' By default it is "U".  Corresponds to column "denominator" in the output. If it
#' is a vector, then all values are done iteratively.
#' @param tr the true relationship of the pairs. Default is "U". (i.e. you are going to
#' get samples of Lambda under their distribution when the true relationship is tr).
#' Operates over all values if a vector. Corresponds to column "true_relat" in the
#' output.
#' @param method the Monte Carlo method to use.  Either "IS" for importance sampling,
#' "vanilla" for vanilla Monte Carlo---regular Monte Carlo without importance sampling---or
#' or "both". The method that was used for any row of the output is reported in the
#' column "mc_method".
#' @param pstar the relationship used for the importance sampling distribution.  By
#' default this is the same relationship as \code{nu}. This is reported in column
#' "pstar" in the output.
#' @param fnr the false negative rates at which to evaluate the false positive rates.
#' These are reported in column "fnr" in the output. These should all be between
#' 0 and 1.  By default fnr is c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001).
#' @param rando_miss_n how many loci to mask at random on each
#' iteration.  If this is a vector than it does the whole analysis for each
#' one of the values.  Corresponds to column "rando_miss_n" in the output.
#' Default is 0.
#' @param rando_miss_wts weights to be given to different loci that influence whether they
#' will be one of the rando_miss_n missing loci in any iteration.  These will be recycled
#' (or truncated) to have length equal to the number of loci (or to the first_n_loci if that
#' is in effect), and they will be normalized to sum to one as appropriate (so you can provide
#' them in unnormalized form.)  The idea of this is to be able to use observed rates of
#' missingness amongst loci to mask some loci as missing.  Given as a comma-delimited
#' string in column "rando_miss_wts" in the output.
#' @param first_n_loci  a vector telling how many loci to use.  For each value, n, it
#' will use just the first n loci.  This is useful if you want to duplicate your loci and see
#' how many more loci typical of what you have already collected would be required to achieve
#' a certain accuracy for relationship inference. Corresponds to column "first_n_loci" in the
#' output.  By default it's NA, in which case it has no effect. Incompatible with the option
#' @return A long format data frame.  It will have a column of \code{tot_loci} that gives the total
#' number of loci.
monte_carlo_single <- function(Q,
                               nu,
                               de = "U",
                               tr = "U",
                               method = c("IS", "vanilla", "both")[1],
                               pstar = nu,
                               FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001),
                               rando_miss_n = 0,
                               rando_miss_wts = 1,
                               first_n_loci = NA
) {

  # here test that everything is OK...
  NULL

}
