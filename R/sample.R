
# a variety of functions that all deal with sampling from the big lists of
# quantities that were simulated by, for example, simulate_and_calc_Q

# these functions are where the actual Monte Carlo and Importance Sampling
# takes place.


#### NON-EXPORTED FUNCTIONS #####
# non-exported helper function to format printing of mixing proportions of relationships
#
# If R is just a single name it just prints that name.  If it is a named vector of
# proportions, it prints them as they are
# @param R the relationship.
# @param colchars What characters should go between then relationships in the formatting
# @examples
# format_mixed_r(c(U=.9997, FS=.0001, PO=.00001, HS=.0002))
format_mixed_r <- function(R, colchars = ",") {
  if(length(R) == 1) {
    return(paste(R))
  }

  if(is.null(names(R))) stop("argument R must be a named vector if it has length > 1 in format_mixed_r")

  paste(paste(names(R), R, sep ="="), collapse = colchars)
}


## This is the importance sampling workhorse function...
# Q is for Qvals, nu is the relationship in the numerator of lambda,
# de is the relationship in the denominator of lambda which always be
# U in the current context, tr is the true
# relationship, which will typically be taken to be "U" in these
# contexts, and pstar is the relationship of the importance sampling
# distribution, which in this context will almost always by nu.
# FNRs are the false negative rates you want to investigate.
imp_samp <- function(Q, nu, de, tr, pstar, FNRs) {
  # get the importance weights and the corresponding lambdas when
  # the sample is from pstar
  iw <- dplyr::data_frame(lambda = Q[[pstar]][[nu]] - Q[[pstar]][[de]],
                   impwt = exp(Q[[pstar]][[tr]] - Q[[pstar]][[pstar]])) %>%
    dplyr::arrange(dplyr::desc(lambda)) %>%
    dplyr::mutate(FPR = cumsum(impwt))

  # and now we gotta get the lambdas for the true correct relationship
  trues <- Q[[pstar]][[nu]] - Q[[pstar]][[de]]

  # get the lambda values those correspond to
  cutoffs <- quantile(trues, probs = FNRs)

  # then get the FPRs for each of those
  tmp <- lapply(cutoffs, function(x) {
    sum(iw$impwt[iw$lambda >= x])
  }) %>%
    unlist() %>%
    unname()

    dplyr::data_frame(FNR = FNRs, FPR = tmp, Lambda_star = cutoffs) %>%
      dplyr::arrange(FNR)
}


#### EXPORTED FUNCTIONS #####
#' sample Q values to get and analyze a sample of Lambdas with simple (non-mixture) hypotheses
#'
#' Once you have gotten an output from simulate_and_calc_Q you can pass that
#' to this function along with instructions on what quantities to compute.
#' This version assumes that the the denominator of Lambda and the true relationship can be specified as a
#' a simple, single
#' relationship (typically, and by default, "U"), rather than a mixture of
#' possible relationships. For the latter see \code{\link{mc_sample_mixture}}
#'
#' The output is a long format data frame.
#' @param Q the list that is the output of simulate_and_calc_Q.
#' @param nu the name of the relationship that is in the \strong{nu}merator
#' of the likelihood ratio (Lambda) whose distribution you wish to learn about.
#' It is a string, for example "FS", or "PO", or "U".  The Q values for that
#' relationship must be included in parameter Q.  If this is a vector, then
#' all different values are used in combination with all the values of
#' \code{de}, \code{tr}, and, possibly, \code{pstar}. Corresponds to column "numerator" in the output
#' @param de the relationship that appears in the \strong{de}nominator of Lambda.
#' By default it is "U".  Corresponds to column "denominator" in the output. If it
#' is a vector, then all values are done iteratively in combination with other values as
#' described for \code{nu}.
#' @param tr the true relationship of the pairs. Default is "U". (i.e. you are going to
#' get samples of Lambda under their distribution when the true relationship is tr).
#' Operates over all values if a vector. Corresponds to column "true_relat" in the
#' output.
#' @param method the Monte Carlo method to use.  Either "IS" for importance sampling,
#' "vanilla" for vanilla Monte Carlo---regular Monte Carlo without importance sampling---or
#' or "both". The method that was used for any row of the output is reported in the
#' column "mc_method".
#' @param pstar the relationship used for the importance sampling distribution.
#' If set as NA and importance sampling (method == "IS" or "both") is used, then
#' the value of \code{nu} is used as need be.  If not NA, then this can be a vector
#' of relationships.  Each value will be used in all combinations of pstar, nu, de, and tr.
#' This is reported in column "pstar" in the output.
#' @param fnr the false negative rates at which to evaluate the false positive rates.
#' These are reported in column "fnr" in the output. These should all be between
#' 0 and 1.  By default fnr is c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001).
#' @return A long format data frame.  It will have a column of \code{tot_loci} that gives the total
#' number of loci.
mc_sample_simple <- function(Q,
                               nu,
                               de = "U",
                               tr = "U",
                               method = c("IS", "vanilla", "both")[1],
                               pstar = NA,
                               FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001)
) {

  #### here test that everything is OK and catch input errors  ####
#  stopifnot(length(nu) == 1, length(de) == 1, length(tr) == 1, length(pstar) == 1)
  stopifnot(is.character(nu) == TRUE,
            is.character(de) == TRUE,
            is.character(tr) == TRUE,
            is.na(pstar) || is.character(pstar) == TRUE)
  tr_lack <- setdiff(c(tr, pstar[!is.na(pstar)]), names(Q))
  if(length(tr_lack) > 0) stop("Asking for relationships in tr or pstar that are not available in Q: ",
                               paste(tr_lack, collapse = ", "))
  nu_lack <- setdiff(c(nu, de), unique(unlist(lapply(Q, names))))
  if(length(nu_lack) > 0) stop("Asking for relationships in nu or de that are not available in Q: ",
                               paste(nu_lack, collapse = ", "))
  stopifnot(all(FNRs > 0  & FNRs < 1) == TRUE)
  stopifnot(length(method)==1, method %in% c("IS", "vanilla", "both"))


  #### cycle over different relationships and do the calculations ####

  lapply(nu, function(nu_) {
    lapply(de, function(de_) {
      lapply(tr, function(tr_) {
        is <- NULL  # setting these to default NULL for easy row-binding later if they didn't get set
        van <- NULL
        if(method == "IS" || method == "both") {
          if(is.na(pstar)) {  # just taking care of defaulting pstar to nu_ if pstar is NA
            pstar_tmp <- nu_
          } else {
            pstar_tmp <- pstar
          }
          is <- lapply(pstar_tmp, function(pstar_) {
            tmp <- imp_samp(Q = Q, nu = nu_, de = de_, tr = tr_, pstar_, FNRs)
            tmp$pstar <- pstar_
            tmp
          }) %>%
            dplyr::bind_rows()

          is$mc_method = "IS"
        }
        if(method == "vanilla" || method == "both") {
          NULL; # gotta add in here code for the vanilla method
          # van$mc_method = "vanilla"
        }

        ret <- dplyr::bind_rows(is, van)
        ret$numerator = nu_
        ret$denominator = de_
        ret$true_relat = tr_
        ret
      }) %>%
        dplyr::bind_rows()
    }) %>%
      dplyr::bind_rows()
  }) %>%
    dplyr::bind_rows()
}
