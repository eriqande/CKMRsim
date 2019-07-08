

#' Extract log-likelihood ratios and associated values from a Qij object
#'
#' This function applies the prior weights on the pairs for the numerator (for generality)
#' and denominator and returns the log-likelihood ratios in a tidy format
#' @param Q a Qij object.  Like that returned by \code{\link{simulate_Qij}}
#' @param numer a named vector of weights to be given to the relationships for the numerator.
#' The names are the relationships and the values are the weights.  Relationships that
#' appear in the "Tos" field of the Q that do not appear will get weights of 0.
#' @param denom a named vector of weights to be given to the relationships for the denominator.
#' The names are the relationships and the values are the weights.  Relationships that
#' appear in the "Tos" field of the Q that do not appear will get weights of 0.
#' @export
extract_logls <- function(Q, numer, denom) {

  if(FALSE) {  # a block of stuff for testing
    numer = c(PO=1)
    denom = c(U=.999, FS=.001)
  }
  #### first some error checking
  # make sure numer does not ask for relationships that don't exist
  numer_bads <- setdiff(names(numer), names(Q[[1]]))
  if(length(numer_bads) > 0) {
    stop("Requested numer weights for ", paste(numer_bads, collapse = ", "), " but these do not exist as \"Tos\" in Q")
  }
  # same for denom
  denom_bads <- setdiff(names(denom), names(Q[[1]]))
  if(length(denom_bads) > 0) {
    stop("Requested numer weights for ", paste(denom_bads, collapse = ", "), " but these do not exist as \"Tos\" in Q")
  }

  # make sure there are no negative values in numer and denom
  stopifnot(all(numer >= 0))
  stopifnot(all(denom >= 0))

  # make sure that at least one relationship in Tos has a positive weight
  numer_pos <- numer[numer > 0]
  if(length(intersect(names(numer_pos), names(Q[[1]]))) == 0) {
    stop("Must have at least one positive weight in numer")
  }
  # same for denom
  denom_pos <- denom[denom > 0]
  if(length(intersect(names(denom_pos), names(Q[[1]]))) == 0) {
    stop("Must have at least one positive weight in denom")
  }

  #### Then normalize the weights so they are sure to sum to one
  numer_norm <- numer/sum(numer)
  denom_norm <- denom/sum(denom)

  # explicitly remove 0's, since they don't serve any purpose
  numer_norm_p_unsrt <- numer_norm[numer_norm > 0]
  denom_norm_p_unsrt <- denom_norm[denom_norm > 0]

  # then sort them from most to least weight
  numer_norm_p <- sort(numer_norm_p_unsrt, decreasing = TRUE)
  denom_norm_p <- sort(denom_norm_p_unsrt, decreasing = TRUE)



  # and get strings to summarise these weights
  numer_string <- paste(names(numer_norm_p), numer_norm_p, sep = "=", collapse = ":")
  denom_string <- paste(names(denom_norm_p), denom_norm_p, sep = "=", collapse = ":")

  # collect other information in a convenient place:
  simtype = attributes(Q)$simtype
  PO_sim  = attributes(Q)$PO_sim

  # now we lapply over the Trues, and for each one we compute the weighted logls
  results_df <- lapply(names(Q), function(true_relat) {
    x <- Q[[true_relat]]
    numerx <- x[names(numer_norm_p)]
    numerx_mat <- do.call(rbind, numerx)

    denomx <- x[names(denom_norm_p)]
    denomx_mat <- do.call(rbind,denomx)

    nconst <- max(numerx_mat)
    numer_exp <- exp(numerx_mat - nconst)
    numer_logl <- loggy_reweight(numerx_mat, numer_norm_p)
    denom_logl <- loggy_reweight(denomx_mat, denom_norm_p)
    logls <- numer_logl - denom_logl

    # now, some care in retrieving the log probability of being simulated from the true relationship
    if(simtype == "linked" | !(true_relat %in% names(x))) {
      true_log_prob <- NA
    } else {
      true_log_prob <- x[[true_relat]]
    }
    dplyr::data_frame(simtype = simtype,
                      PO_sim = PO_sim,
                      rando_miss_n = attributes(Q)$rando_miss_n,
                      numer_wts = numer_string,
                      denom_wts = denom_string,
                      true_relat = true_relat,
                      rep = 1:length(logls),
                      numer_logl = numer_logl,
                      denom_logl = denom_logl,
                      logl_ratio = logls,
                      true_log_prob = true_log_prob

    )
  }) %>%
    dplyr::bind_rows()

  results_df

}



#' reweight logs without underflow
#'
#' This is not exposed to the user
#' @param mat a matrix with r rows and c columns of logged values (like probabilities)
#' @param wts a vector of length r that is parallel to the rows in mat.
#' @details This is made for the situation where each column of mat holds log probabilities of a
#' given realization under r different hypotheses.   We want to get the log probability of the
#' given realization under a prior wts on the different hypotheses.  This is super easy when mat
#' consists of probabilities, but is a little trickier when mat holds logged prob values which need
#' to be exponentiated, summed with weights, and then re-logged.  One has to worry about underflow
#' in these situations.  This function deals with that by pulling a constant out of each and putting it
#' back in as a log for each column.
#' @export
#' @examples
#' # this is just a little thing to show that it gets the correct result:
#' tmp <- matrix(runif(120,0,1), nrow = 3)
#' wtt <- c(.23, .67, 1.8)
#' wts <- wtt/sum(wtt)
#' correct <- log(colSums(tmp * wts))
#'
#' # and then compute it with loggy_reweight
#' lw <- loggy_reweight(log(tmp), wts)
#'
#' # see they are give the same thing
#' all.equal(lw, correct)
#'
#' # see that they are only different at machine precision
#' lw - correct
loggy_reweight <- function(mat, wts) {
  logconsts <- colMeans(mat)  # use the mean to put the constant in the middle of everything
  mat2 <- mat - rep(logconsts, each = nrow(mat))
  log(colSums(exp(mat2) * wts)) + logconsts
}
