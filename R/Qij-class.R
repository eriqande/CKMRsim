

#' Constructor function for the class Qij
#'
#' The Qij class holds the output of simulated genotype pairs whose
#' probabilities have been evaluated under different pairwise relationship
#' hypotheses.  It is output by the function \code{\link{simulate_Qij}}.
#'
#' A Qij object is quite simply a list of lists of vectors.  At the first
#' level the components are named by relationships from which data were
#' simulated (the "froms").  Within each of those is a named list of
#' relationships for which the probability of each simulated genotype
#' pair was evaluated (the "tos").
#' @param Q the list of lists to be turned into a Qij object.
#' @param unlinked Logical that says whether simulation was of unlinked markers or not
#' @param forceLinkagePO Logical.  If linked markers, should PO be forced to be simulated via Mendel
#' @inheritParams simulate_Qij
Qij_class <- function(Q, unlinked, forceLinkagePO, miss_mask_mat, rando_miss_wts, rando_miss_n) {
  # first, make sure that the names of each component are the same
  names1 <- names(Q[[1]])

  # ECA commented out the following two lines when he went to the list-based calc_relats option
  #names_correct <- all(sapply(Q, function(x) all(names(x) == names1)))
  #if(!names_correct) stop("Q does not seem to have consistent \"tos\" relationships.")

  # then make sure that the number of reps in all cases is correct
  reps <- length(Q[[1]][[1]])
  reps_correct <- all(reps == unlist(sapply(Q, function(x) sapply(x, length))))
  if(!reps_correct) stop("Q does not seem to have consistent reps across all relationships.")

  class(Q) <- "Qij"
  if(unlinked == FALSE) {
    attr(Q, "simtype") <- "linked"
  } else {
    attr(Q, "simtype") <- "unlinked"
  }

  if(forceLinkagePO == TRUE) {
    attr(Q, "PO_sim") <- "forced_MENDEL"
  } else {
    attr(Q, "PO_sim") <- "not_forced_MENDEL"
  }

  attr(Q, "miss_mask_mat_NULL") <- is.null(miss_mask_mat)
  attr(Q, "rando_miss_wts") <- rando_miss_wts
  attr(Q, "rando_miss_n") <- rando_miss_n


  Q
}


#' A simple function to print out the calc_relats in a list-oriented fashion
#'
#' @param Q an object of class \code{\link{Qij_class}}
print_calc_relats <- function(Q) {
  ll <- lapply(names(Q), function(q) {
    paste0("    ", q, ":  ", paste(names(Q[[q]]), collapse = ", "))
  })
  paste(c("", ll), collapse = "\n")
}


#' format method for Qij class (to print)
#'
#' Just prints relevant information for a quick look.
#' @param x an object of class \code{\link{Qij_class}}.
#' @param ... additional arguments to format. (But nothing implemented).
#' @export
format.Qij <- function(x, ...) {
  ret <- character()
  ret[1] <- paste0("A Qij object with ", length(x[[1]][[1]]), " reps")
  if(attributes(x)$simtype == "unlinked") {
    ret[2] <- paste0("simulated with markers ", attributes(x)$simtype)
  } else {
    ret[2] <- paste0("simulated with markers ", attributes(x)$simtype, " with PO treated as ", attributes(x)$PO_sim)
  }
  ret[3] <- paste0("\"sim_relats\" relationships: ", paste(names(x), collapse = ", "))
  ret[4] <- paste0("\"calc_relats\"   relationships: ", print_calc_relats(x))
  ret[5] <- paste0("rando_miss_n: ", attributes(x)$rando_miss_n)
  ret
}


#' print a Qij object nicely
#'
#' For nice printing.
#' @param x an object of class \code{\link{Qij_class}}
#' @param ... additional arguments to print. (But nothing implemented).
#' @export
print.Qij <- function(x, ...) {
  cat(format(x), sep = "\n")
}
