

#' Constructor function for the class Qij
#'
#' The Qij class holds the output of simulated genotype pairs whose
#' probabilities have been evaluated under different pairwise relationship
#' hypotheses.  It is output by the function \code{\link{simulate_Qij}}.
#'
#' A Qij object is quite simply a list of lists of vectors.  At the first
#' level the components are named by relationships from which data were
#' simualated (the "froms").  Within each of those is a named list of
#' relationships for which the probability of each simulated genotype
#' pair was evaluated (the "tos").
#' @param Q the list of lists to be turned into a Qij object.
Qij_class <- function(Q) {
  # first, make sure that the names of each component are the same
  names1 <- names(Q[[1]])
  names_correct <- all(sapply(Q, function(x) all(names(x) == names1)))
  if(!names_correct) stop("Q does not seem to have consistent \"tos\" relationships.")

  # then make sure that the number of reps in all cases is correct
  reps <- length(Q[[1]][[1]])
  reps_correct <- all(reps == sapply(Q, function(x) sapply(x, length)))
  if(!reps_correct) stop("Q does not seem to have consistent reps across all relationships.")

  class(Q) <- "Qij"
  Q
}



#' format method for Qij class (to print)
#'
#' Just prints relevant information for a quick look.
#' @param Q an object of class \code{\link{Qij_class}}.
format.Qij <- function(Q) {
  ret <- character()
  ret[1] <- paste0("A Qij object with ", length(Q[[1]][[1]]), " reps")
  ret[2] <- paste0("\"Froms\" relationships: ", paste(names(Q), collapse = ", "))
  ret[3] <- paste0("\"Tos\"   relationships: ", paste(names(Q[[1]]), collapse = ", "))
  ret
}


print.Qij <- function(Q) {
  cat(format(Q), sep = "\n")
}
