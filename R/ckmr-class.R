

#' constructor function for the ckmr class
#'
#' This is a constructor function and where we will document the \code{ckmr} class.
#' The \code{ckmr} class is an object that includes all the stuff necessary to
#' sample genotypes of pairs of individuals at different loci.  Primarily this
#' will be the object returned by \code{\link{ckmr_create}}. A \code{ckmr} object
#' is a named list.  At the first level are components \code{orig_data} (the original
#' data frame of markers that went into its construction) and \code{loci}, a component in which
#' each component corresponds to a locus and has the following
#' components
#' \describe{
#'   \item{freqs}{A vector of allele frequencies, named by the alleles.}
#'   \item{X_l}{A list, named by the pairwise relationships investigated, of matrices
#'   giving the probabilities of pairwise genotypes between two individauls having that
#'   relationship. }
#'   \item{C_l}{A matrix of probabilities giving the conditional probability of observing
#'   each possible genotype, given the true underyling state of the genotype in an individual.
#'   This is what encodes the genotyping error model that is \emph{assumed} for computing the Q ratios
#'   (See the paper).}
#'   \item{C_l_true}{A matrix of probabilities giving the conditional probability of observing
#'   each possible genotype, given the true underyling state of the genotype in an individual.
#'   This is what encodes the genotyping error model that is \emph{actually the true model} and which
#'   is used for simulating the genotypes when computing the Q ratios
#'   (See the paper).}
#'   \item{Y_l}{The assumed probability of observing pairwise genotypes between two individuals.  This is a
#'   list of matrices, one for each relationship, as named.  The quantities are obtained by
#'   matrix multiplication of \code{X_l} to \code{C_l}.}
#'   \item{Y_l_true}{The true probability of observing pairwise genotypes between two individuals.  This is a
#'   list of matrices, one for each relationship, as named.  The quantities are obtained by
#'   matrix multiplication of \code{X_l} to \code{C_l_true}.}
#' }
#' @param L the list that you want to turn into a \code{ckmr} object.
ckmr_class <- function(L) {
  stopifnot(names(L) == c("orig_data", "loci"))

  # gonna just do minimal checking here
  check_there <- all(sapply(L$loci, function(x) base::setequal(names(x), c("freqs", "X_l", "C_l_true", "C_l", "Y_l_true", "Y_l") )))
  if(!check_there) stop("L doesn't seem to have the right components to be a ckmr object.")

  class(L) <- "ckmr"
  L
}





#' format method for ckmr class (to print)
#'
#' Just prints relevant information for a quick look.
#' @param C an object of class \code{\link{ckmr_class}}.
format.ckmr <- function(C) {

  CL <- C$loci

  # get range of allele nums
  arange <- range(sapply(CL, function(x) length(x$freqs)))

  relnames <- names(CL[[1]]$X_l)

  # check to see if assumed genotyping error is same as true error.  We just check it on the
  # first relationship
  gerr_same <- all(sapply(CL, function(x) all(x$Y_l[[1]] == x$Y_l_true[[1]])))
  if(gerr_same) {
    sameness <- "the same"
  } else {
    sameness <- "different"
  }

  ret <- character()
  ret[1] <- paste0("A ckmr object with ", length(CL), " loci having between ", arange[1], " and ", arange[2], " alleles.")
  ret[2] <- paste0("Locus names: ", names(CL)[1], ", ", names(CL)[2], ",   ...   , ", names(CL[length(CL)]))
  ret[3] <- paste0("Relationships:  ", paste(relnames, collapse = ", "))
  ret[4] <- paste0("Assumed and true genotyping error models are ", sameness)
  ret
}


#' print method for ckmr class
#'
#' Just wraps a call to the format.ckmr function
#' @param C an object of class \code{\link{ckmr_class}}.
print.ckmr <- function(C) {
  cat(format(C), sep = "\n")
}


#' return a vector of the relationships in a ckmr_class object
#'
#' Just a simple convenience function
#' @param CK a ckmr object
#' @export
ckmr_relats <- function(CK) {
  stopifnot(class(CK) == "ckmr")
  names(CK$loci[[1]]$X_l)
}
