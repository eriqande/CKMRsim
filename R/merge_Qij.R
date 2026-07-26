

#' Merge the results from multiple Qij objects
#'
#' This is a simple function that you will likely only ever need to use when
#' you have split up the simulation burden to do many reps of the same scenarios,
#' typically when you are parallelizing over different pairs of individuals
#' with the same amount of missing data.  It is up to the user to make sure that
#' this is used intelligently. i.e., don't merge Qij objects that were simulated
#' in different fashions (i.e. linked and unlinked).
#' @param Qlist A list of Qij objects
#' @return A Qij object with component names identical to that of the first Qij object
#' on the list, with values that are the concatenated versions of all of the
#' Qij values in the list. The attributes of the return item are the same as the
#' first Qij object on the list.
#' @export
merge_Qij <- function(Qlist) {
  Qnames <- names(Qlist[[1]])
  ret <- lapply(Qnames, function(qn) {  # cycle over the sim_relats
    innerQnames <- names(Qlist[[1]][[qn]])
    inner <- lapply(innerQnames, function(iqn) {
      lapply(Qlist, function(q) {
        q[[qn]][[iqn]]
      }) %>% unlist()
    })
    names(inner) <- innerQnames
    inner
  })
  names(ret) <- Qnames
  attributes(ret) <- attributes(Qlist[[1]])
  ret
}
