

#' break up (or augment) LTpairs into chunks of size reps
#'
#' This is a weird little function that lets us break up the LTpairs
#' that will be used for the miss_mask_mats in the Snakemake missing-data
#' simulation workflow.
#' The idea is that we can get extra parellelization and keep the memory overhead
#' for each job down low by restricting the number of reps done.  If the names
#' of the LTPairs list are like c("120", "121",...), the the new names will
#' be like c("120_1", "120_2", "120_3", 121_1", "121_2", etc...)
#' @param L a list with elements that are the LT_pairs.  Each element is a list
#' with elements r and c.
#' @param R the number of reps desired.  Should typically be about 50,000 to make
#' the memory footprint of each simulation manageable.
break_up_LTpairs <- function(L, R) {
  nested_list <- lapply(L, function(x) {
    l <- length(x$r)
    if(l < R) {
      i <- rep(1:l, length.out = R)
      ret <- list()
      ret[["1"]]$r <- x$r[i]
      ret[["1"]]$c <- x$c[i]
      return(ret)
    } else {
      num_splits <- ceiling(l / R)
      spl <- rep(1:num_splits, each = R, length.out = l)
      rs <- split(x$r, spl)
      cs <- split(x$c, spl)
      ret <- list()
      for(i in 1:length(rs)) {
        ret[[as.character(i)]]$r <- rs[[i]]
        ret[[as.character(i)]]$c <- cs[[i]]

        # if the last one has less than R pairs in it, fill the remainder with
        # pars drawn randomly from all of the pairs
        ll <- length(ret[[as.character(i)]]$r)
        if(ll < R) {
          rem <- R - ll
          idxs <- sample(1:l, size = rem)
          ret[[as.character(i)]]$r <- c(ret[[as.character(i)]]$r, x$r[idxs])
          ret[[as.character(i)]]$c <- c(ret[[as.character(i)]]$c, x$c[idxs])
        }
      }
      return(ret)
    }
  })

  # I used to Flatten and rename these, but now that I am writing them out for
  # snakemake, I do not do that.  We can use the nested list nicely.
  #unnested_list <- setNames(
  #  unlist(nested_list, recursive = FALSE),
  #  nm = unlist(lapply(names(nested_list), function(outer_name) {
  #    inner_names <- names(nested_list[[outer_name]])
  #    paste0(outer_name, "_", inner_names)
  #  }))
  #)

  return(nested_list)
}
