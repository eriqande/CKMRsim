#### Some Helper Functions First ##########

#' Collect Lower Triangle Pairs by Value
#'
#' Given a symmetric square integer matrix `M`, this function extracts all (r, c)
#' index pairs from the lower triangle (r > c) for each unique value `V*` present
#' in that region of the matrix. The result is a named list where each name corresponds
#' to a value `V*`, and each element is a list containing two integer vectors: r and c.
#'
#' @param M A square numeric matrix. Assumed to be symmetric with integer-like values.
#'
#' @return A named list. Each element corresponds to a unique value in the lower triangle of `M`,
#'         and contains a list with components r and c, the row and column indices (with r > c)
#'         where that value occurs.
#'
#' @examples
#' set.seed(1)
#' M <- matrix(sample(80:180, size = 100^2, replace = TRUE), nrow = 100)
#' result <- collect_lower_triangle_pairs(M)
#' str(result[["100"]])
#'
collect_lower_triangle_pairs <- function(M) {
  stopifnot(is.matrix(M))
  stopifnot(nrow(M) == ncol(M))

  n <- nrow(M)

  # Use which() on the lower triangle to get indices
  lower_indices <- which(lower.tri(M), arr.ind = TRUE)
  values <- M[lower_indices]

  # Split the row and column indices by the value in the matrix
  split_indices <- split(as.data.frame(lower_indices), values)

  # Format the result as requested: list named by V*, each element is a list with $r and $c
  result <- lapply(split_indices, function(df) list(r = df$row, c = df$col))

  return(result)
}



#' make a miss_mask_mat for simulate_Qij()
#'
#' This uses the missing data patterns at pairs of individuals
#' to specify the missing data patterns.
#' @param MG a matrix of N rows and L columns with a 0 if the data are missing
#' and a 1 otherwise.
#' @param L a list with two elements, r and c, each one a vector of the same length.
#' Each pair of r,c values gives a pair of individuals (rows in MG) whose missing
#' data patterns are to be used.
#' @param reps Minimum number of Monte Carlo reps to do.
#' @return A miss_mask_mat: A logical matrix with L columns and R rows. The (r,c)-th
#' element is TRUE if the c-th locus should be considered missing in the r-th
#' simulated sample. Note that the number of rows is the number of Monte Carlo
#' reps that should be done.
make_miss_mask_mat <- function(MG, L, reps) {
  NP <- length(L$r)
  Recyc <- 1  # assume only only copy of the patterns
  if(NP < reps) {
    Recyc <- ceiling(reps / NP)
  }
  Mat <- lapply(1:NP, function(i) {
    !as.logical(MG[L$r[i],] *  MG[L$c[i],])
  }) %>%
    do.call(what = rbind, args = .)

  # now, recycle rows as necessary
  recy <- rep(1:NP, times = Recyc)
  ret <- Mat[recy,]
  ret
}


######## END OF HELPER FUNCTIONS #################

#' Simulate Qijs across multiple relevant levels of missing data
#'
#' Some data sets have a lot of missing markers.  If this is the case, it is not
#' OK to just do the simulations as if there is no missing data.  This function
#' wraps up a lot of different steps that can be taken to try to get more accurate
#' "first-pass" FPRs and FNRs for situations with a lot of missing data.
#' The steps are:
#' - Tabulate the distribution of the number of informative (i.e., not missing in
#'   either member of the pair) markers, across all pairs. (Note, this requires
#'   that you have an actual data set that you are trying to do relationship
#'   inference in.)
#' - Estimate missingness rates per locus, and from that calculate the rate of
#'   missingness in pairs, under a simple independence assumption.
#' - Simulate Q_ij values at a series of different numbers of non-missing loci
#'   to calculate FPRs and FNRs for those.
#'
#' Note that this function requires a CKMR object `C` that is suitable for both
#' linked an unlinked simulation.  So, if you don't have that, then you best
#' prepare it, even if it means sprinkling your markers into a pseudogenome
#' as described [here](https://eriqande.github.io/tws-ckmr-2022/kin-finding-lab.html#power-for-kin-finding-while-accounting-for-physical-linkage).
#' @return This function returns a list.  More on that later.
#' @inheritParams create_integer_genotype_matrix
#' @inheritParams simulate_Qij
#' @param num_points the number of different values between the lowest observed
#' number of pairwise non-missing genotypes and the highest, inclusive, that
#' simulations will be performed for.
#' @param num_cores Number of cores to parallelize the simulations over (using mclapply)
#' from the parallel package.  On Windows, parallelization is not available from forking
#' so this must remain equal to 1 on Windows.
#' @param tabulate_and_exit if TRUE, then the function exits after tabulating the
#' number of pairs with different numbers on shared non-missing genotypes, and returns
#' only information about that.  This is useful for getting an overview of what the
#' missing data structure in your data look like.
#' @param simulation_approach this determines how simulations should be done. The
#' options are:
#' - `simple_miss_freq`: do simulations at `num_points` different values of the number
#' of pairwise shared genotypes, drawing the pattern of missingness from a very simple
#' independent random sampling model according to the frequency of missing data in all
#' loci.
#' - `pairwise_patterns`: do a separate simulation for each number of shared non-missing
#' loci observed amongst all pairs.  The simulation will do at least one rep for each
#' pair within a category of the number of shared non-missing loci, and it will use
#' the observed patterns of missing loci for all pairs. If the number Np of pairs in a category
#' is less than `reps`, then the actual number of reps done will be the smallest integer
#' multiple of Np >= `reps`, and the patterns of missing data will be recycled.
#' @inheritDotParams simulate_Qij sim_relats calc_relats reps
#'
#' @export
#' @examples
#' # this is just here for testing at the moment
#' LG <- read_rds("/tmp/LG.rds")
#' C <- read_rds("/tmp/C.rds")
simulate_missing_data_array <- function(
    LG,
    C,
    num_points = 11,
    num_cores = 1,
    tabulate_and_exit = FALSE,
    simulation_approach = "simple_miss_freq",
    ...
) {

  ret <- list()

  #### Step 1. Tabulate missing data quantities across pairs  ####
  AF <- C$orig_data
  MG <- create_integer_genotype_matrix(LG, AF)
  L <- ncol(MG)  # number of loci.  Good to have

  # MG is N rows by L columns.  We want to turn values >=0 to 1 and = -1 to 0
  MG[MG >= 0] <- 1
  MG[MG == -1] <- 0
  big <- MG %*% t(MG)

  pairwise_non_miss_counts <- tibble(
    num_non_miss = big[lower.tri(big)]
  ) %>%
    count(num_non_miss)

  pairwise_non_miss_counts_plot <- ggplot(pairwise_non_miss_counts, aes(x = num_non_miss, y = n)) +
    geom_col() +
    xlab("Number of loci missing in neither pair member") +
    ylab("Number of pairs")

  # while we are at it, let's count up the number of non-missing loci in each individual
  non_miss_counts_by_indiv <- tibble(
    num_non_miss = rowSums(MG)
  ) %>%
    count(num_non_miss)

  non_miss_counts_by_indiv_plot <- ggplot(non_miss_counts_by_indiv, aes(x = num_non_miss, y = n)) +
    geom_col() +
    xlab("Number of non-missing loci") +
    ylab("Number of individuals")


  #### Step 1.5*.  Collect the pairs that have the same number of shared loci (if needed)
  if(simulation_approach == "pairwise_patterns") {
    LTpairs <- collect_lower_triangle_pairs(big)
  }


  #### Step 2. Estimate missing data rate per locus   ####

  # this is simple
  miss_rates_by_locus <- colMeans(MG == 0)

  # the rate at which things are non-missing in pairs is
  pairwise_miss_rates_by_locus = 1 - (1 - miss_rates_by_locus)^2


  #### Stop for a moment and start filling up the return list ####
  ret$background$values$pairwise_non_miss_counts <- pairwise_non_miss_counts
  ret$background$values$non_miss_counts_by_indiv <- non_miss_counts_by_indiv
  ret$background$values$miss_rates_by_locus <- miss_rates_by_locus
  ret$background$values$pairwise_miss_rates_by_locus <- pairwise_miss_rates_by_locus

  ret$background$plots$non_miss_counts_by_indiv_plot <- non_miss_counts_by_indiv_plot
  ret$background$plots$pairwise_non_miss_counts_plot <- pairwise_non_miss_counts_plot

  # bail out if only tabulating and exiting
  if(tabulate_and_exit) {
    return(ret)
  }



  #### If we are simulating using pairwise patterns, we can just do it here
  if(simulation_approach == "pairwise_patterns") {
    # deal with ... stuff
    dotL <- list(...)
    bad_params_logi <- !(names(dotL) %in% c("sim_relats", "calc_relats", "reps"))
    if(sum(bad_params_logi) > 1) {
      bad_params <- names(dotL)[bad_params_logi]
      message("Ignoring ... params: ", paste(bad_params, collapse = ", "))
    }
    good_params <- dotL[!bad_params_logi]

    # first, we simulate the Qijs unlinked
    Qijs_unlinked <- parallel::mclapply(
      LTpairs,  #### Just do first 20 for testing!
      function(x) {
        rp <- good_params$reps
        MMM <- make_miss_mask_mat(MG, x, rp)
        gp_unlinked <- good_params
        gp_unlinked$reps <- nrow(MMM)
        plist1 <- list(
          C = C,
          miss_mask_mat = MMM
        )
        plist <- c(plist1, gp_unlinked)
        do.call(simulate_Qij, plist)
      },
      mc.cores = num_cores
    )

    # then we simulate the Qijs linked
    Qijs_linked <- parallel::mclapply(
      LTpairs,
      function(x) {
        rp <- good_params$reps
        MMM <- make_miss_mask_mat(MG, x, rp)
        gp_linked <- good_params
        gp_linked$reps <- nrow(MMM)
        plist1 <- list(
          C = C,
          miss_mask_mat = MMM,
          unlinked = FALSE,
          pedigree_list = pedigrees
        )
        plist <- c(plist1, gp_linked)
        do.call(simulate_Qij, plist)
      },
      mc.cores = num_cores
    )


    # then organize those into some list columns
    Qtib <- tibble(
      num_non_missing_loci = as.integer(names(Qijs_unlinked)),
      Qijs_unlinked = Qijs_unlinked,
      Qijs_linked = Qijs_linked
    )


    ret$Qij <- Qtib %>%
      mutate(num_missing_loci = L - num_non_missing_loci, .after = num_non_missing_loci)

    return(ret)
  }





  #### Look at the range of pairwise non-missing values  ####
  rg <- range(pairwise_non_miss_counts$num_non_miss)
  span <- rg[2] - rg[1]
  step <- span / (num_points - 1)
  sim_L_vals <- round(seq(from = rg[1], to = rg[2], by = step))

  # get that as the number of missing loci too
  sim_Miss_vals <- L - sim_L_vals
  names(sim_Miss_vals) <- sim_Miss_vals


  #### Run simulations at the multiple sim_Miss_vals, in parallel ####
  if(simulation_approach == "simple_miss_freq") {
    # deal with ... stuff
    dotL <- list(...)
    bad_params_logi <- !(names(dotL) %in% c("sim_relats", "calc_relats", "reps"))
    if(sum(bad_params_logi) > 1) {
      bad_params <- names(dotL)[bad_params_logi]
      message("Ignoring ... params: ", paste(bad_params, collapse = ", "))
    }
    good_params <- dotL[!bad_params_logi]

    # first, we simulate the Qijs unlinked
    Qijs_unlinked <- parallel::mclapply(
      sim_Miss_vals,
      function(x) {
        plist1 <- list(
          C = C,
          rando_miss_wts = pairwise_miss_rates_by_locus,
          rando_miss_n = x
        )
        plist <- c(plist1, good_params)
        do.call(simulate_Qij, plist)
      },
      mc.cores = num_cores
    )

    # then we simulate the Qijs linked
    Qijs_linked <- parallel::mclapply(
      sim_Miss_vals,
      function(x) {
        plist1 <- list(
          C = C,
          rando_miss_wts = pairwise_miss_rates_by_locus,
          rando_miss_n = x,
          unlinked = FALSE,
          pedigree_list = pedigrees
        )
        plist <- c(plist1, good_params)
        do.call(simulate_Qij, plist)
      },
      mc.cores = num_cores
    )


    # then organize those into some list columns
    Qtib <- tibble(
      num_missing_loci = as.integer(names(Qijs_unlinked)),
      Qijs_unlinked = Qijs_unlinked,
      Qijs_linked = Qijs_linked
    )


    ret$Qij <- Qtib %>%
      mutate(num_non_missing_loci = L - num_missing_loci, .before = num_missing_loci)

    return(ret)

  } else {
    stop("Not a recognized simulation approach: ", simulation_approach)
  }
}
