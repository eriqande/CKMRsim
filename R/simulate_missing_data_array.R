

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
#' This function helps estimate error rates (FPRs/FNRs) under realistic missing data.
#' Some data sets have a lot of missing markers.  If this is the case, it is not
#' OK to just do the simulations as if there is no missing data.  This function
#' wraps up a lot of different steps that can be taken to try to get more accurate
#' "first-pass" FPRs and FNRs for situations with a lot of missing data.
#'
#' Note that this function requires a path to the output of `summarize_missing_data()`.
#' Rather than passing things in directly, we write an intermediate, uncompressed
#' file to a temp directory of the necessary inputs so that we can easily parallelize
#' this on a cluster if needed.

#'
#' Steps:
#' 1. Summarize missingness across individuals and pairs
#' 2. Estimate per-locus missing rates, then derive pairwise missing probabilities.
#' 3. Record background stats and plots
#' 4. Clean and prepare dot args
#' 5. Simulate Qij values using either random missingness (when `simulation_approach == "simple_miss_freq"`)
#'    or observed pairwise patterns of shared loci (when `simulation_approach == "pairwise_patterns"`).
#'
#' @param input_path the absolute path to the file of inputs written out by
#' `summarize_missing_data()`.
#' @inheritParams simulate_Qij
#' @param num_points Number of points between the observed min and max non-missing loci to simulate.
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
#' @param .parallel_plan Future strategy (e.g., "sequential", "multisession").
#' @param .workers Number of workers for parallel backend.
#' @param CKMRsim.discard_stderr logical.  Set to true by default to discard stderr,
#' @param CKMRsim.discard_stdout logical. Set to true by default to discard stdout
#' @inheritDotParams simulate_Qij sim_relats calc_relats reps
#' @importFrom future plan sequential
#' @importFrom future.apply future_lapply
#' @export
simulate_missing_data_array <- function(
    input_path,
    num_points = 11,
    tabulate_and_exit = FALSE,
    simulation_approach = "simple_miss_freq",
    ...,
    .parallel_plan = "sequential",
    .workers = NULL,
    CKMRsim.discard_stderr = TRUE,
    CKMRsim.discard_stdout = TRUE
) {

  # Set up future parallel plan
  if (is.character(.parallel_plan)) {
    if (.parallel_plan == "sequential") {
      future::plan(future::sequential)
    } else if (.parallel_plan == "multisession") {
      future::plan(future::multisession, workers = .workers %||% num_cores)
    } else if (.parallel_plan == "multicore") {
      future::plan(future::multicore, workers = .workers %||% num_cores)
    } else {
      stop("Unknown .parallel_plan string: ", .parallel_plan)
    }
  } else {
    future::plan(.parallel_plan)
  }


  ####  clean and prepare dot args ####
  dotL <- list(...)
  bad_params_logi <- !(names(dotL) %in% c("sim_relats", "calc_relats", "reps"))
  if (sum(bad_params_logi) > 1) {
    message("Ignoring ... params: ", paste(names(dotL)[bad_params_logi], collapse = ", "))
  }
  good_params <- dotL[!bad_params_logi]


  #### Get the lower-triangle pairs list needed to know what to future_lapply() over
  inList <- read_rds(input_path)
  LTpairs <- inList$LTpairs
  L <- ncol(inList$IG)  # total number of loci
  pairwise_non_miss_counts <- inList$pairwise_non_miss_counts
  ret <- list()  # will return things in a list (more flexible, if we want to add return items later)
  rm(inList)


  #### Step 3a: simulate using observed pairwise missing patterns ####
  if (simulation_approach == "pairwise_patterns") {
    Qijs_unlinked <- future.apply::future_lapply(
      LTpairs,
      function(x) {
        options(
          CKMRsim.discard_stderr = TRUE,
          CKMRsim.discard_stdout = TRUE
        )
        Inputs <- read_rds(input_path)
        MG <- Inputs$IG
        C <- Inputs$C
        rp <- good_params$reps
        MMM <- make_miss_mask_mat(MG, x, rp)
        gp_unlinked <- good_params
        gp_unlinked$reps <- nrow(MMM)
        plist1 <- list(C = C, miss_mask_mat = MMM)
        plist <- c(plist1, gp_unlinked)
        do.call(simulate_Qij, plist)
      },
      future.seed = TRUE
    )

    Qijs_linked <- future.apply::future_lapply(
      LTpairs,
      function(x) {
        options(
          CKMRsim.discard_stderr = TRUE,
          CKMRsim.discard_stdout = TRUE
        )
        Inputs <- read_rds(input_path)
        MG <- Inputs$IG
        C <- Inputs$C
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
      future.seed = TRUE
    )

    Qtib <- tibble(
      num_non_missing_loci = as.integer(names(Qijs_unlinked)),
      Qijs_unlinked = Qijs_unlinked,
      Qijs_linked = Qijs_linked
    )

    ret$Qij <- Qtib %>%
      mutate(num_missing_loci = L - num_non_missing_loci, .after = num_non_missing_loci)

    # Reset plan (recommended in packages)
    future::plan(future::sequential)

    return(ret)
  }

  #### Step 2b: simulate using random missingness frequencies ####
  rg <- range(pairwise_non_miss_counts$num_non_miss)
  span <- rg[2] - rg[1]
  step <- span / (num_points - 1)
  sim_L_vals <- round(seq(from = rg[1], to = rg[2], by = step))
  sim_Miss_vals <- L - sim_L_vals
  names(sim_Miss_vals) <- sim_Miss_vals

  if (simulation_approach == "simple_miss_freq") {
    Qijs_unlinked <- future.apply::future_lapply(
      sim_Miss_vals,
      function(x) {
        options(
          CKMRsim.discard_stderr = TRUE,
          CKMRsim.discard_stdout = TRUE
        )
        Inputs <- read_rds(input_path)
        C <- Inputs$C
        pairwise_miss_rates_by_locus = Inputs$pairwise_miss_rates_by_locus
        plist1 <- list(
          C = C,
          rando_miss_wts = pairwise_miss_rates_by_locus,
          rando_miss_n = x
        )
        plist <- c(plist1, good_params)
        do.call(simulate_Qij, plist)
      },
      future.seed = TRUE
    )

    Qijs_linked <- future.apply::future_lapply(
      sim_Miss_vals,
      function(x) {
        options(
          CKMRsim.discard_stderr = TRUE,
          CKMRsim.discard_stdout = TRUE
        )
        Inputs <- read_rds(input_path)
        C <- Inputs$C
        pairwise_miss_rates_by_locus = Inputs$pairwise_miss_rates_by_locus
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
      future.seed = TRUE
    )

    Qtib <- tibble(
      num_missing_loci = as.integer(names(Qijs_unlinked)),
      Qijs_unlinked = Qijs_unlinked,
      Qijs_linked = Qijs_linked
    )

    ret$Qij <- Qtib %>%
      mutate(num_non_missing_loci = L - num_missing_loci, .before = num_missing_loci)

    # Reset plan (recommended in packages)
    future::plan(future::sequential)

    return(ret)
  } else {
    stop("Not a recognized simulation approach: ", simulation_approach)
  }


}
