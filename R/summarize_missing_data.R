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
#' @export
#' @examples
#' set.seed(1)
#' M <- matrix(sample(80:180, size = 100^2, replace = TRUE), nrow = 100)
#' result <- collect_lower_triangle_pairs(M)
#' str(result[["100"]])
#'
collect_lower_triangle_pairs <- function(M) {
  stopifnot(is.matrix(M))
  stopifnot(nrow(M) == ncol(M))

  # Use which() on the lower triangle to get indices
  lower_indices <- which(lower.tri(M), arr.ind = TRUE)
  values <- M[lower_indices]

  # Split the row and column indices by the value in the matrix
  split_indices <- split(as.data.frame(lower_indices), values)

  # Format the result as requested: list named by V*, each element is a list with $r and $c
  result <- lapply(split_indices, function(df) list(r = df$row, c = df$col))

  return(result)
}




######## END OF HELPER FUNCTIONS #################



#' Summarize the missing data and make some plots.
#'
#' In default mode, this will just summarize the missing data and make and
#' store some plots in the list output. If you provide `snakemake_dir`, this
#' function writes a self-contained Snakemake arena for running the missing-data
#' simulations in separate R processes.
#'
#' Note that C must be a CKMR object that is suitable for *both* linked and unlinked
#' simulation.  So, if you don't have that, then you best
#' prepare it, even if it means sprinkling your markers into a pseudogenome
#' as described [here](https://eriqande.github.io/tws-ckmr-2022/kin-finding-lab.html#power-for-kin-finding-while-accounting-for-physical-linkage).
#'
#' @return A list containing summary plots and statistics, or the absolute path
#' to the Snakemake directory when `snakemake_dir` is provided.
#' @inheritParams create_integer_genotype_matrix
#' @param C A CKMR object created by `create_ckmr()`.
#' @param snakemake_dir Name of directory to create in order to write out the
#' materials for running the simulations via a Snakefile.
#' @param snake_rep_split If snakemake_dir is non-NA, this is the number of simulation
#' reps to be done for each snakemake job.  If a partition has fewer pairs that snake_rep_split
#' of if snake_rep_split is not perfectly divisible by the number of pairs
#' then the pairs are recycled so that every job has exactly snake_rep_split reps.
#' @export
summarize_missing_data <- function(
    LG,
    C,
    snakemake_dir = NA,
    snake_rep_split = 5e4
) {
  ret <- list()

  #### Step 1: summarize missingness across individuals and pairs ####
  AF <- C$orig_data
  MG <- create_integer_genotype_matrix(LG, AF)
  L <- ncol(MG)  # number of loci

  MG[MG >= 0] <- 1  # recode genotypes as 1 (non-missing)
  MG[MG == -1] <- 0 # missing values become 0
  big <- MG %*% t(MG)

  pairwise_non_miss_counts <- tibble(
    num_non_miss = big[lower.tri(big)]
  ) %>%
    count(num_non_miss)

  pairwise_non_miss_counts_plot <- ggplot(pairwise_non_miss_counts, aes(x = num_non_miss, y = n)) +
    geom_col() +
    xlab("Number of loci missing in neither pair member") +
    ylab("Number of pairs")

  non_miss_counts_by_indiv <- tibble(
    num_non_miss = rowSums(MG)
  ) %>%
    count(num_non_miss)

  non_miss_counts_by_indiv_plot <- ggplot(non_miss_counts_by_indiv, aes(x = num_non_miss, y = n)) +
    geom_col() +
    xlab("Number of non-missing loci") +
    ylab("Number of individuals")

  # Collect the lower triangle so Snakemake can run each shared-locus-count
  # category separately.
  LTpairs <- collect_lower_triangle_pairs(big)
  rm(big)

  #### Step 2: estimate missingness rates ####
  miss_rates_by_locus <- colMeans(MG == 0)
  pairwise_miss_rates_by_locus <- 1 - (1 - miss_rates_by_locus)^2

  #### Step 3: record background stats and plots ####
  ret$background$values$pairwise_non_miss_counts <- pairwise_non_miss_counts
  ret$background$values$non_miss_counts_by_indiv <- non_miss_counts_by_indiv
  ret$background$values$miss_rates_by_locus <- miss_rates_by_locus
  ret$background$values$pairwise_miss_rates_by_locus <- pairwise_miss_rates_by_locus
  ret$background$values$lower_triangle_pairs <- LTpairs
  ret$background$values$indiv_names_tibble <- tibble::tibble(Indiv = rownames(MG)) %>%
    mutate(idx = 1:n(), .before = Indiv)

  ret$background$plots$non_miss_counts_by_indiv_plot <- non_miss_counts_by_indiv_plot
  ret$background$plots$pairwise_non_miss_counts_plot <- pairwise_non_miss_counts_plot

  # return from here if that is all we are doing.
  if (is.na(snakemake_dir)) return(ret)

  # in this case we create a directory and write some files out to it that will
  # make it easy to run all the simulations with Snakemake
  dir.create(snakemake_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(snakemake_dir, "resources"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(snakemake_dir, "resources", "LTpairs"), recursive = TRUE, showWarnings = FALSE)

  # now, break the LTpairs into reasonably sized chunks
  chunks <- break_up_LTpairs(LTpairs, R = snake_rep_split)

  # Save each element of chunks to resources/LTpairs and capture the names of
  # all the chunks for Snakemake as well.
  chunks_tibble <- lapply(
    names(chunks),
    function(n) {
      # write out the rds files with the pairs in them
      write_rds(
        chunks[[n]],
        file = file.path(
          snakemake_dir, "resources", "LTpairs",
          paste(n, ".rds", collapse = "", sep = "")
        )
      )

      # return tibbles with the number of loci and the splits in it
      # to give to Snakemake later
      nsplit <- length(chunks[[n]])
      tibble(
        num_loci = rep(n, nsplit)
      ) %>%
        mutate(splits = 1:nsplit)
    }) %>%
    bind_rows()


  # then save the other things that we will need
  write_rds(C, file.path(snakemake_dir, "resources", "ckmr_object.rds"))
  write_rds(MG, file.path(snakemake_dir, "resources", "integer_genotype_matrix.rds"))
  write_rds(ret$background$values$indiv_names_tibble, file.path(snakemake_dir, "resources", "indiv_names_tibble.rds"))
  write_tsv(chunks_tibble, file = file.path(snakemake_dir, "resources", "chunks_tibble.tsv"))

  #### then write the Snakefile and the scripts that it will use ####
  relnames <- paste(paste0('"', names(C$loci[[1]]$X_l), '"'), collapse = ", ")  # get the relationships that are in the CKMR object, C
                                                              # this variable gets glued into the Snakefile

  # 1. Locate and read the Snakefile-skeleton
  skeleton_path <- system.file("snake-stuff/Snakefile-skeleton", package = "CKMRsim")
  skeleton_contents <- readLines(skeleton_path) %>%
    paste(collapse = "\n")

  # 2. Replace the line with glue interpolation
  rendered_contents <- glue::glue_collapse(glue::glue(skeleton_contents, .open = "<<", .close = ">>"), sep = "\n")

  # 3. Write the result to a Snakefile in snakemake_dir
  snakefile_path <- file.path(snakemake_dir, "Snakefile")
  writeLines(rendered_contents, con = snakefile_path)

  # 4. Copy the scripts directory into snakemake_dir
  scripts_dir <- system.file("snake-stuff/scripts", package = "CKMRsim")
  fs::dir_copy(scripts_dir, file.path(snakemake_dir, "scripts"), overwrite = TRUE)

  # go ahead and return the directory
  return(normalizePath(snakemake_dir))


}
