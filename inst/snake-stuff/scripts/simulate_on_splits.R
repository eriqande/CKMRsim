if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
}


# Detect whether running inside Snakemake
snakemake_exists <- exists("snakemake")

if (snakemake_exists) {
  pairs_path   <- snakemake@input$pairs
  ckmr_path   <- snakemake@input$ckmr
  imat_path   <- snakemake@input$imat
  output_path <- snakemake@output$outf
  split_id    <- snakemake@params$spl
  relats      <- snakemake@params$relats
} else {
  # Local test case
  pairs_path   <- "resources/LTpairs/97.rds"
  ckmr_path   <- "resources/ckmr_object.rds"
  imat_path   <- "resources/integer_genotype_matrix.rds"
  output_path <- "results/logls/97/splits/11.rds"
  split_id    <- "11"
  relats <- c("PO", "FS", "HS", "HAN", "HFC", "U")
}


library(tidyverse)
library(CKMRsim)


# Placeholder for loading input data
all_pairs <- read_rds(pairs_path)
spl_pairs <- all_pairs[[split_id]]
rm(all_pairs)

ckmr <- read_rds(ckmr_path)
imat <- read_rds(imat_path)

message("Running simulation for split ", split_id)
message("Using pairs file: ", pairs_path)
message("Saving result to: ", output_path)
message("Using as relats: ", paste(relats, collapse = ", "))


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
make_miss_mask_mat_no_recyc <- function(MG, L) {
  NP <- length(L$r)
  Mat <- lapply(1:NP, function(i) {
    !as.logical(MG[L$r[i],] *  MG[L$c[i],])
  }) %>%
    do.call(what = rbind, args = .)
  Mat
}


# get the miss_mask_matrix
MMM <-  make_miss_mask_mat_no_recyc(imat, spl_pairs)

nreps <- nrow(MMM)

# simulate the unlinked version
Q <- simulate_Qij(
  C = ckmr,
  sim_relats = relats,
  calc_relats = relats,
  reps = nreps,
  miss_mask_mat = MMM
)


# simulate the unlinked version
options(
  CKMRsim.discard_stderr = TRUE,
  CKMRsim.discard_stdout = TRUE
)
Q_linked <- simulate_Qij(
  C = ckmr,
  sim_relats = relats,
  calc_relats = relats,
  reps = nreps,
  miss_mask_mat = MMM,
  unlinked = FALSE,
  pedigree_list = pedigrees
)

write_rds(
  list(
    Q_linked = Q_linked,
    Q_unlinked = Q
  ),
  file = output_path,
  compress = "xz"
)
