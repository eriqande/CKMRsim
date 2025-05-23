if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
}


# Detect whether running inside Snakemake
snakemake_exists <- exists("snakemake")

if (snakemake_exists) {
  split_files <- snakemake@input$infs
  outf_linked <- snakemake@output$outf_linked
  outf_unlinked <- snakemake@output$outf_unlinked
} else {
  # Local test case
  split_files <- list.files("results/logls/160/splits", pattern = "\\.rds$", full.names = TRUE)
  outf_linked <- "results/logls/160/Q_linked.rds"
  outf_unlinked <- "results/logls/160/Q_unlinked.rds"
}

library(tidyverse)
library(CKMRsim)

message("Aggregating ", length(split_files), " files.")
message("Output will be written to: ", outf_linked, " and ", outf_unlinked)

# Do it
all_splits <- lapply(split_files, read_rds)


write_rds(
  merge_Qij(Qlist = lapply(all_splits, function(x) x$Q_linked)),
  file = outf_linked,
  compress = "xz"
)

write_rds(
  merge_Qij(Qlist = lapply(all_splits, function(x) x$Q_unlinked)),
  file = outf_unlinked,
  compress = "xz"
)

