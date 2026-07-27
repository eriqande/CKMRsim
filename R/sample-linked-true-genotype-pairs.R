#' simulate true genotype-pairs from linked markers using Rcpp
#'
#' @param df A data frame in the format of \code{\link{long_markers}}.
#' @param ped The pedigree to be simulating from.
#' @param num Number of reps of gene-dropping to do. Default is 1000.
#' @param cM_per_Mb Recombination-rate conversion used to create a Map column
#' from Pos when Map is not present. The value is centiMorgans per megabase.
#' @param min_crossovers Minimum number of crossovers in the typed marker span
#' of each chromosome for each meiosis.
#' @keywords internal
sample_linked_true_genotype_pairs <- function(
    df,
    ped,
    num = 1000,
    cM_per_Mb = 1,
    min_crossovers = 0
) {
  inp <- linked_marker_cpp_inputs(df, cM_per_Mb = cM_per_Mb)

  out <- linked_gene_drop_true_genotypes(
    kid = as.integer(ped$Kid),
    pa = as.integer(ped$Pa),
    ma = as.integer(ped$Ma),
    n_alleles = as.integer(inp$n_alleles),
    map_pos = inp$map_pos,
    chrom_id = as.integer(inp$chrom_id),
    chrom_start = as.integer(inp$chrom_start),
    chrom_end = as.integer(inp$chrom_end),
    freqs = inp$freqs,
    reps = as.integer(num),
    min_crossovers = as.integer(min_crossovers)
  )

  lapply(out, function(x) {
    dimnames(x) <- list(NULL, locus = inp$list_name)
    x
  })
}
