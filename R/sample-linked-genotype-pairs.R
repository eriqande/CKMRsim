#' simulate genotype-pairs from linked markers
#'
#' This is the main function to use.  Pass it a data frame of markers (indexed in order)
#' and it will drop the genes, do genotyping error and then return
#' a list of vectors which hold the genotype index that you can use to subscript the
#' joint probability vectors with.  This function assumes that the loci in df are ordered
#' appropriately (i.e have been run through reindex_markers()) and that the components in
#' C are named Chrom.Locus.Pos, as is typical.  Obviously the C list should correspond
#' exactly to the markers/alleles in df.
#' @param df  A data frame in the format of \code{\link{long_markers}}.
#' @param ped  The pedigree to be simulating from
#' @param C a list whose elements contain, at a minimum, the "C-matrices" which give the
#' probability of observed genotypes given true genotypes.  If this is NULL (the default), then the
#' function will assume no genotyping error.
#' @param num Number of reps of gene-dropping to do. Default is 1000
#' @param useMendel If TRUE, use the external Mendel program for linked gene
#' dropping. If FALSE, use CKMRsim's internal Rcpp simulator.
#' @param cM_per_Mb Recombination-rate conversion used to create a Map column
#' from Pos when Map is not present. The value is centiMorgans per megabase.
#' @param min_crossovers Minimum number of crossovers in the typed marker span
#' of each chromosome for each meiosis.
#' @return  This returns a list named by Chrom.Locus.Pos (the names of C), in which each
#' component is a vector of length num that is the integer index of the simulated pair of genotypes.
sample_linked_genotype_pairs <- function(
    df,
    ped,
    C = NULL,
    num = 1000,
    useMendel = FALSE,
    cM_per_Mb = 1,
    min_crossovers = 0
) {

  # first get the number of alleles at each locus
  alle_nums <- df %>%
   dplyr::mutate(list_name = paste(Chrom, Locus, Pos, sep = ".")) %>%
   dplyr::mutate(list_name = factor(list_name, levels = unique(list_name))) %>%
    dplyr::group_by(list_name) %>%
    dplyr::tally()

  # check to make sure that the names are correct here
  if(!all(names(C) == alle_nums$list_name)) stop("Mismatch between Chrom.Locus.Pos names in df and in C")

  if(useMendel == TRUE) {
    message("Simulating linked markers with MENDEL")
    outgenos <- sample_linked_true_genotype_pairs_mendel(df = df, ped = ped, num = num)
  } else {
    message("Simulating linked markers with internal Rcpp gene dropper")
    outgenos <- sample_linked_true_genotype_pairs(
      df = df,
      ped = ped,
      num = num,
      cM_per_Mb = cM_per_Mb,
      min_crossovers = min_crossovers
    )
  }

  # at this juncture, we have the genotypes of each individual at all the loci, but we
  # need to apply the true genotyping errors to them.  This we do by
  # lapplying over the loci (as the names of C).
  G1 <- outgenos$indiv1
  G2 <- outgenos$indiv2
  locs <- names(C)
  names(locs) <- locs


  obs_geno_pairs <- lapply(locs, function(n) {
    g1 <- G1[, n] # simulated genotypes of indiv1
    g2 <- G2[, n]
    Clt <- C[[n]]$C_l_true
    g1e <- samp_from_mat(Clt[g1,])   # Ctl[g1,] is a matrix where each row corresponds to the probs of observed genos given the true geno of indiv1
    g2e <- samp_from_mat(Clt[g2,])

    ## NOTE: IF I WANTED TO HAVE A FEATURE THAT ALLOWED WRITING OUT A FILE OF GENOTYPES
    ## SIMULATED WITH PHYSICAL LINKAGE, THEN THIS WOULD BE THE PLACE TO DO IT...

    # here is the number of genos
    nG <- nrow(Clt)

    # so, the vector-position of the genotype of the two individuals in a matrix of
    # joint probabilities will be  nG * (g2 - 1) + g1
    nG * (g2e - 1) + g1e
  })

  obs_geno_pairs
}
