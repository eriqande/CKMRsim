
#' turns long format genotypes into a matrix of integer-code genotypes
#'
#' @param LG the genotypes in long format.  It must have the columns
#' Indiv (unique IDs of the individuals), Locus, gene_copy (must be 1 or 2
#' denoting which of the two gene copies in a diploid each allele is), and
#' Allele, which must be a character. If there are any missing genotypes
#' in the data frame, they must appear as NAs in the Allele column.
#' @param AF The data frame of allele frequencies that was passed
#' to create_ckmr() to create the ckmr object that will be used for computing
#' genotype probabilities.  It must have columns
#' Chrom, Pos, Locus, Allele, LocIdx, AlleIdx, and Freq.
#' @export
create_integer_genotype_matrix <- function(LG, AF) {

  # first check that all the columns are right
  if(length(intersect(names(LG), c("Indiv", "Locus", "gene_copy", "Allele"))) != 4) {
    stop("Missing columns in LG.  Should have Indiv, Locus, gene_copy, and Allele")
  }
  if(length(setdiff(names(LG), c("Indiv", "Locus", "gene_copy", "Allele"))) != 0) {
    stop("Extra columns in LG.  Should only have Indiv, Locus, gene_copy, and Allele")
  }

  if(length(intersect(names(AF), c("Chrom", "Pos", "Locus", "Allele", "LocIdx", "AlleIdx", "Freq"))) != 7) {
    stop("Missing columns in AF.  Should have Chrom, Pos, Locus, Allele, LocIdx, AlleIdx, Freq")
  }
  if(length(setdiff(names(AF), c("Chrom", "Pos", "Locus", "Allele", "LocIdx", "AlleIdx", "Freq"))) != 0) {
    stop("Missing columns in AF.  Should have Chrom, Pos, Locus, Allele, LocIdx, AlleIdx, Freq")
  }


  # filter out any missing data in it.
  LG2 <- LG %>%
    dplyr::filter(!is.na(Allele))

  # now, make a data frame of genotype indexes
  LGI_wide_frame <- AF %>%
    dplyr::select(Locus, Allele, LocIdx, AlleIdx) %>%
    dplyr::group_by(Locus) %>%
    dplyr::mutate(NumA = dplyr::n()) %>%  # get the number of alleles at each locus
    dplyr::ungroup() %>%
    dplyr::left_join(LG2, ., by = c("Locus", "Allele")) %>%  # join the alle_idx's onto the actual genotype data
    dplyr::select(Indiv, Locus, gene_copy, LocIdx, NumA, AlleIdx) %>%
    tidyr::spread(key = gene_copy, value = AlleIdx) %>%
    dplyr::mutate(GenoIdx = index_ab(a = `1`, b = `2`, A = NumA)) %>%  # get the genotype indexes
    dplyr::ungroup() %>%
    dplyr::select(Indiv, LocIdx, GenoIdx) %>%
    tidyr::spread(data = ., key = LocIdx, value = GenoIdx) # now spread it into a data frame that looks like the matrix we want

  # now, turn that into a base-0 integer matrix
  mat <- as.matrix(LGI_wide_frame[, -1])
  rownames(mat) <- LGI_wide_frame$Indiv
  mat[is.na(mat)] <- 0
  mat <- mat - 1
  storage.mode(mat) <-  "integer"

  mat
}
