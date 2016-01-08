

#' re-index loci and or alleles  given a data frame of markers
#'
#' The data frame should be of the format of \code{\link{long_markers}}.
#' This function will reset the values of AlleIdx and LocIdx in case
#' some markers (or some alleles) were removed.  It also rescales the
#' allele frequencies to sum to one, and, of course, sorts the loci into
#' map order along the chromosomes.
#' @param M a data frame of markers.
#' @examples
#' # subsample alleles from the long_markers data frame.
#' # Note, you would never subsample alleles like this but it
#' # helps to demonstrate the function.
#' set.seed(5);  # for reproducibility
#' A_few_markers <- long_markers %>%
#'   sample_n(1000) %>%     # subsample 4000 alleles
#'   group_by(Chrom, Locus) %>%
#'   mutate(NumAlle = n()) %>%
#'   filter(NumAlle > 1)  %>%  # chuck out loci only one allele now...
#'   select(-NumAlle)  # get rid of the NumAlle column we created
#'
#'   reindex_markers(A_few_markers)   # then reindex them
reindex_markers <- function(M) {
  M %>%
    ungroup %>%
    arrange(Chrom, Pos, desc(Freq)) %>%
    group_by(Chrom) %>%
    mutate(locidx = as.integer(factor(Locus, levels = unique(Locus)))) %>%
    group_by(Chrom, Locus) %>%
    mutate(alleidx = as.integer(factor(Allele, levels = unique(Allele))),
           newfreq = Freq / sum(Freq)
           ) %>%
    select(-AlleIdx, -LocIdx, -Freq) %>%
    rename(Freq = newfreq, AlleIdx = alleidx, LocIdx = locidx) %>%
    ungroup
}
