

#' re-index loci and or alleles  given a data frame of markers
#'
#' The data frame should be of the format of \code{\link{long_markers}}.
#' This function will reset the values of AlleIdx and LocIdx in case
#' some markers (or some alleles) were removed.  It also rescales the
#' allele frequencies to sum to one, and, of course, sorts the loci into
#' map order along the chromosomes.
#' @param M a data frame of markers.
#' @export
#' @examples
#' # subsample alleles from the long_markers data frame.
#' # Note, you would never subsample alleles like this but it
#' # helps to demonstrate the function.
#' set.seed(5);  # for reproducibility
#' A_few_markers <- long_markers %>%
#'   dplyr::sample_n(1000) %>%     # subsample 4000 alleles
#'   dplyr::group_by(Chrom, Locus) %>%
#'   dplyr::mutate(NumAlle = n()) %>%
#'   dplyr::filter(NumAlle > 1)  %>%  # chuck out loci only one allele now...
#'   dplyr::select(-NumAlle)  # get rid of the NumAlle column we created
#'
#'   reindex_markers(A_few_markers)   # then reindex them
reindex_markers <- function(M) {
  M %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Chrom, Pos, desc(Freq)) %>%
    dplyr::group_by(Chrom) %>%
    dplyr::mutate(locidx = as.integer(factor(Locus, levels = unique(Locus)))) %>%
    dplyr::group_by(Chrom, Locus) %>%
    dplyr::mutate(alleidx = as.integer(factor(Allele, levels = unique(Allele))),
                  newfreq = Freq / sum(Freq)
    ) %>%
    dplyr::select(-AlleIdx, -LocIdx, -Freq) %>%
    rename(Freq = newfreq, AlleIdx = alleidx, LocIdx = locidx) %>%
    dplyr::ungroup()
}
