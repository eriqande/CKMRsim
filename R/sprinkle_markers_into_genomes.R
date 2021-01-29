
#' Randomly distribute markers along the chromosomes of a genome
#'
#' If you have a bunch of markers, but you don't know where in the genome
#' they are, BUT you have an idea of how big the genome is and how many
#' chromosomes are in it, this will place your markers uniformly within
#' the length of the genome.  This will let you assess how much of an effect
#' the physical linkage might have on pairwise relationship inference.
#' @param markers A tibble of markers and allele frequencies like `markers40`. It
#' should have the columns Chrom, Locus, Pos, Allele, LocIdx, AlleIdx, and Freq,
#' AND it must have been run through `reindex_markers()`.
#' @param genome A tibble like that in the
#' `chromo_lengths` component of the list returned by `geometric_chromo_lengths()`. It
#' must have the columns idx, chrom, scaled_length, num_bases.
#' @return This returns a tibble that is of the same format as `markers`, but the positions
#' of the markers have been simulated uniformly into the length of the genome in `genome`.
#' Note that the `Locus` column in `markers` is retained as the name of each locus.
#' @export
#' @examples
#' # make a fake genome
#' genome = geometric_chromo_lengths(n = 41, L = 7.4, sl = 0.2)$chrom_lengths
#' # get the markers to sprinkle into that genome
#' markers = markers40
sprinkle_markers_into_genome <- function(markers, genome) {

  # record how many markers there are. This technically is the number
  # of distinct Locus names on each Chrom, summed across Chrom's. (Note
  # that when things get reindexed the LocIdx starts over on every chromosome).
  Ltib <- markers %>%
    dplyr::group_by(Chrom) %>%
    dplyr::summarise(n = dplyr::n_distinct(Locus))

  L <- sum(Ltib$n)

  # simulate the number of markers that will be on each chromosome,
  # and filter out any chromosomes that don't have any markers on them,
  # and then simulate positions in each chromomes uniformly,
  # and, finally, randomly sample an AbsoluateIdex to each of those for joining
  # the allele frequencies onto it.
  G <- genome %>%
    dplyr::mutate(
      num_markers = rmultinom(
        n = 1,
        size = L,
        prob = scaled_length
      )[,1]
    ) %>%
    dplyr::filter(num_markers > 0) %>%
    dplyr::mutate(
      pos = purrr::map2(
        .x = num_markers,
        .y = num_bases,
        .f = function(x, y) {
          floor(
            runif(  # note: we assume here that no two will be in the same position...
              n = x,
              min = 0,
              max = y
            )
          )
        }
      )
    ) %>%
    dplyr::select(chrom, pos) %>%
    tidyr::unnest(pos) %>%
    dplyr::mutate(AbsoluteIndex = sample(x = 1:L, size = L))

  # now, join the markers on there according to LocIdx to get the new genome
  # positions
  M <- markers %>%
    dplyr::select(Chrom, Locus, Allele, LocIdx, AlleIdx, Freq) %>%
    dplyr::mutate(AbsoluteIndex = as.integer(as.factor(
      paste0(Chrom, Locus)
    ))) %>%
    dplyr::select(-Chrom) %>%
    dplyr::left_join(G, by = "AbsoluteIndex")

  # Now, we just move the names around and reindex it, and return that.
  M2 <- M %>%
    dplyr::rename(Chrom = chrom, Pos = pos) %>%
    dplyr::select(Chrom, Locus, Pos, dplyr::everything()) %>%
    reindex_markers(.) %>%
    dplyr::select(-AbsoluteIndex)

  M2

}
