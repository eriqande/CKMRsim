#' simulate a distribution of fake chromosome lengths
#'
#' Sometimes you don't know the genome of your organism and you don't
#' know where in the genome your markers sit.  But you might have
#' an idea of how big the genome is and how many chromosomes it is
#' organized into.  In that case, you can get a sense for how much
#' an effect physical linkage will have on false positive rates, etc.,
#' by placing the markers randomly within a genome that approximates
#' your organisms genome.  The first step to that is simulating chromosome
#' sizes.  That is what this function does.  It is very simple.
#' @param n number of chromosomes
#' @param L length of genome in Gigabases
#' @param sl length of the smallest chromosome as a fraction of the largest chromosome
#' @return A list with three components:
#'   - `chrom_lengths`: a tibble with columns idx (chromsome index), chrom (chromosome name),
#'   scaled_length (length of the chromosome as a fraction of the whole genome), and num_bases.
#'   Note that the names of the chromosomes are like "fc01", "fc02", and so on.  fc stands for
#'   "fake chromosome."
#'   - `chrom_length_plot`: a ggplot object that shows the simulated lengths of the chromosomes.
#'   - `chrom_length_parameters`: the parameters used in the function.
#' @export
#' @examples
#' gcl <- geometric_chromo_lengths(
#'     n = 40,
#'     L = 4.6,
#'     sl = 0.5
#' )
geometric_chromo_lengths <- function(n, L, sl) {

  # we use a geometric-like function
  # basically, for chromosome 1, length is set at 1, and for chromosome i we
  # set length to (1-p)^d(i-1).

  # So, first, we find p so that the last chromosome has length sl.
  p <- 1 - exp(log(sl) / (n - 1))

  # now, make a tibble of lengths
  lengths <- tibble::tibble(
    idx = 1:n,
    chrom = paste0("fc", sprintf("%02d", idx)),
    scaled_length = (1 - p) ^ (idx - 1),
  ) %>%
    dplyr::mutate(num_bases = floor(L * 1e9 * (scaled_length / sum(scaled_length))))

  # make a plot of the lengths, too:
  chrom_length_plot <- ggplot2::ggplot(lengths, ggplot2::aes(x = idx, y = num_bases / 1e6)) +
    ggplot2::geom_col() +
    ggplot2::xlab("Chromosome number") +
    ggplot2::ylab("Megabases")

  list(
    chrom_lengths = lengths,
    chrom_length_plot = chrom_length_plot,
    chrom_length_parameters = c(
      n = n,
      L = L,
      sl = sl
    )
  )
}
