#' prepare linked marker inputs for the Rcpp gene dropper
#'
#' @param df A data frame in the format of \code{\link{long_markers}}.
#' @param cM_per_Mb Recombination-rate conversion used to create a Map column
#' from Pos when Map is not present. The value is centiMorgans per megabase.
#' @keywords internal
linked_marker_cpp_inputs <- function(df, cM_per_Mb = 1) {
  needed <- c("Chrom", "Locus", "Pos", "AlleIdx", "LocIdx", "Freq")
  miss <- setdiff(needed, names(df))
  if(length(miss) > 0) {
    stop("Missing columns in df for linked simulation: ", paste(miss, collapse = ", "))
  }
  if(!("Map" %in% names(df))) {
    df$Map <- df$Pos / 1e6 * cM_per_Mb
  }

  locs <- df %>%
    dplyr::arrange(Chrom, LocIdx, AlleIdx) %>%
    dplyr::group_by(Chrom, Locus, Pos, LocIdx) %>%
    dplyr::summarise(
      Map = dplyr::first(Map),
      map_n = dplyr::n_distinct(Map),
      n_alleles = dplyr::n(),
      freqs = list(Freq),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Chrom, Map, LocIdx)

  if(any(locs$map_n != 1)) {
    stop("Each allele row for a locus must have the same Map value.")
  }

  list_name <- paste(locs$Chrom, locs$Locus, locs$Pos, sep = ".")
  if(anyDuplicated(list_name)) {
    stop("Chrom.Locus.Pos names are not unique.")
  }
  if(any(!is.finite(locs$Map))) {
    stop("Map positions must be finite.")
  }

  chrom_fac <- factor(locs$Chrom, levels = unique(locs$Chrom))
  split_idx <- split(seq_len(nrow(locs)), chrom_fac)
  chrom_start <- vapply(split_idx, min, integer(1))
  chrom_end <- vapply(split_idx, max, integer(1))
  bad_chroms <- vapply(split_idx, function(i) any(diff(locs$Map[i]) < 0), logical(1))
  if(any(bad_chroms)) {
    stop("Map positions must be nondecreasing within each chromosome.")
  }

  list(
    list_name = list_name,
    n_alleles = locs$n_alleles,
    map_pos = locs$Map,
    chrom_id = as.integer(chrom_fac),
    chrom_start = chrom_start,
    chrom_end = chrom_end,
    freqs = locs$freqs
  )
}
