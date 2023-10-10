#' Identify loci in PO pairs that are Mendelian incompatible
#'
#' This is a handy function for taking parent-offspring pairs, typically identified
#' using CKMRsim, and then tagging the loci at which the observed/recorded genotypes
#' are not compatible with Mendelian inheritance.  This version is currently very
#' simple and merely looks at pairs of individuals. (i.e., it does not attempt to identify
#' Mendelian incompatible loci on the basis of pa-ma-kid trios.)  Noetheless, it should
#' be useful for identifying poor-performing loci.
#' @param po a tibble of the parent-offspring pairs.  This will typically be obtained
#' by filtering the output for `pairwise_kin_logl_ratios()` and it must have the format
#' of that output: columns of `D2_indiv`, `D1_indiv`, `logl_ratio`, and  `num_loc`.
#' @param geno a long format data frame / tibble of genotype data for both the candidate
#' offspring and the candidate parents.  Typically this is what you get if you use
#' `bind_rows()` on the `D1` and `D2` inputs to `pairwise_kin_logl_ratios()` for a parent-offspring
#' type of analysis.  Must have the
#' columns: `Indiv`, `Locus`, `gene_copy`, and `Allele`. Obviously, the IDs in Indiv must
#' correspond to those used in `D1_indiv` and `D2_indiv`
#' @return This returns a tibble with the columns:
#' - `D2_indiv`: the ID of the D2_indiv
#' - `D1_indiv`: the ID of the D1_indiv
#' - `logl_ratio`: the log-likelihood ratio for parent-offspring for the pair
#' - `Locus`: the name of the locus
#' - `D1a_1, D1a_2, D2a_1, D2a_2`: D1_indiv's allelic type of gene copy 1 and 2, and then
#' D2_indiv's allelic type of gene copy 1 and 2, respectively
#' - `is_MI`: logical vector.  TRUE if the pair is Mendelian incompatible at the Locus. Otherwise
#' FALSE (or NA if either member of pair was NA)
#' @details Note that if you want to count up the relative frequency (across all pairs) of
#' Mendelian incompatibilities at each locus you can do like this:
#' ```r
#' tag_mendelian_incompatibilities(po, geno) %>%
#'  group_by(Locus) %>%
#'  summarise(mean_incompat = mean(mend_incompat, na.rm = TRUE)) %>%
#'  arrange(desc(mean_incompat))
#' ```
#' @export
tag_mendelian_incompatibilities <-  function(po, geno) {
  ret <- left_join(po, geno, by = c("D2_indiv" = "Indiv"), relationship = "many-to-many") %>%
    rename(D2a = Allele) %>%
    select(-num_loc) %>%
    left_join(geno, by = c("D1_indiv" = "Indiv", "Locus", "gene_copy")) %>%
    rename(D1a = Allele) %>%
    pivot_wider(names_from = gene_copy, values_from = c(D1a, D2a), names_sep = "_") %>%
    mutate(
      is_MI = D1a_1 != D2a_1 & D1a_1 != D2a_2 & D1a_2 != D2a_1 & D1a_2 != D2a_2
    )

  ret
}
