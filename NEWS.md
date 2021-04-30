# CKMRsim 0.1.1

## Changes

* Pull request #4 from Brage Førland (4/30/21):
    - Bugfix: `find_close_matching_samples()` didn't return the tibble!
    - Enhancement: With `denom` left as NULL, allow `pairwise_kin_logl_ratios()`
      to return the pairwise log-likleihoods (instead of always log-likelihood ratios).
    - Improvement: Modify `create_integer_genotype_matrix()` to account for loci
      included in the allele freqs tibble, but missing in a collection of compared individuals.
 
