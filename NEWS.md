# CKMRsim

# Version 0.1.3

*   Update pkgdown stuff
*   Add Mendelian incompatibility-tagging function
*   Use tibble instead of tbl_df
*   Also remove summarise that returns more than one row in some error checking of mendel-interface.R
*   Also count up some mendelian incompatilities in an example.
*   Added option to not write mendel standard output or standard error to screen.
    This can most easily be used by doing: `options(CKMRsim.discard_stderr = TRUE)` and/or
    `options(CKMRsim.discard_stdout = TRUE)` before running simulations using physical
    linkage.
*   Also, if you want to suppress the other messages coming from RCpp when using
    the linked marker simulation process, you can do `options(CKMRsim.linkage_verbosity = 0)`.
    Note that the other messages coming out of `simulate_Qij()` can be suppressed with
    `suppressMessages()`
* Passes CRAN checks.

# Version 0.1.2

*   Bug fixes in some of the RCpp types.
*   Rearranging column order in return tibble from close_matching_samples()
*   Fixing a C-stack issue that seemed to cause problems compiling the Example-1 vignette
*   Add the half-aunt-niece relationship
*   Add downloading of Mendel via R function.  (Big improvement in reproducibility)


# Version 0.1.1

## Changes

* Pull request #4 from Brage Førland (4/30/21):
    - Bugfix: `find_close_matching_samples()` didn't return the tibble!
    - Enhancement: With `denom` left as NULL, allow `pairwise_kin_logl_ratios()`
      to return the pairwise log-likleihoods (instead of always log-likelihood ratios).
    - Improvement: Modify `create_integer_genotype_matrix()` to account for loci
      included in the allele freqs tibble, but missing in a collection of compared individuals.
 
