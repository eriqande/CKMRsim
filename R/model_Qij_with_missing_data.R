
#' Model the results of `simulate_missing_data_array()` to interpolate FNR and FPR values
#'
#' More later
#' @param Qtib The tibble output as the `Qij` component of the return list from
#' `simulate_missing_data_array()`.  This is a tibble with three columns:
#' `num_non_missing_loci` (int), `num_missing_loci` (int), and `Qijs` (a tibble
#' of Qij data frames).
#' @inheritParams mc_sample_simple
#' @export
#' @examples
#' # here for testing at the moment
#' Qtib <- read_rds("/tmp/Qtib.rds")
#'
model_Qij_with_missing_data <- function(
  Qtib,
  nu,
  de = "U",
  tr = "U",
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01)
) {


  MC <- Qtib %>%
    mutate(
      mc = map2(
        .x = Qijs_unlinked,
        .y = Qijs_linked,
        .f = function(x, y) {
          mc_sample_simple(
            Q = x,
            nu = nu,
            de = de,
            tr = tr,
            FNRs = FNRs,
            Q_for_fnrs = y
          )}
      )
    )

  # now, unnest that stuff
  MCu <- MC %>%
    select(-Qijs_linked, -Qijs_unlinked) %>%
    unnest(cols = mc) %>%
    mutate(
      log10_FPR = log10(FPR),
      log10_FPR_plus_2_se = log10(FPR + 2 * se),
      log10_FPR_minus_2_se = log10(FPR - 2 * se),
      .after = FPR
    )

  g <- ggplot(MCu, aes(x = num_non_missing_loci, y = log10_FPR, colour = factor(FNR))) +
    geom_point()

  g

  # let's do a simple linear model on each FNR
  simp <- MCu %>%
    group_by(FNR) %>%
    nest() %>%
    mutate(
      simplm = map(.x = data, .f = function(x) lm(log10_FPR ~ num_non_missing_loci, data = x)),
      tidy = map(simplm, broom::tidy)
    )

  s2 <- simp %>%
    select(FNR, tidy) %>%
    unnest(cols = tidy) %>%
    mutate(term = str_replace_all(term, "[^a-zA-Z]", "")) %>%
    select(FNR, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate)


  g +
    geom_abline(
      data = s2,
      mapping = aes(intercept = Intercept, slope = numnonmissingloci, colour = factor(FNR))
    )


}
