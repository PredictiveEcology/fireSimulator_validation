## Does the simulation put the years in the right order?

#' Spearman of the model median against the observation, with a permutation floor
#'
#' Spearman correlation, across scored years, between the observed value and the median of the simulated reps, with a
#' one-sided permutation p-value (years shuffled: the N0 floor). A `ref` series (e.g. an annual fire-weather index)
#' gives the higher bar: the same correlation for `ref`.
#'
#' @param d `data.table(year, obs, simMedian)`.
#' @param ref Optional `data.table(year, value)`: a reference predictor of the same years.
#' @param nPerm Permutations.
#' @return One-row `data.table(nYears, rho, pPerm, rhoRef)`.
scoreRanking <- function(d, ref = NULL, nPerm = 9999, seed = 1) {
  set.seed(seed)
  rho <- suppressWarnings(stats::cor(d$obs, d$simMedian, method = "spearman"))
  perm <- replicate(nPerm, suppressWarnings(stats::cor(sample(d$obs), d$simMedian, method = "spearman")))
  rhoRef <- if (is.null(ref)) NA_real_ else
    suppressWarnings(stats::cor(d$obs, ref$value[match(d$year, ref$year)], method = "spearman", use = "complete.obs"))
  data.table::data.table(nYears = nrow(d), rho = rho, pPerm = (1 + sum(perm >= rho, na.rm = TRUE)) / (nPerm + 1), rhoRef = rhoRef)
}
