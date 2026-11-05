#' Share of hindcast years whose observed value lies inside the runs' central range
#'
#' For each per-year metric of one ELF (area burned by fires of `MINHA` ha or more, ignitions, escapes): the share of
#' scored years in which the observed value is inside the central `level` range of that year's runs ([insideInterval()]).
#' A calibrated model gives about `level`. Years without an observed value are left out.
#'
#' @param elf One ELF input ([runWave1()]).
#' @param level The central range, e.g. 0.95 (2.5% to 97.5% of the runs).
#' @return Named numeric (0-1, NA with no scored year): `area`, `ignitions`, `escapes`.
yearCoverage <- function(elf, level = 0.95) {
  yrs <- scoredYears(elf)
  vapply(c(area = "area50", ignitions = "ignitions", escapes = "escapes"), function(col) {
    sim <- simBy(elf, col, yrs); obs <- obsBy(elf, col, yrs); ok <- is.finite(obs)
    if (!any(ok)) NA_real_ else mean(mapply(insideInterval, sim[ok], obs[ok], MoreArgs = list(level = level)))
  }, 0)
}

## the yearCoverage() value each metric figure of the page is mapped by; fireSizes has none (not a per-year metric)
coverageOfMetric <- c(areaByYear = "area", annualArea = "area", ignitions = "ignitions", escapes = "escapes")
