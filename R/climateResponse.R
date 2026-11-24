## Climate response: does simulated annual area burned respond to the year's climate as the observed does?
## The fits' spread climate covariate (`elf$covariate`: the ELF median of the year's value, higher = drier, e.g. CMD) is standardised over the scored years.
## Definitions are those of the climate comparison of the extreme-fire-years analysis (climCompare.R); a run is one simulated history of the scored years.

climateMetrics <- c(slope = "climate slope", dry10 = "driest 10% of years", wetHalf = "wetter half of years")

#' Climate response of one ELF
#'
#' Three ratios of the simulated to the observed record over the scored years, each computed for every run alone. Area is `area50` (fires of the escape size or
#' more, as the area row of the page); a year a run (or the record) lacks counts as no area. The covariate (`elf$covariate`, higher = drier) is standardised over the scored years,
#' so the observed slope is per standard deviation of the covariate over the scored years.
#' \itemize{
#'   \item `slope`: the slope of log annual area burned on the standardised covariate, years with fire only; sim / obs.
#'   \item `dry10`: the area burned in the driest 10% of the scored years (years at or above the 90% quantile of the covariate over the scored years); sim / obs.
#'   \item `wetHalf`: the area burned in the wetter half of the scored years (covariate at or below its median); sim / obs.
#' }
#' The area ratios use [zeroSafeRatio()]. A slope is NA with fewer than 3 years with fire, and a ratio of slopes is NA when the observed slope is 0 (to within 1.5e-8) or NA. The single runs are summarised by [runSummary()].
#'
#' @param elf An ELF input list ([runWave1()]); without a `covariate` of the scored years there is no result (`NULL`).
#' @return `data.frame`, one row per metric: `level` ("ELF"), `ELF`, `metric`, `observed` (the observed slope, or area in ha), `simMedian` (median over runs of the run's value), `ratio` (median over runs of the
#'   run's ratio to the observed), `runLo`, `runHi` (2.5% and 97.5% of the single runs' ratio), `obsRank` (share of runs below the observed), `elfLo`, `elfHi` (NA here), `nELF` (1), `nNonFinite` (0), `nRuns`, `nYears`.
scoreClimateResponse <- function(elf) {
  years <- scoredYears(elf)
  cv <- covariateOf(elf, years)
  if (is.null(elf$covariate) || anyNA(cv) || length(years) < 3 || stats::sd(cv) == 0) return(NULL)
  z <- (cv - mean(cv)) / stats::sd(cv)
  dry <- z >= stats::quantile(z, 0.9); wet <- z <= stats::median(z)
  slopeOf <- function(area) { k <- !is.na(area) & area > 0; if (sum(k) < 3) NA_real_ else unname(stats::coef(stats::lm(log(area[k]) ~ z[k]))[2]) }
  values <- function(area) c(slope = slopeOf(area), dry10 = sum(area[dry], na.rm = TRUE), wetHalf = sum(area[wet], na.rm = TRUE))   # area: per scored year, in `years` order
  obs <- values(obsBy(elf, "area50", years))
  ratioTo <- function(v) c(slope = if (isTRUE(abs(obs[["slope"]]) > sqrt(.Machine$double.eps))) v[["slope"]] / obs[["slope"]] else NA_real_,
                           dry10 = zeroSafeRatio(v[["dry10"]], obs[["dry10"]]), wetHalf = zeroSafeRatio(v[["wetHalf"]], obs[["wetHalf"]]))
  runs <- lapply(split(elf$simYears, elf$simYears$run), function(r) values(r$area50[match(years, r$year)]))
  v <- do.call(rbind, runs); single <- vapply(runs, ratioTo, numeric(3)); rs <- runSummary(single)
  data.frame(level = "ELF", ELF = elf$elf, metric = names(climateMetrics), observed = unname(obs), simMedian = apply(v, 2, stats::median, na.rm = TRUE),
             ratio = apply(single, 1, medFinite), runLo = rs$runLo, runHi = rs$runHi, obsRank = rs$obsRank,
             elfLo = NA_real_, elfHi = NA_real_, nELF = 1L, nNonFinite = 0L, nRuns = length(runs), nYears = length(years), row.names = NULL)
}

#' Climate response of every ELF, and pooled
#'
#' @param elfs List of ELF inputs ([runWave1()]). @return The rows of [scoreClimateResponse()] of every ELF with a covariate, then per metric a pooled row (`level` "pooled", `ELF` "all") from [poolRatios()]:
#'   `ratio` the median over ELFs of the ELF ratios, `elfLo`, `elfHi` their range, `nELF` the ELFs, `nNonFinite` those left out as not finite. `NULL` when no ELF has a covariate.
climateTable <- function(elfs) {
  x <- do.call(rbind, lapply(elfs, scoreClimateResponse))
  if (is.null(x)) return(NULL)
  pooled <- do.call(rbind, lapply(names(climateMetrics), function(m) {
    pr <- poolRatios(x[x$metric == m, ])
    data.frame(level = "pooled", ELF = "all", metric = m, observed = NA_real_, simMedian = NA_real_, ratio = pr$ratio, runLo = NA_real_, runHi = NA_real_, obsRank = pr$obsRank,
               elfLo = pr$elfLo, elfHi = pr$elfHi, nELF = pr$nELF, nNonFinite = pr$nNonFinite, nRuns = NA_integer_, nYears = pr$nYears)
  }))
  out <- rbind(x, pooled); rownames(out) <- NULL
  out
}

## `climate` as rows of the long `results` table of runWave1(): test "C", the metric in `null`, measure "ratio"; `lo` and `hi` are the single-run range (ELF rows) or the range over ELFs (pooled rows)
climateResultRows <- function(climate) {
  if (is.null(climate)) return(NULL)
  p <- climate$level == "pooled"
  longRows(climate, "C", climate$metric, "ratio", climate$ratio, ifelse(p, climate$elfLo, climate$runLo), ifelse(p, climate$elfHi, climate$runHi))
}

## the rows of the page's climate table: the median of the ELFs, then one per ELF; each metric is a ratio cell, `obsSlope` the observed slope
climatePageRows <- function(climate) {
  if (is.null(climate)) return(list())
  rowOf <- function(label, x) {
    r <- list(elf = label)
    for (m in names(climateMetrics)) { z <- x[x$metric == m, ]; r[[m]] <- rangeCell(z, z$level == "pooled") }
    r$obsSlope <- signif(x$observed[x$metric == "slope"], 3)
    r
  }
  c(list(rowOf("all: median of ELFs", climate[climate$level == "pooled", ])),
    lapply(unique(climate$ELF[climate$level == "ELF"]), function(id) rowOf(id, climate[climate$level == "ELF" & climate$ELF == id, ])))
}
