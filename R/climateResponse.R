## Climate response: does simulated annual area burned respond to the year's climate as the observed does?
## The fits' spread climate covariate (`elf$covariate`: the ELF median of the year's value, higher = drier, e.g. CMD) is standardised over the scored years.
## Definitions are those of the climate comparison of the extreme-fire-years analysis (climCompare.R); a run is one simulated history of the scored years.

climateMetrics <- c(slope = "climate slope", dry10 = "driest 10% of years", wetHalf = "wetter half of years")

## the scored years and the covariate standardised over them (`list(years, z)`), or NULL when the ELF has no usable covariate
climateInputs <- function(elf) {
  years <- scoredYears(elf)
  cv <- covariateOf(elf, years)
  if (is.null(elf$covariate) || anyNA(cv) || length(years) < 3 || stats::sd(cv) == 0) return(NULL)
  list(years = years, z = (cv - mean(cv)) / stats::sd(cv))
}

## the slope of `y` (a log value; years with a non-finite value are left out) on `z`: NA with fewer than 3 years
logSlope <- function(y, z) { k <- is.finite(y); if (sum(k) < 3) NA_real_ else unname(stats::coef(stats::lm(y[k] ~ z[k]))[2]) }

#' Climate response of one ELF
#'
#' Three comparisons of the simulated with the observed record over the scored years, each computed for every run alone. Area is `area50` (fires of the escape size or
#' more, as the area row of the page); a year a run (or the record) lacks counts as no area. The covariate (`elf$covariate`, higher = drier) is standardised over the scored years,
#' so the observed slope is per standard deviation of the covariate over the scored years.
#' \itemize{
#'   \item `slope`: the slope of log annual area burned on the standardised covariate, years with fire only; reported as the difference sim - obs, in the units of the slope (change in log annual area per standard deviation of the covariate; 0 is the observed response). A ratio of slopes is unstable when the observed slope is near 0, so there is none.
#'   \item `dry10`: the area burned in the driest 10% of the scored years (years at or above the 90% quantile of the covariate over the scored years); sim / obs.
#'   \item `wetHalf`: the area burned in the wetter half of the scored years (covariate at or below its median); sim / obs.
#' }
#' The area ratios use [zeroSafeRatio()]. A slope is NA with fewer than 3 years with fire (the difference is then NA); an observed slope of 0 gives a finite difference. The single runs are summarised by [runSummary()]; the observed rank of the slope is the share of runs with a slope below the observed (difference below 0).
#'
#' @param elf An ELF input list ([runWave1()]); without a `covariate` of the scored years there is no result (`NULL`).
#' @return `data.frame`, one row per metric: `level` ("ELF"), `ELF`, `metric`, `observed` (the observed slope, or area in ha), `simMedian`, `simLo`, `simHi` (median, 2.5% and 97.5% over runs of the run's value), `ratio` (median over runs of the
#'   run's value relative to the observed: the difference sim - obs for the slope, the ratio for the areas), `runLo`, `runHi` (2.5% and 97.5% of the single runs' `ratio`), `obsRank` (share of runs below the observed), `elfLo`, `elfHi` (NA here), `nELF` (1), `nNonFinite` (0), `nRuns`, `nYears`.
scoreClimateResponse <- function(elf) {
  ci <- climateInputs(elf)
  if (is.null(ci)) return(NULL)
  years <- ci$years; z <- ci$z
  dry <- z >= stats::quantile(z, 0.9); wet <- z <= stats::median(z)
  values <- function(area) c(slope = logSlope(log(area), z), dry10 = sum(area[dry], na.rm = TRUE), wetHalf = sum(area[wet], na.rm = TRUE))   # area: per scored year, in `years` order
  obs <- values(obsBy(elf, "area50", years))
  ratioTo <- function(v) c(slope = v[["slope"]] - obs[["slope"]],
                           dry10 = zeroSafeRatio(v[["dry10"]], obs[["dry10"]]), wetHalf = zeroSafeRatio(v[["wetHalf"]], obs[["wetHalf"]]))
  runs <- lapply(split(elf$simYears, elf$simYears$run), function(r) values(r$area50[match(years, r$year)]))
  v <- do.call(rbind, runs); single <- vapply(runs, ratioTo, numeric(3)); rs <- runSummary(single, at = c(0, 1, 1))
  data.frame(level = "ELF", ELF = elf$elf, metric = names(climateMetrics), observed = unname(obs), simMedian = apply(v, 2, stats::median, na.rm = TRUE),
             simLo = apply(v, 2, stats::quantile, 0.025, na.rm = TRUE, names = FALSE), simHi = apply(v, 2, stats::quantile, 0.975, na.rm = TRUE, names = FALSE),
             ratio = apply(single, 1, medFinite), runLo = rs$runLo, runHi = rs$runHi, obsRank = rs$obsRank,
             elfLo = NA_real_, elfHi = NA_real_, nELF = 1L, nNonFinite = 0L, nRuns = length(runs), nYears = length(years), row.names = NULL)
}

#' Climate response of every ELF, and pooled
#'
#' @param elfs List of ELF inputs ([runWave1()]). @return The rows of [scoreClimateResponse()] of every ELF with a covariate, then per metric a pooled row (`level` "pooled", `ELF` "all") from [poolRatios()]:
#'   `ratio` the median over ELFs of the ELF values (the slope's difference, the area ratios), `elfLo`, `elfHi` their 2.5-97.5% range over ELFs, `nELF` the ELFs, `nNonFinite` those left out as not finite. `NULL` when no ELF has a covariate.
climateTable <- function(elfs) withPooled(do.call(rbind, lapply(elfs, scoreClimateResponse)))

## the rows of `x` (one per ELF and metric, as scoreClimateResponse() gives) and per metric a pooled row (`level` "pooled", `ELF` "all") from poolRatios(); NULL for NULL
withPooled <- function(x) {
  if (is.null(x)) return(NULL)
  pooled <- do.call(rbind, lapply(unique(x$metric), function(m) {
    pr <- poolRatios(x[x$metric == m, ])
    data.frame(level = "pooled", ELF = "all", metric = m, observed = NA_real_, simMedian = NA_real_, simLo = NA_real_, simHi = NA_real_, ratio = pr$ratio, runLo = NA_real_, runHi = NA_real_, obsRank = pr$obsRank,
               elfLo = pr$elfLo, elfHi = pr$elfHi, nELF = pr$nELF, nNonFinite = pr$nNonFinite, nRuns = NA_integer_, nYears = pr$nYears)
  }))
  out <- rbind(x, pooled); rownames(out) <- NULL
  out
}

## `climate` as rows of the long `results` table of runWave1(): test "C", the metric in `null`, measure "ratio" ("difference" for the slope); `lo` and `hi` are the single-run range (ELF rows) or the range over ELFs (pooled rows)
climateResultRows <- function(climate) {
  if (is.null(climate)) return(NULL)
  p <- climate$level == "pooled"
  longRows(climate, "C", climate$metric, ifelse(climate$metric == "slope", "difference", "ratio"), climate$ratio, ifelse(p, climate$elfLo, climate$runLo), ifelse(p, climate$elfHi, climate$runHi))
}

## the rows of the page's climate table: the median of the ELFs, then one per ELF; `slope` is the difference sim - obs, `simSlope` the simulated slope (median, single-run range; no rank), `obsSlope` the observed slope; the areas are ratio cells
climatePageRows <- function(climate) {
  if (is.null(climate)) return(list())
  rowOf <- function(label, x) {
    r <- list(elf = label)
    for (m in names(climateMetrics)) { z <- x[x$metric == m, ]; r[[m]] <- rangeCell(z, z$level == "pooled") }
    s <- x[x$metric == "slope", ]
    r$simSlope <- ratioCell(s$simMedian, s$simLo, s$simHi, NA_real_)       # no pooled simulated slope: the cell is empty
    r$obsSlope <- signif(s$observed, 3)
    r
  }
  c(list(rowOf("all: median of ELFs", climate[climate$level == "pooled", ])),
    lapply(unique(climate$ELF[climate$level == "ELF"]), function(id) rowOf(id, climate[climate$level == "ELF" & climate$ELF == id, ])))
}

## Which channel is the missing climate response in? log(area) = log(ignitions) + log(escapes / ignitions) + log(area / escapes), so the slope of log annual area on the standardised covariate is the sum of three slopes.
slopeSplitMetrics <- c(slopeIgnitions = "ignitions", slopeEscapeRate = "escape rate", slopeSize = "mean escaped-fire size")

## the three slopes of one record (`ignitions`, `escapes`, `area` per scored year, NA for a year it lacks). Only years with at least one ignition, one escape and some area enter: a year with no escape has no
## log(area / escapes), so it is left out of all three, and the three then add up exactly to the slope of log(area) on those years (least-squares slopes are linear in the response).
slopeParts <- function(ignitions, escapes, area, z) {
  k <- (ignitions > 0 & escapes > 0 & area > 0) %in% TRUE
  y <- function(v) ifelse(k, v, NA_real_)
  c(slopeIgnitions = logSlope(y(log(ignitions)), z), slopeEscapeRate = logSlope(y(log(escapes / ignitions)), z), slopeSize = logSlope(y(log(area / escapes)), z))
}

#' The climate slope split into its ignition, escape-rate and fire-size parts
#'
#' The slope of log annual area burned on the covariate ([scoreClimateResponse()]) is the sum of three slopes: of log ignitions, of log escape rate (escapes / ignitions) and of log mean escaped-fire size
#' (`area50` / escapes), each on the covariate standardised over the scored years. Each is computed for the observed record and for every run, over the same years: the scored years with at least one ignition,
#' one escape and some area (a year without an escape is left out of all three parts, as it is of the total slope, so the parts add up to the total slope on those years). A slope is NA with fewer than 3 such years.
#' Each part is reported as the difference sim - obs, in log units per standard deviation of the covariate; the three differences add up to the total slope difference.
#'
#' @param elf An ELF input list ([runWave1()]); without a `covariate` of the scored years there is no result (`NULL`).
#' @return `data.frame`, one row per part (`metric` `slopeIgnitions`, `slopeEscapeRate`, `slopeSize`), with the columns of [scoreClimateResponse()].
scoreSlopeSplit <- function(elf) {
  ci <- climateInputs(elf)
  if (is.null(ci)) return(NULL)
  years <- ci$years
  obs <- slopeParts(obsBy(elf, "ignitions", years), obsBy(elf, "escapes", years), obsBy(elf, "area50", years), ci$z)
  runs <- lapply(split(elf$simYears, elf$simYears$run), function(r) { i <- match(years, r$year); slopeParts(r$ignitions[i], r$escapes[i], r$area50[i], ci$z) })
  v <- do.call(rbind, runs); single <- vapply(runs, function(s) s - obs, numeric(3)); rs <- runSummary(single, at = 0)
  data.frame(level = "ELF", ELF = elf$elf, metric = names(slopeSplitMetrics), observed = unname(obs), simMedian = apply(v, 2, stats::median, na.rm = TRUE),
             simLo = apply(v, 2, stats::quantile, 0.025, na.rm = TRUE, names = FALSE), simHi = apply(v, 2, stats::quantile, 0.975, na.rm = TRUE, names = FALSE),
             ratio = apply(single, 1, medFinite), runLo = rs$runLo, runHi = rs$runHi, obsRank = rs$obsRank,
             elfLo = NA_real_, elfHi = NA_real_, nELF = 1L, nNonFinite = 0L, nRuns = length(runs), nYears = length(years), row.names = NULL)
}

#' The slope split of every ELF, and pooled
#'
#' @param elfs List of ELF inputs ([runWave1()]). @return The rows of [scoreSlopeSplit()] of every ELF with a covariate, then per part a pooled row from [withPooled()] (as [climateTable()]). `NULL` when no ELF has a covariate.
slopeSplitTable <- function(elfs) withPooled(do.call(rbind, lapply(elfs, scoreSlopeSplit)))

## `slopeSplit` as rows of the long `results` table: test "C", the part in `null`, measure "difference"
slopeSplitResultRows <- function(slopeSplit) {
  if (is.null(slopeSplit)) return(NULL)
  p <- slopeSplit$level == "pooled"
  longRows(slopeSplit, "C", slopeSplit$metric, "difference", slopeSplit$ratio, ifelse(p, slopeSplit$elfLo, slopeSplit$runLo), ifelse(p, slopeSplit$elfHi, slopeSplit$runHi))
}

## the rows of the page's slope-split table, as climatePageRows(): per part `<m>` the difference cell, `sim<m>` the simulated slope cell (no pooled one), `obs<m>` the observed slope
slopeSplitPageRows <- function(slopeSplit) {
  if (is.null(slopeSplit)) return(list())
  rowOf <- function(label, x) {
    r <- list(elf = label)
    for (m in names(slopeSplitMetrics)) {
      z <- x[x$metric == m, ]
      r[[m]] <- rangeCell(z, z$level == "pooled")
      r[[paste0("sim", m)]] <- ratioCell(z$simMedian, z$simLo, z$simHi, NA_real_)
      r[[paste0("obs", m)]] <- signif(z$observed, 3)
    }
    r
  }
  c(list(rowOf("all: median of ELFs", slopeSplit[slopeSplit$level == "pooled", ])),
    lapply(unique(slopeSplit$ELF[slopeSplit$level == "ELF"]), function(id) rowOf(id, slopeSplit[slopeSplit$level == "ELF" & slopeSplit$ELF == id, ])))
}
