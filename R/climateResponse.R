## Climate response: does simulated annual area burned respond to the year's climate as the observed does?
## The fits' spread climate covariate (`elf$covariate`: the ELF median of the year's value, higher = drier, e.g. CMD) is standardised over the scored years.
## Definitions are those of the climate comparison of the extreme-fire-years analysis (climCompare.R); a run is one simulated history of the scored years.

climateMetrics <- c(slope = "climate slope", dry10 = "driest 10% of years", wetHalf = "wetter half of years")

#' Climate response of one ELF
#'
#' Three ratios of the simulated to the observed record over the scored years, each computed for every run alone:
#' \itemize{
#'   \item `slope`: the slope of log annual area burned (`area50`, years with fire only) on the standardised climate covariate; sim / obs.
#'   \item `dry10`: the area burned in the driest 10% of the scored years (covariate at or above its 90% quantile); sim / obs.
#'   \item `wetHalf`: the area burned in the wetter half of the scored years (covariate at or below its median); sim / obs.
#' }
#' The area ratios use [zeroSafeRatio()]. A slope is NA with fewer than 3 years with fire, and a ratio of slopes is NA when the observed slope is 0 or NA.
#'
#' @param elf An ELF input list ([runWave1()]); without a `covariate` of the scored years there is no result (`NULL`).
#' @return `data.frame`, one row per metric: `level` ("ELF"), `ELF`, `metric`, `observed` (the observed slope, or area in ha), `simMedian` (median over runs of the run's value), `ratio` (median over runs of the
#'   run's ratio to the observed), `runLo`, `runHi` (2.5% and 97.5% of the single runs' ratio), `elfLo`, `elfHi` (NA here), `nELF` (1), `nRuns`, `nYears`.
scoreClimateResponse <- function(elf) {
  years <- sort(unique(elf$simYears$year))
  cv <- elf$covariate$value[match(years, elf$covariate$year)]
  if (is.null(elf$covariate) || anyNA(cv) || length(years) < 3 || stats::sd(cv) == 0) return(NULL)
  z <- (cv - mean(cv)) / stats::sd(cv)
  dry <- z >= stats::quantile(z, 0.9); wet <- z <= stats::median(z)
  slopeOf <- function(area) { k <- area > 0; if (sum(k) < 3) NA_real_ else unname(stats::coef(stats::lm(log(area[k]) ~ z[k]))[2]) }
  values <- function(area) c(slope = slopeOf(area), dry10 = sum(area[dry]), wetHalf = sum(area[wet]))   # area: per scored year, in `years` order
  obs <- values(elf$obsYears$area50[match(years, elf$obsYears$year)])
  ratioTo <- function(v) c(slope = if (isTRUE(obs[["slope"]] != 0)) v[["slope"]] / obs[["slope"]] else NA_real_,
                           dry10 = zeroSafeRatio(v[["dry10"]], obs[["dry10"]]), wetHalf = zeroSafeRatio(v[["wetHalf"]], obs[["wetHalf"]]))
  runs <- lapply(split(elf$simYears, elf$simYears$run), function(r) values(r$area50[match(years, r$year)]))
  v <- do.call(rbind, runs); rt <- t(vapply(runs, ratioTo, numeric(3)))
  qs <- function(x, p) { x <- x[is.finite(x)]; if (length(x)) unname(stats::quantile(x, p)) else NA_real_ }
  data.frame(level = "ELF", ELF = elf$elf, metric = names(climateMetrics), observed = unname(obs), simMedian = apply(v, 2, stats::median, na.rm = TRUE),
             ratio = apply(rt, 2, qs, 0.5), runLo = apply(rt, 2, qs, 0.025), runHi = apply(rt, 2, qs, 0.975),
             elfLo = NA_real_, elfHi = NA_real_, nELF = 1L, nRuns = length(runs), nYears = length(years), row.names = NULL)
}

#' Climate response of every ELF, and pooled
#'
#' @param elfs List of ELF inputs ([runWave1()]). @return The rows of [scoreClimateResponse()] of every ELF with a covariate, then per metric a pooled row (`level` "pooled", `ELF` "all"):
#'   `ratio` the median over ELFs of the ELF ratios, `elfLo`, `elfHi` their range, `nELF` the ELFs with a finite ratio. `NULL` when no ELF has a covariate.
climateTable <- function(elfs) {
  x <- do.call(rbind, lapply(elfs, scoreClimateResponse))
  if (is.null(x)) return(NULL)
  pooled <- do.call(rbind, lapply(names(climateMetrics), function(m) {
    r <- x$ratio[x$metric == m & is.finite(x$ratio)]
    data.frame(level = "pooled", ELF = "all", metric = m, observed = NA_real_, simMedian = NA_real_, ratio = if (length(r)) stats::median(r) else NA_real_,
               runLo = NA_real_, runHi = NA_real_, elfLo = if (length(r)) min(r) else NA_real_, elfHi = if (length(r)) max(r) else NA_real_,
               nELF = length(r), nRuns = NA_integer_, nYears = sum(x$nYears[x$metric == m]))
  }))
  out <- rbind(x, pooled); rownames(out) <- NULL
  out
}

## `climate` as rows of the long `results` table of runWave1(): test "C", the metric in `null`, measure "ratio"; `lo` and `hi` are the single-run range (ELF rows) or the range over ELFs (pooled rows)
climateResultRows <- function(climate) {
  if (is.null(climate)) return(NULL)
  p <- climate$level == "pooled"
  data.frame(test = "C", null = climate$metric, level = climate$level, ELF = climate$ELF, measure = "ratio", value = climate$ratio,
             lo = ifelse(p, climate$elfLo, climate$runLo), hi = ifelse(p, climate$elfHi, climate$runHi), nELF = climate$nELF,
             nUnits = ifelse(p, climate$nYears, climate$nRuns), verdict = NA_character_)
}

## the rows of the page's climate table: the median of the ELFs, then one per ELF; each metric is `c(ratio, low, high, NA)` (ratioCell), `obsSlope` the observed slope
climatePageRows <- function(climate) {
  if (is.null(climate)) return(list())
  rowOf <- function(label, x) {
    r <- list(elf = label)
    for (m in names(climateMetrics)) { z <- x[x$metric == m, ]; p <- z$level == "pooled"; r[[m]] <- ratioCell(z$ratio, if (p) z$elfLo else z$runLo, if (p) z$elfHi else z$runHi, NA) }
    r$obsSlope <- signif(x$observed[x$metric == "slope"], 3)
    r
  }
  c(list(rowOf("all: median of ELFs", climate[climate$level == "pooled", ])),
    lapply(unique(climate$ELF[climate$level == "ELF"]), function(id) rowOf(id, climate[climate$level == "ELF" & climate$ELF == id, ])))
}
