## Fire-size distribution: Anderson-Darling distance of simulated vs observed sizes.

#' Anderson-Darling statistic of simulated vs observed fire sizes, against size nulls
#'
#' Per scored year: the AD statistic (`adFun`, fireSenseUtils' `adStatistic()`, the one the spread objective
#' minimises) between the pooled simulated sizes and the observed sizes, both of fires >= `minHa`. Each null provides
#' `nNull` draws of the AD statistic against the same observation, with as many fires as the simulation produced:
#' `lognormal`, a left-truncated log-normal fitted to the sizes of ALL years (no year information);
#' `history`, a resample, with replacement, of the observed sizes of the years in `hist`. With `leaveOut = TRUE`
#' (the prototype's behaviour) the scored year is left out of that pool; with `leaveOut = FALSE` it stays in (the
#' in-sample null of wave 1).
#' Lower AD is better. `pctile` is the share of null draws with AD at or below the model's (small = the model
#' beats the null); `skill` = 1 - AD_model / median AD_null (> 0 = better than the typical null draw).
#'
#' @param sim `data.table(year, ha)`: simulated fire sizes, all reps pooled.
#' @param obs `data.table(year, ha)`: observed fire sizes.
#' @param hist `data.table(year, ha)`: the historical observed fire sizes of every year.
#' @param minHa Smallest size scored (escape size).
#' @param nNull Null draws per year.
#' @param nulls Which nulls to draw: `"lognormal"`, `"history"` or both.
#' @param leaveOut Leave the scored year out of the `history` pool?
#' @param adFun Function of two samples returning the AD statistic.
#' @return A list. `byYear`: `data.table(year, nSim, nObs, adModel, null, adNullMedian, pctile, skill)`, one
#'   row per year and null (NA when the year has no observed or no simulated fire >= `minHa`). `pooled`: the same per
#'   null for the sum of AD over the years scored. `draws`: the null draws behind `pooled`.
scoreSizesVsNull <- function(sim, obs, hist, minHa = 50, nNull = 999, seed = 1,
                             nulls = c("lognormal", "history"), leaveOut = TRUE,
                             adFun = fireSenseUtils:::adStatistic) {
  nulls <- match.arg(nulls, several.ok = TRUE)
  set.seed(seed)
  sim <- sim[sim$ha >= minHa, ]; obs <- obs[obs$ha >= minHa, ]; hist <- hist[hist$ha >= minHa, ]
  if ("lognormal" %in% nulls) {
    lh <- log(hist$ha)
    nll <- function(p) { if (p[2] <= 0) return(Inf)
      -sum(stats::dlnorm(hist$ha, p[1], p[2], log = TRUE) - stats::plnorm(minHa, p[1], p[2], lower.tail = FALSE, log.p = TRUE)) }
    fit <- stats::optim(c(mean(lh), max(stats::sd(lh), 0.5)), nll)$par
    rLogn <- function(n) stats::qlnorm(stats::runif(n, stats::plnorm(minHa, fit[1], fit[2]), 1), fit[1], fit[2])
  }
  years <- sort(unique(sim$year))
  one <- lapply(years, function(y) {
    s <- sim$ha[sim$year == y]; o <- obs$ha[obs$year == y]
    ho <- if (leaveOut) hist$ha[hist$year != y] else hist$ha
    if (!length(s) || !length(o)) return(list(ad = NA_real_, nSim = length(s), nObs = length(o),
                                              draws = stats::setNames(lapply(nulls, function(nm) rep(NA_real_, nNull)), nulls)))
    draws <- lapply(stats::setNames(nulls, nulls), function(nm) switch(nm,
      lognormal = replicate(nNull, adFun(rLogn(length(s)), o)),
      history = replicate(nNull, adFun(ho[sample.int(length(ho), length(s), TRUE)], o))))
    list(ad = adFun(s, o), nSim = length(s), nObs = length(o), draws = draws)
  })
  summ <- function(ad, d) data.table::data.table(adModel = ad, adNullMedian = stats::median(d),
    pctile = if (is.na(ad)) NA_real_ else mean(d <= ad), skill = 1 - ad / stats::median(d))
  byYear <- data.table::rbindlist(lapply(seq_along(years), function(i) data.table::rbindlist(lapply(nulls, function(nm)
    cbind(data.table::data.table(year = years[i], nSim = one[[i]]$nSim, nObs = one[[i]]$nObs, null = nm),
          summ(one[[i]]$ad, one[[i]]$draws[[nm]]))))))
  ok <- !vapply(one, function(z) is.na(z$ad), TRUE)
  draws <- lapply(stats::setNames(nulls, nulls), function(nm) rowSums(vapply(one[ok], function(z) z$draws[[nm]], numeric(nNull))))
  adSum <- sum(vapply(one[ok], function(z) z$ad, 0))
  pooled <- data.table::rbindlist(lapply(nulls, function(nm)
    cbind(data.table::data.table(null = nm, nYears = sum(ok)), summ(adSum, draws[[nm]]))))
  list(byYear = byYear, pooled = pooled, draws = draws)
}

#' Median AD of the model's sizes against other years' observations (the shuffled-years floor)
#'
#' @param sim,obs As in [scoreSizesVsNull()]. @param minHa,adFun As there.
#' @return `data.table(year, adShuffled)`: for each scored year with observed and simulated fires, the median AD between
#'   the observation and the simulated sizes of every other scored year that has simulated fires.
shuffledSizesAD <- function(sim, obs, minHa = 50, adFun = fireSenseUtils:::adStatistic) {
  sim <- sim[sim$ha >= minHa, ]; obs <- obs[obs$ha >= minHa, ]
  years <- sort(unique(sim$year))
  data.table::rbindlist(lapply(years, function(y) {
    o <- obs$ha[obs$year == y]
    other <- setdiff(years, y)
    ad <- if (length(o)) vapply(other, function(z) adFun(sim$ha[sim$year == z], o), 0) else NA_real_
    data.table::data.table(year = y, adShuffled = stats::median(ad))
  }))
}
