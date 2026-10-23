## Skill tests S1-S5. Each takes one ELF input (see runWave1()) and returns per-year tables, never a verdict:
## verdicts are pooled over ELFs by runWave1(). Losses are CRPS on log(1 + x) unless stated.

MINHA <- 50                                   # escape size: size and area metrics use fires of this size or more

## losses of the model and of N0 / N1 / N2 for each scored year, plus the model's calibration
## simEns: list (per scored year) of simulated ensembles; obs: observed per year; nulls: named list (N1, N2) of lists of ensembles
scoreEnsembles <- function(years, simEns, obs, nulls, nShuffle = 50L, lossFun = crpsEnsemble) {
  lossModel <- vapply(seq_along(years), function(i) lossFun(simEns[[i]], obs[i]), 0)
  n0 <- nullShuffled(simEns, obs, lossFun, nShuffle)
  scores <- rbind(data.frame(year = years, null = "N0", lossModel = lossModel, lossNull = n0),
                  do.call(rbind, lapply(names(nulls), function(nm)
                    data.frame(year = years, null = nm, lossModel = lossModel,
                               lossNull = vapply(seq_along(years), function(i) lossFun(nulls[[nm]][[i]], obs[i]), 0)))))
  cal <- data.frame(year = years, obs = obs,
                    cover50 = vapply(seq_along(years), function(i) insideInterval(simEns[[i]], obs[i], 0.5), TRUE),
                    cover90 = vapply(seq_along(years), function(i) insideInterval(simEns[[i]], obs[i], 0.9), TRUE))
  list(scores = scores, cal = cal)
}

covariateOf <- function(elf, years) elf$covariate$value[match(years, elf$covariate$year)]
scoredYears <- function(elf) sort(unique(elf$simYears$year))
simBy <- function(elf, col, years) lapply(years, function(y) elf$simYears[[col]][elf$simYears$year == y])
obsBy <- function(elf, col, years) elf$obsYears[[col]][match(years, elf$obsYears$year)]

#' S1: ignitions per year
#'
#' CRPS on log(1 + ignitions). N1: the ELF's own history, all years. N2: negative binomial GLM on the ELF-mean covariate.
#' Calibration: coverage of the 50% and 90% intervals, and the log score of the observation under a distribution fitted to
#' the simulated ensemble (`logScore`, negative binomial).
#' @param elf ELF input list. @param K Ensemble size of the fitted nulls. @return `list(scores, cal)`.
scoreS1 <- function(elf, K = 100L, nShuffle = 50L) {
  yrs <- scoredYears(elf); h <- elf$obsYears
  sim <- simBy(elf, "ignitions", yrs); obs <- obsBy(elf, "ignitions", yrs)
  r <- scoreEnsembles(yrs, lapply(sim, log1p), log1p(obs),
                      list(N1 = nullHistory(log1p(h$ignitions), length(yrs)),
                           N2 = nullNegBinGLM(h$ignitions, covariateOf(elf, h$year), covariateOf(elf, yrs), K)), nShuffle)
  r$cal$logScore <- vapply(seq_along(yrs), function(i) fitPredictive(sim[[i]], "count")(obs[i]), 0)
  r
}

#' S2: escape rate given the observed number of ignitions
#'
#' CRPS on escapes / ignitions, for the scored years with at least one observed ignition. The simulated ensemble is the
#' escape rate of each run (runs with no ignition dropped). N1: beta-binomial fitted to the ELF's history. N2: logistic
#' regression on the covariate.
#' @inheritParams scoreS1 @return `list(scores, cal)`.
scoreS2 <- function(elf, K = 100L, nShuffle = 50L) {
  h <- elf$obsYears; yrs <- intersect(scoredYears(elf), h$year[h$ignitions > 0])
  s <- elf$simYears
  sim <- lapply(yrs, function(y) { x <- s[s$year == y & s$ignitions > 0, ]; x$escapes / x$ignitions })
  nObs <- obsBy(elf, "ignitions", yrs); obs <- obsBy(elf, "escapes", yrs) / nObs
  keep <- lengths(sim) > 0; yrs <- yrs[keep]; sim <- sim[keep]; nObs <- nObs[keep]; obs <- obs[keep]
  r <- scoreEnsembles(yrs, sim, obs,
                      list(N1 = nullBetaBinomial(h$escapes, h$ignitions, nObs, K),
                           N2 = nullLogistic(h$escapes, h$ignitions, covariateOf(elf, h$year), nObs, covariateOf(elf, yrs), K)), nShuffle)
  r
}

#' S3: annual area of fires of 50 ha or more
#'
#' CRPS on log(1 + area). N1: the ELF's own history. N2: regression of log(1 + area) on the covariate.
#' @inheritParams scoreS1 @return `list(scores, cal)`.
scoreS3 <- function(elf, K = 100L, nShuffle = 50L) {
  yrs <- scoredYears(elf); h <- elf$obsYears
  sim <- simBy(elf, "area50", yrs); obs <- obsBy(elf, "area50", yrs)
  r <- scoreEnsembles(yrs, lapply(sim, log1p), log1p(obs),
                      list(N1 = nullHistory(log1p(h$area50), length(yrs)),
                           N2 = nullLogLm(h$area50, covariateOf(elf, h$year), covariateOf(elf, yrs), K)), nShuffle)
  r$cal$logScore <- vapply(seq_along(yrs), function(i) fitPredictive(sim[[i]], "area")(obs[i]), 0)
  r
}

#' S4: fire sizes
#'
#' (a) Anderson-Darling distance between the pooled simulated and the observed sizes of fires of 50 ha or more, per
#' scored year; the loss of N1 is the median AD of draws from the sizes pooled over all years (scored year included).
#' N0 is the model's sizes from other years. (b) CRPS of the largest fire of the year on log(1 + ha), 0 when no fire; N1 is
#' the ELF's annual maxima over all years. N2 does not apply to sizes.
#' @param elf ELF input list. @param nNull Draws per year for the AD null. @param adFun The AD statistic.
#' @return `list(ad, largest)`: each `list(scores, cal)`-like; `ad$scores` has `year, null, lossModel, lossNull`.
scoreS4 <- function(elf, nNull = 199, nShuffle = 50L, adFun = fireSenseUtils:::adStatistic) {
  yrs <- scoredYears(elf)
  sim <- elf$simSizes[elf$simSizes$ha >= MINHA, ]; hist <- elf$obsSizes[elf$obsSizes$ha >= MINHA, ]
  a <- scoreSizesVsNull(sim, hist, hist, MINHA, nNull, nulls = "history", leaveOut = FALSE, adFun = adFun)$byYear
  sh <- shuffledSizesAD(sim, hist, MINHA, adFun)
  ad <- rbind(data.frame(year = a$year, null = "N1", lossModel = a$adModel, lossNull = a$adNullMedian),
              data.frame(year = sh$year, null = "N0", lossModel = a$adModel[match(sh$year, a$year)], lossNull = sh$adShuffled))
  ## (b) the largest fire per year, per run
  runs <- sort(unique(elf$simYears$run))
  maxOf <- function(sz, y, r) {
    keep <- sz$year == y
    if (!is.null(r)) keep <- keep & sz$run == r
    if (any(keep)) max(sz$ha[keep]) else 0
  }
  simMax <- lapply(yrs, function(y) vapply(runs, function(r) maxOf(sim, y, r), 0))
  obsMax <- vapply(yrs, function(y) maxOf(hist, y, NULL), 0)
  allYears <- elf$obsYears$year
  annualMax <- log1p(vapply(allYears, function(y) maxOf(hist, y, NULL), 0))
  lg <- scoreEnsembles(yrs, lapply(simMax, log1p), log1p(obsMax), list(N1 = nullHistory(annualMax, length(yrs))), nShuffle)
  list(ad = list(scores = ad), largest = lg)
}

#' S5: ranking of years
#'
#' Spearman of the model's mean annual area (fires of 50 ha or more) against the observed area, with the permutation floor
#' (N0) and the index's own Spearman (N2: the covariate against the observed area).
#' @param elf ELF input list. @param nPerm Permutations.
#' @return One-row `data.frame(nYears, rho, pPerm, rhoIndex, rhoDiff)`, and the per-year table as attribute `years`.
scoreS5 <- function(elf, nPerm = 9999) {
  yrs <- scoredYears(elf)
  d <- data.frame(year = yrs, obs = obsBy(elf, "area50", yrs),
                  simMean = vapply(simBy(elf, "area50", yrs), mean, 0), ref = covariateOf(elf, yrs))
  r <- as.data.frame(scoreRanking(data.table::as.data.table(d), elf$covariate, nPerm))
  names(r)[names(r) == "rhoRef"] <- "rhoIndex"
  r$rhoDiff <- r$rho - r$rhoIndex
  attr(r, "years") <- d
  r
}
