## Pooled verdicts: one bootstrap over ELFs and, within each, over years.

#' Resample ELFs, then years within each ELF
#'
#' @param d `data.frame` with columns `ELF` and `year` and whatever `stat` needs.
#' @param stat `function(df)` returning a numeric vector; `df` has the resampled rows and an `ELF` column in which
#'   every drawn ELF has its own label (a copy drawn twice counts twice).
#' @param nBoot Resamples. @param seed Seed (the caller's RNG state is restored).
#' @return List: `est` (stat on `d`), `lo`, `hi` (2.5% and 97.5% of the resamples, per element of `stat`), `boot`
#'   (matrix of resamples, one column each).
bootstrapPooled <- function(d, stat, nBoot = 999, seed = 1) {
  g <- split(d, d$ELF)
  est <- stat(d)
  boot <- withSeed(seed, vapply(seq_len(nBoot), function(b) {
    pick <- sample(names(g), length(g), replace = TRUE)
    stat(do.call(rbind, lapply(seq_along(pick), function(i) {
      x <- g[[pick[i]]]; x <- x[sample.int(nrow(x), replace = TRUE), , drop = FALSE]; x$ELF <- paste0(pick[i], "#", i); x })))
  }, numeric(length(est))))
  boot <- matrix(boot, nrow = length(est))
  q <- apply(boot, 1, stats::quantile, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  list(est = est, lo = q[1, ], hi = q[2, ], boot = boot)
}

#' Pooled CRPS skill (or any loss skill) with its bootstrap interval
#'
#' @param d `data.frame(ELF, year, lossModel, lossNull)`, one row per ELF-year.
#' @inheritParams bootstrapPooled
#' @return One-row `data.frame(skill, lo, hi, nELF, nUnits)`.
bootstrapSkill <- function(d, nBoot = 999, seed = 1) {
  d <- d[is.finite(d$lossModel) & is.finite(d$lossNull), , drop = FALSE]
  if (!nrow(d)) return(data.frame(skill = NA_real_, lo = NA_real_, hi = NA_real_, nELF = 0L, nUnits = 0L))
  b <- bootstrapPooled(d, function(df) skillScore(df$lossModel, df$lossNull), nBoot, seed)
  data.frame(skill = b$est, lo = b$lo, hi = b$hi, nELF = length(unique(d$ELF)), nUnits = nrow(d))
}

#' Pooled summary of a per-unit value with its bootstrap interval: the median (log score), or the share for 0/1 values (coverage)
#'
#' @param d `data.frame(ELF, year, value)`. @inheritParams bootstrapPooled
## the share of 0/1 values (coverage: the fraction of units inside a band), else the median (central tendency is
## always the median in these tests)
pooledValue <- function(x) if (all(x %in% c(0, 1))) sum(x) / length(x) else stats::median(x)

bootstrapMean <- function(d, nBoot = 999, seed = 1) {
  d <- d[is.finite(d$value), , drop = FALSE]
  if (!nrow(d)) return(data.frame(value = NA_real_, lo = NA_real_, hi = NA_real_, nELF = 0L, nUnits = 0L))
  b <- bootstrapPooled(d, function(df) pooledValue(df$value), nBoot, seed)
  data.frame(value = b$est, lo = b$lo, hi = b$hi, nELF = length(unique(d$ELF)), nUnits = nrow(d))
}

#' Pooled Spearman of the model and of an index, and the paired difference, with a bootstrap interval
#'
#' The statistic is the median over ELFs of the ELF's Spearman correlations. Model and index are resampled on the same
#' years, so the difference is paired.
#'
#' @param d `data.frame(ELF, year, obs, simMedian, ref)`. @inheritParams bootstrapPooled
#' @return `data.frame(measure = c("rhoModel", "rhoIndex", "rhoDiff"), value, lo, hi)`.
bootstrapRankDiff <- function(d, nBoot = 999, seed = 1) {
  rho <- function(a, b) if (length(a) > 2 && stats::sd(a) > 0 && stats::sd(b) > 0) stats::cor(a, b, method = "spearman") else NA_real_
  stat <- function(df) {
    r <- do.call(rbind, lapply(split(df, df$ELF), function(x) c(rho(x$obs, x$simMedian), rho(x$obs, x$ref))))
    c(stats::median(r[, 1], na.rm = TRUE), stats::median(r[, 2], na.rm = TRUE), stats::median(r[, 1] - r[, 2], na.rm = TRUE))
  }
  b <- bootstrapPooled(d[is.finite(d$ref), , drop = FALSE], stat, nBoot, seed)
  data.frame(measure = c("rhoModel", "rhoIndex", "rhoDiff"), value = b$est, lo = b$lo, hi = b$hi)
}
