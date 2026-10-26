## Scoring rules on ensembles. Plain numeric vectors in, numbers out; nothing here knows what a fire is.

#' Ensemble CRPS
#'
#' `E|X - y| - 0.5 E|X - X'|` for an ensemble `x` and one observation `y` (lower is better). Same formula as the
#' hindcast evaluation prototype (`eval_boundfix/new/stats.R`).
#'
#' @param x Numeric ensemble. @param y One observed value.
#' @return One number.
crpsEnsemble <- function(x, y) mean(abs(x - y)) - 0.5 * mean(abs(outer(x, x, "-")))

#' CRPS skill: how much of a reference's loss the model removes
#'
#' @param lossModel,lossNull Losses (CRPS or another loss), one per unit scored.
#' @return `1 - sum(lossModel) / sum(lossNull)`: 1 = perfect, 0 = as good as the null, negative = worse. `NA` when
#'   the null's loss is 0 (nothing to improve on).
skillScore <- function(lossModel, lossNull) {
  s <- sum(lossNull)
  if (!is.finite(s) || s <= 0) NA_real_ else 1 - sum(lossModel) / s
}

#' Evenly spaced quantiles of a distribution, used as a deterministic ensemble
#'
#' @param qfun A quantile function (vectorised in its first argument). @param K Ensemble size.
quantileEnsemble <- function(qfun, K = 100L) qfun((seq_len(K) - 0.5) / K)

#' Was the observation inside the central interval of the ensemble?
#'
#' @param level 0.5 or 0.9 (any value in (0, 1)); the interval runs from quantile `(1 - level) / 2` to `(1 + level) / 2`.
#' @return `TRUE`/`FALSE`.
insideInterval <- function(x, y, level) {
  q <- stats::quantile(x, c((1 - level) / 2, (1 + level) / 2), names = FALSE)
  y >= q[1] & y <= q[2]
}

#' Log predictive score of an observation under a distribution fitted to an ensemble
#'
#' Smooth, so a value outside the sample never scores -Inf. `"count"`: negative binomial by moments (Poisson when the
#' sample is not over-dispersed), mean shrunk by half a pseudo-count so a sample of zeros still gives a positive mean.
#' `"area"`: a point mass at 0 (half a pseudo-observation) plus a log-normal for the positive values; with fewer than
#' 2 positive values the log-normal is N(log of the value or 0, sd 2). The mass-plus-density mix is only meaningful
#' for differences between two fits of this family.
#'
#' @param x Numeric sample (>= 0).
#' @param type `"count"` or `"area"`.
#' @return A function of `y` returning the log-probability (count) or log-mass/density (area).
fitPredictive <- function(x, type = c("count", "area")) {
  type <- match.arg(type)
  n <- length(x)
  if (type == "count") {
    mu <- (sum(x) + 0.5) / (n + 1)
    v <- if (n > 1) stats::var(x) else mu
    size <- if (v > mu * 1.0001) mu^2 / (v - mu) else Inf
    return(function(y) if (is.finite(size)) stats::dnbinom(y, size = size, mu = mu, log = TRUE)
                       else stats::dpois(y, mu, log = TRUE))
  }
  p0 <- (sum(x == 0) + 0.5) / (n + 1)
  pos <- log(x[x > 0])
  m <- if (length(pos)) mean(pos) else 0
  s <- if (length(pos) >= 2) max(stats::sd(pos), 0.1) else 2
  function(y) ifelse(y == 0, log(p0), log1p(-p0) + stats::dlnorm(pmax(y, 1e-12), m, s, log = TRUE))
}
