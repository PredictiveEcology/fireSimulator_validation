## Null forecasts. Each returns one ensemble per target year (a list of numeric vectors), fitted in-sample on
## the years given (the scored year included), so a model has to beat a null that has seen the answer's year.
## N0: the model's own ensembles assigned to the wrong years (a floor). N1: the ELF's own history, all years.
## N2: a one-line regression on a covariate.

#' N0: shuffled years
#'
#' Scores the model's ensembles against other years' observations. Returns, for each year, the mean loss over
#' `nShuffle` random assignments of ensembles to years (a derangement is not forced; with few years some shuffles
#' leave a year in place, which only raises the floor slightly).
#'
#' @param simEns List (one element per year) of ensembles. @param obs Observed value per year.
#' @param lossFun `function(ensemble, y)`. @param nShuffle Number of shuffles.
#' @return Mean loss per year.
nullShuffled <- function(simEns, obs, lossFun = crpsEnsemble, nShuffle = 50L, seed = 1L) {
  n <- length(obs)
  withSeed(seed, {
    perms <- replicate(nShuffle, sample.int(n), simplify = FALSE)
    vapply(seq_len(n), function(i) mean(vapply(perms, function(p) lossFun(simEns[[p[i]]], obs[i]), 0)), 0)
  })
}

withSeed <- function(seed, expr) {
  if (exists(".Random.seed", globalenv())) { old <- get(".Random.seed", globalenv()); on.exit(assign(".Random.seed", old, globalenv())) }
  else on.exit(if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv()))
  set.seed(seed); expr
}

#' N1: the ELF's own history, every year
#'
#' @param history Observed values (already transformed) of every year of the record, scored years included.
#' @param nTarget Number of target years.
#' @return List of `nTarget` copies of `history` (the empirical distribution is the ensemble).
nullHistory <- function(history, nTarget) rep(list(history[is.finite(history)]), nTarget)

#' N2 for counts: negative binomial GLM on a covariate
#'
#' Falls back to a Poisson GLM when `MASS::glm.nb` fails or the data are not over-dispersed.
#'
#' @param y Counts of the record years. @param x Covariate of the record years. @param xNew Covariate of the target years.
#' @param K Ensemble size. @param trans Transform applied to the ensemble members (default `log1p`).
#' @return List of ensembles, one per element of `xNew`.
nullNegBinGLM <- function(y, x, xNew, K = 100L, trans = log1p) {
  ok <- is.finite(y) & is.finite(x)
  d <- data.frame(y = y[ok], x = x[ok])
  fit <- tryCatch(suppressWarnings(MASS::glm.nb(y ~ x, data = d)), error = function(e) NULL)
  theta <- if (is.null(fit)) Inf else fit$theta
  if (is.null(fit) || !is.finite(theta) || theta > 1e6) {
    fit <- stats::glm(y ~ x, family = stats::poisson(), data = d); theta <- Inf
  }
  mu <- stats::predict(fit, data.frame(x = xNew), type = "response")
  lapply(mu, function(m) trans(quantileEnsemble(function(p) if (is.finite(theta)) stats::qnbinom(p, size = theta, mu = m)
                                                             else stats::qpois(p, m), K)))
}

#' N2 for annual area: linear regression of log(1 + area) on a covariate
#'
#' @param y Areas of the record years (untransformed). @param x Covariate of the record years. @param xNew,K As [nullNegBinGLM()].
#' @return List of ensembles on the log(1 + area) scale: fitted value plus normal quantiles with the residual sd.
nullLogLm <- function(y, x, xNew, K = 100L) {
  ok <- is.finite(y) & is.finite(x)
  fit <- stats::lm(log1p(y) ~ x, data = data.frame(y = y[ok], x = x[ok]))
  s <- suppressWarnings(summary(fit)$sigma)
  mu <- stats::predict(fit, data.frame(x = xNew))
  lapply(mu, function(m) m + s * stats::qnorm((seq_len(K) - 0.5) / K))
}

## beta-binomial pmf and the rate ensemble k / n at evenly spaced quantiles of a discrete distribution
bbLogPmf <- function(k, n, a, b) lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b)
rateEnsemble <- function(pmf, n, K) {
  cdf <- pmin(cumsum(pmf) / sum(pmf), 1); cdf[length(cdf)] <- 1
  (0:n)[findInterval((seq_len(K) - 0.5) / K, cdf, left.open = TRUE) + 1L] / n
}

#' N1 for the escape rate: beta-binomial fitted to the historical escapes and ignitions
#'
#' Maximum likelihood of the beta-binomial over the record years with ignitions. A rate is forecast for each
#' target year given its (observed) number of ignitions.
#'
#' @param k,n Escapes and ignitions of the record years. @param nNew Ignitions of the target years (> 0). @param K Ensemble size.
#' @return List of rate ensembles (escapes / ignitions), one per element of `nNew`.
nullBetaBinomial <- function(k, n, nNew, K = 100L) {
  ok <- is.finite(k) & is.finite(n) & n > 0; k <- k[ok]; n <- n[ok]
  p0 <- (sum(k) + 0.5) / (sum(n) + 1)
  nll <- function(th) { m <- stats::plogis(th[1]); phi <- exp(th[2])     # phi = a + b
    -sum(bbLogPmf(k, n, m * phi, (1 - m) * phi)) }
  o <- stats::optim(c(stats::qlogis(p0), log(10)), nll, control = list(maxit = 500))
  m <- stats::plogis(o$par[1]); phi <- min(exp(o$par[2]), 1e6)
  lapply(nNew, function(nn) rateEnsemble(exp(bbLogPmf(0:nn, nn, m * phi, (1 - m) * phi)), nn, K))
}

#' N2 for the escape rate: logistic regression on a covariate
#'
#' @param k,n,x Escapes, ignitions and covariate of the record years. @param nNew,xNew Ignitions and covariate of the target years. @param K Ensemble size.
#' @return List of rate ensembles, one per target year.
nullLogistic <- function(k, n, x, nNew, xNew, K = 100L) {
  ok <- is.finite(k) & is.finite(n) & is.finite(x) & n > 0
  fit <- suppressWarnings(stats::glm(cbind(k, n - k) ~ x, family = stats::binomial(), data = data.frame(k = k[ok], n = n[ok], x = x[ok])))
  p <- stats::predict(fit, data.frame(x = xNew), type = "response")
  Map(function(pp, nn) rateEnsemble(stats::dbinom(0:nn, nn, pp), nn, K), p, nNew)
}
