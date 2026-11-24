## Analysis of variance of hindcast scores: which tested factors (modelling choices, ELF) move the scores.
## Model-agnostic: it sees score folders ([runWave1()] outputs, one per arm) and a table of factors per arm, nothing else.

#' Per-ELF scores of several arms as one long table, with the arms' factors
#'
#' Reads `wave1_results.csv` of each score folder (written by `mode = "score"`; it is the `results` of `wave1.rds`).
#' A metric is `<test>_<null>_<measure>` of the rows of level `"ELF"`, e.g. `S1_N1_skill`. There is one value per arm x ELF
#' x metric: the score folders hold per-year losses (`perELF` of `wave1.rds`) and run ranges of ratios (`gates`), not scores of single runs.
#' Folders with the same arm name (the ELFs of one arm scored in several folders) are combined; an ELF scored twice in one arm stops the run.
#'
#' @param scores Character, score folders. The arm of a folder is its name in this vector, else its basename.
#' @param factors `data.frame` or csv path: `arm` and one column per factor.
#' @param metrics Character, the metrics kept. NULL: every `*_N1_skill` present (skill against each ELF's own history).
#' @return `data.frame(arm, ELF, metric, value, <factors>)`; non-finite values are dropped.
readAnovaScores <- function(scores, factors, metrics = NULL) {
  if (is.character(factors)) factors <- utils::read.csv(factors, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"arm" %in% names(factors)) stop("the factor table has no column `arm`")
  if (anyDuplicated(factors$arm)) stop("the factor table has an arm twice: ", paste(unique(factors$arm[duplicated(factors$arm)]), collapse = ", "))
  arms <- if (is.null(names(scores))) basename(scores) else ifelse(nzchar(names(scores)), names(scores), basename(scores))
  if (length(none <- setdiff(arms, factors$arm))) stop("no row in the factor table for arm ", paste(none, collapse = ", "))
  long <- do.call(rbind, Map(function(f, arm) {
    r <- utils::read.csv(file.path(f, "wave1_results.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    r <- r[r$level == "ELF" & r$measure == "skill", , drop = FALSE]
    data.frame(arm = arm, ELF = as.character(r$ELF), metric = paste(r$test, r$null, r$measure, sep = "_"), value = r$value, stringsAsFactors = FALSE)
  }, scores, arms))
  if (is.null(metrics)) metrics <- grep("_N1_skill$", unique(long$metric), value = TRUE)
  if (length(none <- setdiff(metrics, long$metric))) stop("no ELF row in the scores for metric ", paste(none, collapse = ", "))
  long <- long[long$metric %in% metrics & is.finite(long$value), , drop = FALSE]
  if (anyDuplicated(long[c("arm", "ELF", "metric")])) stop("an ELF is scored twice in one arm: ", paste(unique(long$arm[duplicated(long[c("arm", "ELF", "metric")])]), collapse = ", "))
  merge(long, factors, by = "arm", sort = FALSE)
}

#' Partition the variation of one response among ELF, factors and their interactions
#'
#' Model `value ~ ELF + factors + ELF:factor + factor:factor` (every pair of factors) with sum-to-zero contrasts, terms in that order.
#' Per term: degrees of freedom, `pctSS` (sequential sum of squares, in that order, as a percentage of the total about the mean),
#' `dropOnePctSS` (the increase in residual sum of squares when the term's columns alone are left out of the full model, as a percentage of the total:
#' the share the term holds apart from the others, for unbalanced designs), `F` and `p` (sequential, against the residual).
#' A term whose columns are a combination of earlier terms' (a factor constant over the data, two factors that vary together, an
#' interaction with cells that no arm holds) has fewer degrees of freedom than columns; `note` says which earlier terms it is aliased with. It is never dropped silently.
#' With one value per arm x ELF the residual is what the terms leave, the interactions not in the model, not run-to-run noise.
#'
#' @param d `data.frame` with `value`, `ELF` and one column per factor.
#' @param factors Character, the factor columns of `d`. A factor with a single level is reported as constant and left out.
#' @return `data.frame(term, df, pctSS, dropOnePctSS, F, p, note)`, one row per term and a `Residual` row.
anovaTerms <- function(d, factors) {
  d <- d[is.finite(d$value), , drop = FALSE]
  fx <- c("ELF", factors)
  one <- fx[vapply(fx, function(f) length(unique(d[[f]])) < 2L, TRUE)]
  fx <- setdiff(fx, one)
  for (f in fx) d[[f]] <- factor(d[[f]])
  f2 <- setdiff(fx, "ELF")
  labs <- c(fx, if ("ELF" %in% fx) paste0("ELF:", f2), if (length(f2) > 1L) apply(utils::combn(f2, 2), 2, paste, collapse = ":"))
  X <- stats::model.matrix(stats::reformulate(labs), d, contrasts.arg = lapply(d[fx], function(x) "contr.sum"))
  tl <- c("(Intercept)", labs)[attr(X, "assign") + 1L]   # the term of each column; the order of `labs` is that of the columns
  y <- d$value; tot <- sum((y - mean(y))^2)
  fit <- stats::lm.fit(X, y); rss <- sum(fit$residuals^2); resDf <- length(y) - fit$rank
  eff <- fit$effects[seq_len(fit$rank)]^2; used <- fit$qr$pivot[seq_len(fit$rank)]
  aliased <- setdiff(seq_len(ncol(X)), used)
  dropOne <- function(t) { k <- tl != t; sum(stats::lm.fit(X[, k, drop = FALSE], y)$residuals^2) - rss }
  rows <- do.call(rbind, lapply(labs, function(t) {
    ss <- sum(eff[tl[used] == t]); df <- sum(tl[used] == t); a <- aliased[tl[aliased] == t]
    with <- unique(unlist(lapply(a, function(j) {   # earlier columns an aliased column is a combination of
      prev <- setdiff(used[used < j], 1L); b <- stats::lm.fit(X[, c(1L, prev), drop = FALSE], X[, j])$coefficients
      tl[c(1L, prev)][!is.na(b) & abs(b) > 1e-8 & c(1L, prev) != 1L & tl[c(1L, prev)] != t]
    })))
    note <- if (length(a)) paste0(length(a), " of ", sum(tl == t), " columns aliased", if (length(with)) paste0(" with ", paste(with, collapse = ", ")) else " with the intercept",
                                  if (!df) ": term not estimable")
    d1 <- dropOne(t)
    if (df && d1 < 1e-8 * tot && ss > 1e-8 * tot)
      note <- paste(c(note, "its sum of squares is shared with aliased terms: nothing is lost by dropping it from the full model"), collapse = "; ")
    Fv <- if (df && resDf > 0 && rss > 0) (ss / df) / (rss / resDf)
    data.frame(term = t, df = df, pctSS = 100 * ss / tot, dropOnePctSS = 100 * d1 / tot, F = if (is.null(Fv)) NA_real_ else Fv,
               p = if (is.null(Fv)) NA_real_ else stats::pf(Fv, df, resDf, lower.tail = FALSE), note = if (is.null(note)) NA_character_ else note)
  }))
  res <- data.frame(term = "Residual", df = resDf, pctSS = 100 * rss / tot, dropOnePctSS = NA_real_, F = NA_real_, p = NA_real_,
                    note = "one value per arm x ELF: the residual is the interactions left out of the model, not run-to-run noise")
  rbind(if (length(one)) data.frame(term = one, df = 0L, pctSS = 0, dropOnePctSS = 0, F = NA_real_, p = NA_real_, note = "constant: one level in the data, left out"),
        rows, res)
}

#' Cell medians of the scores
#'
#' The median (never the mean) of the arms of each metric x ELF x combination of factor levels, and the number of arms in it.
#' @param long [readAnovaScores()] table. @param factors Character, the factor columns.
#' @return `data.frame(metric, ELF, <factors>, n, median)`.
anovaCells <- function(long, factors) {
  by <- long[c("metric", "ELF", factors)]
  x <- stats::aggregate(long["value"], by, function(v) c(n = length(v), median = stats::median(v)))
  cbind(x[names(by)], n = x$value[, "n"], median = x$value[, "median"])
}

#' Stacked percent of the sum of squares per metric
#' @param terms `terms` of [computeAnova()]: `metric`, `term`, `pctSS`. @param baseSize Base font size of the theme.
#' @return A ggplot.
plotAnova <- function(terms, baseSize = 15) {
  lv <- unique(terms$term)
  ggplot(terms, aes(x = metric, y = pctSS, fill = factor(term, lv))) + geom_col(width = 0.7) +
    scale_fill_manual(values = stats::setNames(rep_len(okabeIto, length(lv)), lv), name = NULL) +
    labs(x = NULL, y = "% of total sum of squares (sequential)") + theme_bw(base_size = baseSize)
}

#' Partition the variation of each metric of the scores
#'
#' @param long [readAnovaScores()] table.
#' @return List: `terms` (`metric`, then [anovaTerms()]), `cells` ([anovaCells()]), `long` (the values analysed).
computeAnova <- function(long) {
  fac <- setdiff(names(long), c("arm", "ELF", "metric", "value"))
  terms <- do.call(rbind, lapply(unique(long$metric), function(m) cbind(metric = m, anovaTerms(long[long$metric == m, ], fac))))
  rownames(terms) <- NULL
  list(terms = terms, cells = anovaCells(long, fac), long = long)
}

#' Write the two tables of [computeAnova()]
#'
#' @param res A [computeAnova()] result. @param outDir Folder of the two csv files.
#' @return `res`, invisibly. Written: `anova_terms.csv`, `anova_cells.csv`.
writeAnova <- function(res, outDir) {
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res$terms, file.path(outDir, "anova_terms.csv"), row.names = FALSE)
  utils::write.csv(res$cells, file.path(outDir, "anova_cells.csv"), row.names = FALSE)
  invisible(res)
}
