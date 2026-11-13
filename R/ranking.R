## Does the simulation put the years in the right order?

#' Spearman of the model median against the observation, with a permutation floor
#'
#' Spearman correlation, across scored years, between the observed value and the median of the simulated reps, with a
#' one-sided permutation p-value (years shuffled: the N0 floor). A `ref` series (e.g. an annual fire-weather index)
#' gives the higher bar: the same correlation for `ref`.
#'
#' @param d `data.table(year, obs, simMedian)`.
#' @param ref Optional `data.table(year, value)`: a reference predictor of the same years.
#' @param nPerm Permutations.
#' @return One-row `data.table(nYears, rho, pPerm, rhoRef)`.
scoreRanking <- function(d, ref = NULL, nPerm = 9999, seed = 1) {
  set.seed(seed)
  rho <- suppressWarnings(stats::cor(d$obs, d$simMedian, method = "spearman"))
  perm <- replicate(nPerm, suppressWarnings(stats::cor(d$obs[sample.int(nrow(d))], d$simMedian, method = "spearman")))
  rhoRef <- if (is.null(ref)) NA_real_ else
    suppressWarnings(stats::cor(d$obs, ref$value[match(d$year, ref$year)], method = "spearman", use = "complete.obs"))
  data.table::data.table(nYears = nrow(d), rho = rho, pPerm = (1 + sum(perm >= rho, na.rm = TRUE)) / (nPerm + 1), rhoRef = rhoRef)
}

#' Rank arms by their area skill against N1
#'
#' The arms' `runWave1()` score tables bound together (with an `arm` column) give the ranking: the median, over the ELFs that every
#' arm has scored, of the area skill against N1 (`S1_N1_skill`; median, not mean). Arms with no skill sort last.
#'
#' @param scores The `scores` of several [runWave1()] runs bound with an `arm` column: `arm`, `level` (`"pooled"` or `"ELF"`), `ELF`, `key`, `value`, `lo`, `hi`.
#'   Common ELFs are those with a finite ELF-level `S1_N1_skill` in every arm that has any.
#' @param summary [summariseGates()] output with an `arm` column, one row per ELF x arm (for the worst area and the largest of the 10 largest fires).
#' @return One row per arm, best first: `rank`, `arm`, `medianSkill`, `commonELFs` (comma-separated), `skill_<ELF>` for every ELF in `scores`,
#'   `pooledSkill`, `pooledLo`, `pooledHi`, `worstArea` (the `area_ratio` furthest from 1 on the log scale over the arm's ELFs), `largestTop10`
#'   (the largest `top10_median_ratio`), `nELF`, `ELFs` (the arm's ELFs in `summary`, comma-separated).
rankArms <- function(scores, summary) {
  skill <- scores[scores$key == "S1_N1_skill", , drop = FALSE]
  arms <- unique(as.character(scores$arm))
  perELF <- skill[skill$level == "ELF", , drop = FALSE]
  elfs <- unique(as.character(perELF$ELF))
  val <- function(a, e) { v <- perELF$value[perELF$arm == a & as.character(perELF$ELF) == e]; if (length(v)) v[1] else NA_real_ }
  mat <- vapply(elfs, function(e) vapply(arms, val, 0, e = e), numeric(length(arms)))
  mat <- matrix(mat, nrow = length(arms), dimnames = list(arms, elfs))
  scored <- rowSums(is.finite(mat)) > 0
  common <- elfs[colSums(!is.finite(mat[scored, , drop = FALSE])) == 0 & any(scored)]
  pooled <- function(a, col) { r <- skill[skill$level == "pooled" & skill$arm == a, col]; if (length(r)) r[1] else NA_real_ }
  worst <- function(x) { x <- x[is.finite(x) & x > 0]; if (length(x)) x[which.max(abs(log(x)))] else NA_real_ }
  big <- function(x) { x <- x[is.finite(x)]; if (length(x)) max(x) else NA_real_ }
  out <- do.call(rbind, lapply(arms, function(a) {
    s <- summary[as.character(summary$arm) == a, , drop = FALSE]
    data.frame(arm = a, medianSkill = if (length(common)) stats::median(mat[a, common]) else NA_real_,
               commonELFs = paste(common, collapse = ", "),
               stats::setNames(as.data.frame(as.list(mat[a, ]), check.names = FALSE), paste0("skill_", elfs)),
               pooledSkill = pooled(a, "value"), pooledLo = pooled(a, "lo"), pooledHi = pooled(a, "hi"),
               worstArea = worst(s$area_ratio), largestTop10 = big(s$top10_median_ratio),
               nELF = nrow(s), ELFs = paste(unique(as.character(s$ELF)), collapse = ", "), stringsAsFactors = FALSE, check.names = FALSE)
  }))
  out <- out[order(-out$medianSkill, -out$pooledSkill, na.last = TRUE), , drop = FALSE]
  cbind(rank = seq_len(nrow(out)), out, row.names = NULL)
}
