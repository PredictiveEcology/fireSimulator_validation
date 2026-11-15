## Weight tuning: which objective-term weights pick the spread-fit members that hindcast best.
## Inputs are the terms table (one row per member x ELF) and the bySet scores ([scoreBySet()]) of the same members.

#' Choose objective weights on some ELFs and report on others
#'
#' For each weight vector `w` of a grid fixed in advance and each ELF, the member picked is the one that minimises
#' `sum_t w_t * T(member, t)` over the ELF's members (first on ties). Its validation indicators are its row of the `bySet`
#' scores ([scoreBySet()]: medians of the ratios with single-run 2.5-97.5% ranges, skill, gate results).
#' The grid's weight vector with the lowest loss on the tune ELFs is selected, and only the selected vector's rows on the
#' report ELFs are the result (`report`). Rows of the tune ELFs are returned apart (`tuning`), labelled `set = "tune"`.
#' With `leaveOneOut`, each report ELF is held out in turn: the weights are selected on the tune ELFs without it, and the
#' held-out ELF's rows are its result. Years are held out by scoring them separately: `tuneBySet` holds scores of the tune years.
#'
#' @param terms `data.frame`: `ELF`, `member` (the set id: `paramSet` of the bySet table) and one numeric column per raw objective term.
#' @param bySet `data.frame` of [scoreBySet()] (`ELF`, `paramSet`, ...): the scores the report rows come from.
#' @param grid `data.frame`, one row per weight vector (at most 50), one column per term it weights; terms it does not name are
#'   not in the objective. An optional `weightId` column names the vectors.
#' @param tuneELFs,reportELFs ELFs the weights are selected on and reported on. NULL `tuneELFs`: the ELFs of `terms` that are not reported.
#'   `reportELFs` NULL is allowed only with `leaveOneOut` (then the tune ELFs).
#' @param leaveOneOut Hold out each report ELF in turn.
#' @param loss Column of `bySet` the choice minimises, as the median over the tune ELFs: for a `_ratio` column the distance of the ratio
#'   from 1 on the log scale, else minus the value (skill).
#' @param tuneBySet Scores the selection uses, default `bySet`. Needed to tune and report on the same ELFs (other years).
#' @return List: `report` (`fold`, `set = "report"`, `weightId`, the weights, `ELF`, `member` and the `bySet` columns),
#'   `tuning` (the same for every grid vector on the tune ELFs, `set = "tune"`), `losses` (`fold`, `weightId`, `loss`, `selected`).
tuneWeights <- function(terms, bySet, grid, tuneELFs = NULL, reportELFs = NULL, leaveOneOut = FALSE, loss = "area_ratio", tuneBySet = NULL) {
  if (nrow(grid) > 50) stop("the weight grid has ", nrow(grid), " vectors: at most 50, fixed in advance")
  if (is.null(grid$weightId)) grid$weightId <- seq_len(nrow(grid))
  tcols <- setdiff(names(grid), "weightId")
  if (length(miss <- setdiff(c("ELF", "member", tcols), names(terms)))) stop("terms table has no column ", paste(miss, collapse = ", "))
  if (!loss %in% names(bySet)) stop("loss column '", loss, "' is not in the bySet scores")
  honest <- !is.null(tuneBySet); if (is.null(tuneBySet)) tuneBySet <- bySet
  allELFs <- unique(as.character(terms$ELF))
  if (is.null(reportELFs) && !leaveOneOut) stop("reportELFs needed unless leaveOneOut")
  if (is.null(tuneELFs)) tuneELFs <- if (leaveOneOut) allELFs else setdiff(allELFs, reportELFs)
  if (is.null(reportELFs)) reportELFs <- tuneELFs
  if (!leaveOneOut && !honest && length(both <- intersect(tuneELFs, reportELFs)))
    stop("ELFs ", paste(both, collapse = ", "), " are tuned and reported: use leaveOneOut, or tuneBySet scored on other years")
  pickRows <- function(elfs, scores, set) do.call(rbind, lapply(elfs, function(e) {
    tm <- terms[as.character(terms$ELF) == e, , drop = FALSE]
    tm <- tm[tm$member %in% scores$paramSet[as.character(scores$ELF) == e], , drop = FALSE]
    if (!nrow(tm)) stop("ELF ", e, ": no member is in both the terms table and the bySet scores")
    do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
      w <- unlist(grid[i, tcols, drop = FALSE])
      m <- tm$member[which.min(as.matrix(tm[, tcols, drop = FALSE]) %*% w)]
      cbind(set = set, grid[i, , drop = FALSE],
            scores[as.character(scores$ELF) == e & scores$paramSet == m, , drop = FALSE][1, , drop = FALSE], member = m)
    }))
  }))
  lossOf <- function(x) { v <- x[[loss]]; if (grepl("_ratio$", loss)) abs(log(v)) else -v }
  tuning <- pickRows(tuneELFs, tuneBySet, "tune")
  folds <- if (leaveOneOut) reportELFs else "all"
  sel <- do.call(rbind, lapply(folds, function(f) {
    use <- if (leaveOneOut) setdiff(tuneELFs, f) else tuneELFs
    if (!length(use)) stop("no tune ELF left to select weights for ", f)
    x <- tuning[tuning$ELF %in% use, , drop = FALSE]
    l <- tapply(lossOf(x), x$weightId, stats::median)
    l <- l[as.character(grid$weightId)]
    data.frame(fold = f, weightId = grid$weightId, loss = as.numeric(l), selected = seq_along(l) == which.min(l))
  }))
  report <- do.call(rbind, lapply(folds, function(f) {
    id <- sel$weightId[sel$fold == f & sel$selected]
    r <- pickRows(if (leaveOneOut) f else reportELFs, bySet, "report")
    cbind(fold = f, r[r$weightId == id, , drop = FALSE])
  }))
  tuning <- cbind(fold = "all", tuning)
  rownames(report) <- rownames(tuning) <- NULL
  list(report = report, tuning = tuning, losses = sel)
}

#' Read the three tables of weight tuning, run [tuneWeights()] and write its tables
#'
#' @param termsFile Terms table, csv: `ELF`, `member` and one column per raw objective term (what `fireSense_spreadFit` writes as `spreadFitTerms_*.csv`).
#' @param bySetFile `wave1_bySet.csv`. @param gridFile Weight grid, csv, one row per vector. @param outDir Folder of the three result csv files.
#' @param tuneBySetFile Optional `bySet` csv of the tune years. @param ... Passed to [tuneWeights()].
#' @return The [tuneWeights()] list, invisibly. Written: `weights_report.csv` (the held-out rows: the result), `weights_tuning.csv`, `weights_losses.csv`.
runWeights <- function(termsFile, bySetFile, gridFile, outDir, tuneBySetFile = NA_character_, ...) {
  rd <- function(f) utils::read.csv(f, check.names = FALSE)
  res <- tuneWeights(rd(termsFile), rd(bySetFile), rd(gridFile), tuneBySet = if (!is.na(tuneBySetFile)) rd(tuneBySetFile), ...)
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  for (n in c("report", "tuning", "losses")) utils::write.csv(res[[n]], file.path(outDir, paste0("weights_", n, ".csv")), row.names = FALSE)
  invisible(res)
}
