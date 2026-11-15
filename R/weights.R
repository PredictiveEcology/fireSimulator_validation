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
#' @param terms `data.frame`: `ELF`, `member` (the set id: it must equal `paramSet` of the bySet table; a member with no bySet row stops the run) and one numeric
#'   column per raw objective term. A `rep` column (replicate scores of a member, as `spreadFitTerms_*.csv`) is reduced to the median per member before selection.
#'   NA or infinite terms stop the run; a NA loss is left out of the median with a warning naming the ELF and weightId.
#' @param bySet `data.frame` of [scoreBySet()] (`ELF`, `paramSet`, ...): the scores the report rows come from.
#' @param grid `data.frame`, one row per weight vector (at most 50), one column per objective term it weights (term names, never parameter names); terms it does not name are
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
  bad <- terms[!stats::complete.cases(terms[, c("ELF", "member", tcols)]) | !is.finite(rowSums(terms[, tcols, drop = FALSE])), , drop = FALSE]
  if (nrow(bad)) stop("terms are NA or not finite: ELF ", paste(unique(bad$ELF), collapse = ", "), ", member ", paste(unique(bad$member), collapse = ", "))
  if ("rep" %in% names(terms))      # replicates of one member: the median, before any selection (the minimum over replicates would pick the luckiest)
    terms <- stats::aggregate(terms[tcols], terms[c("ELF", "member")], stats::median)
  terms <- terms[order(as.character(terms$ELF), terms$member), , drop = FALSE]   # ties go to the lowest member, not to the file order
  honest <- !is.null(tuneBySet); if (is.null(tuneBySet)) tuneBySet <- bySet
  allELFs <- unique(as.character(terms$ELF))
  if (is.null(reportELFs) && !leaveOneOut) stop("reportELFs needed unless leaveOneOut")
  if (is.null(tuneELFs)) tuneELFs <- if (leaveOneOut) allELFs else setdiff(allELFs, reportELFs)
  if (is.null(reportELFs)) reportELFs <- tuneELFs
  if (!leaveOneOut && !honest && length(both <- intersect(tuneELFs, reportELFs)))
    stop("ELFs ", paste(both, collapse = ", "), " are tuned and reported: use leaveOneOut, or tuneBySet scored on other years")
  pickRows <- function(elfs, scores, set) do.call(rbind, lapply(elfs, function(e) {
    tm <- terms[as.character(terms$ELF) == e, , drop = FALSE]
    if (length(none <- setdiff(tm$member, scores$paramSet[as.character(scores$ELF) == e])))
      stop("ELF ", e, ": ", length(none), " of ", nrow(tm), " terms members have no bySet row (member = paramSet): ", paste(utils::head(none, 10), collapse = ", "))
    do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
      w <- unlist(grid[i, tcols, drop = FALSE])
      m <- tm$member[which.min(as.matrix(tm[, tcols, drop = FALSE]) %*% w)]
      cbind(set = set, grid[i, , drop = FALSE],
            scores[as.character(scores$ELF) == e & scores$paramSet == m, , drop = FALSE][1, , drop = FALSE], member = m)
    }))
  }))
  lossOf <- function(x) { v <- x[[loss]]; v <- if (grepl("_ratio$", loss)) abs(log(v)) else -v; v[!is.finite(v)] <- NA_real_; v }
  tuning <- pickRows(tuneELFs, tuneBySet, "tune")
  folds <- if (leaveOneOut) reportELFs else "all"
  sel <- do.call(rbind, lapply(folds, function(f) {
    use <- if (leaveOneOut) setdiff(tuneELFs, f) else tuneELFs
    if (!length(use)) stop("no tune ELF left to select weights for ", f)
    x <- tuning[tuning$ELF %in% use, , drop = FALSE]
    x$loss <- lossOf(x)
    if (anyNA(x$loss)) warning("loss '", loss, "' is NA or not finite, left out of the median: ELF ", paste(unique(x$ELF[is.na(x$loss)]), collapse = ", "),
                               ", weightId ", paste(unique(x$weightId[is.na(x$loss)]), collapse = ", "))
    l <- tapply(x$loss, x$weightId, stats::median, na.rm = TRUE)
    l <- l[as.character(grid$weightId)]
    if (all(is.na(l))) stop("no finite loss '", loss, "' for any weight vector (fold ", f, ")")
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

#' Read the terms files of the ELFs
#'
#' `fireSense_spreadFit` writes `spreadFitTerms_<runName>.csv` for one ELF (`member`, `rep`, parameters, terms; no `ELF` column).
#' `<runName>` is the run's output folder below `outputs/` with `/` as `_`, so it holds the ELF as a token like `5.2`.
#'
#' @param files Character vector of csv files. Named by ELF: those names are the ELFs. Unnamed: the ELF is the one token of the file name
#'   (split at `_`) that looks like `<number>.<number>`; none or several stops. A file that has an `ELF` column keeps it.
#' @return `data.frame` of all files, with an `ELF` column.
readTermsFiles <- function(files) {
  elfOf <- function(f, nm) {
    if (nzchar(nm)) return(nm)
    tok <- grep("^[0-9]+\\.[0-9]+$", strsplit(tools::file_path_sans_ext(basename(f)), "_", fixed = TRUE)[[1]], value = TRUE)
    if (length(unique(tok)) != 1L) stop("cannot tell the ELF from the file name ", basename(f), " (", length(unique(tok)), " candidates): name the files by ELF")
    tok
  }
  nms <- if (is.null(names(files))) rep("", length(files)) else names(files)
  do.call(rbind, lapply(seq_along(files), function(i) {
    d <- utils::read.csv(files[[i]], check.names = FALSE)
    if (!"ELF" %in% names(d)) d <- cbind(ELF = elfOf(files[[i]], nms[i]), d)
    d
  }))
}

#' Read the three tables of weight tuning, run [tuneWeights()] and write its tables
#'
#' @param termsFile One or several terms files ([readTermsFiles()]: one per ELF as `fireSense_spreadFit` writes them, or one table with an `ELF` column).
#' @param bySetFile `wave1_bySet.csv`. @param gridFile Weight grid, csv, one row per vector. @param outDir Folder of the three result csv files.
#' @param tuneBySetFile Optional `bySet` csv of the tune years. @param ... Passed to [tuneWeights()].
#' @return The [tuneWeights()] list, invisibly. Written: `weights_report.csv` (the held-out rows: the result), `weights_tuning.csv`, `weights_losses.csv`.
runWeights <- function(termsFile, bySetFile, gridFile, outDir, tuneBySetFile = NA_character_, ...) {
  rd <- function(f) utils::read.csv(f, check.names = FALSE)
  res <- tuneWeights(readTermsFiles(termsFile), rd(bySetFile), rd(gridFile), tuneBySet = if (!is.na(tuneBySetFile)) rd(tuneBySetFile), ...)
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  for (n in c("report", "tuning", "losses")) utils::write.csv(res[[n]], file.path(outDir, paste0("weights_", n, ".csv")), row.names = FALSE)
  invisible(res)
}
