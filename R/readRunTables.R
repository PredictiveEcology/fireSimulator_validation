## A simulator-neutral reader: per-run tables in, the ELF input list of runWave1() out.

#' One ELF's input from plain per-run tables
#'
#' For a simulator that does not write FireSense hindcast folders (a fast hindcast, another model): the sizes of the fires
#' and the ignitions per run-year are enough. A run-year without fires of `minHa` or more is 0 area and 0 escapes.
#'
#' @param elf ELF name (any string; it names the figures and the rows of the page).
#' @param simSizes `data.frame(year, run, ha)`, one row per simulated fire; fires under `minHa` are dropped.
#' @param simIgnitions `data.frame(year, run, ignitions)`, one row per run-year, ignitions in the whole ELF.
#' @param obsYears,obsSizes,covariate,fit As in [runWave1()]; `covariate` and `fit` may stay `NULL` and `list()` when the
#'   skill tests that use them are not run.
#' @param years Simulated years. @param runs Number of runs. @param minHa Smallest fire (ha) counted as an escape.
#' @return The list [runWave1()] documents.
runTablesToELF <- function(elf, simSizes, simIgnitions, obsYears, obsSizes, years = sort(unique(simIgnitions$year)),
                           runs = max(simIgnitions$run), covariate = NULL, fit = list(), minHa = 50) {
  big <- data.table::as.data.table(simSizes)[ha >= minHa, .(year = as.integer(year), run, ha)]
  sy <- data.table::CJ(year = as.integer(years), run = seq_len(runs))
  ig <- data.table::as.data.table(simIgnitions)[, .(year = as.integer(year), run, ignitions)]
  sy <- merge(sy, ig, by = c("year", "run"), all.x = TRUE)
  sy <- merge(sy, big[, .(escapes = .N, area50 = sum(ha)), by = .(year, run)], by = c("year", "run"), all.x = TRUE)
  sy[is.na(ignitions), ignitions := 0L][is.na(escapes), `:=`(escapes = 0L, area50 = 0)]
  list(elf = elf, simYears = as.data.frame(sy), simSizes = as.data.frame(big), obsYears = as.data.frame(obsYears),
       obsSizes = as.data.frame(obsSizes), covariate = covariate, fit = fit)
}
