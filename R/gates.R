## Gates: pass/fail checks against fixed thresholds, no null. Each returns a one-row data.frame
## `(gate, pass, value, detail)`; `pass` is NA when the inputs needed are missing (not evaluated).

gateRow <- function(gate, pass, value, detail) data.frame(gate = gate, pass = pass, value = value, detail = detail)
inRatio <- function(r, lo = 0.25, hi = 4) is.finite(r) & r >= lo & r <= hi

#' G1: the scoring pipeline reproduces the fit's objective within replicate noise
#'
#' @param pipeline Objective values computed by the scoring pipeline, one per replicate (parameter set or rep).
#' @param ledger The fit's recorded objective (`objFunVal` of the ledger): one value, or one per kept set.
#' @return Passes when the ledger's median lies within the 2.5-97.5% range of the pipeline's replicates. `value` is
#'   `median(pipeline) - median(ledger)`. (Medians, not means, throughout the tests: one runaway replicate moves a mean.)
gateG1 <- function(pipeline, ledger) {
  pipeline <- pipeline[is.finite(pipeline)]; ledger <- ledger[is.finite(ledger)]
  if (length(pipeline) < 2L || !length(ledger)) return(gateRow("G1", NA, NA_real_, "pipeline objective not available"))
  d <- stats::median(pipeline) - stats::median(ledger); q <- stats::quantile(pipeline, c(0.025, 0.975), names = FALSE)
  gateRow("G1", stats::median(ledger) >= q[1] && stats::median(ledger) <= q[2], d,
          sprintf("difference %.3g; replicate 2.5-97.5%% range %.4g-%.4g", d, q[1], q[2]))
}

## one value per run: each run's simulated sizes or totals, so the median is taken over runs
perRun <- function(x, value, fun = sum) {
  if (!"run" %in% names(x)) x$run <- 1L
  if (!nrow(x)) return(numeric(0))
  as.numeric(tapply(x[[value]], x$run, fun))
}

#' G2: whole-ELF area burned over the scored years within 0.25-4x of the observation
#'
#' Each run's total area over the scored years; the median run is compared with the observed total over the same years.
#'
#' @param simYears `data.frame(year, run, area50)`: simulated annual area (ha, fires >= 50 ha) per run.
#' @param obsYears `data.frame(year, area50)`: observed annual area.
#' @param years Years compared (default: the simulated years).
gateG2 <- function(simYears, obsYears, years = unique(simYears$year)) {
  sim <- stats::median(perRun(simYears[simYears$year %in% years, , drop = FALSE], "area50"))
  obs <- sum(obsYears$area50[obsYears$year %in% years])
  r <- if (obs > 0) sim / obs else if (sim == 0) 1 else Inf
  gateRow("G2", inRatio(r), r, sprintf("median run %.0f ha, observed %.0f ha, over %d years", sim, obs, length(years)))
}

#' G3: the fit's fire set carries at least half of the observed area
#'
#' @param fitAreaHa Area (ha) of the fires the fit used. @param totalAreaHa Area (ha) of all observed fires (NBAC).
gateG3 <- function(fitAreaHa, totalAreaHa) {
  if (!is.finite(fitAreaHa) || !is.finite(totalAreaHa) || totalAreaHa <= 0) return(gateRow("G3", NA, NA_real_, "areas not available"))
  gateRow("G3", fitAreaHa / totalAreaHa >= 0.5, fitAreaHa / totalAreaHa,
          sprintf("fit set %.0f ha of %.0f ha observed", fitAreaHa, totalAreaHa))
}

#' G4: the large simulated fires are within 10x of the large observed fires
#'
#' Each run's median of its 10 largest fires over the scored years (`topFires()`); the median run is compared with the
#' same statistic of the observed fires in those years (all observed years when `obsSizes` has no `year`). One largest
#' fire is too noisy a target (Eliot, 2026-10-26); it stays in `detail` for information.
#'
#' @param simSizes `data.frame(year, run, ha)`: simulated fire sizes. @param obsSizes `data.frame(year, ha)`: observed sizes.
#' @param k Number of largest fires.
#' @return `value` = median run's median of its `k` largest / observed median of its `k` largest; the share of simulated
#'   area in fires over 100 kha is in `detail` and in the extra column `share100kha`.
gateG4 <- function(simSizes, obsSizes, k = 10L) {
  if ("year" %in% names(simSizes) && "year" %in% names(obsSizes)) obsSizes <- obsSizes[obsSizes$year %in% unique(simSizes$year), , drop = FALSE]
  simTop <- if (nrow(simSizes)) stats::median(perRun(simSizes, "ha", function(x) topFires(x, k))) else 0
  obsTop <- topFires(obsSizes$ha, k)
  simMax <- if (nrow(simSizes)) stats::median(perRun(simSizes, "ha", max)) else 0
  r <- simTop / obsTop
  share <- if (sum(simSizes$ha) > 0) sum(simSizes$ha[simSizes$ha > 1e5]) / sum(simSizes$ha) else 0
  out <- gateRow("G4", is.finite(r) && r <= 10, r, sprintf(paste("median of the %d largest fires: median run %.0f ha, observed %.0f ha;",
                  "largest fire of the median run %.0f ha, largest observed %.0f ha; %.1f%% of simulated area in fires over 100 kha"),
                  k, simTop, obsTop, simMax, max(obsSizes$ha, 0), 100 * share))
  out$share100kha <- share; out
}

#' Median of the `k` largest values (all of them when there are fewer than `k`; 0 when none)
#'
#' @param x Fire sizes. @param k Number of largest.
#' @return One number.
topFires <- function(x, k = 10L) {
  if (!length(x)) return(0)
  stats::median(utils::head(sort(x, decreasing = TRUE), k))
}

#' G5: the power-law threshold p_c is below p50 in every kept parameter set
#'
#' @param pc `p_c` of each kept set (ledger `maxAsymptote`). @param p50 The ELF's calibrated threshold (ledger `powerLawPc50_spread`).
gateG5 <- function(pc, p50) {
  if (!length(pc) || !length(p50) || !is.finite(p50[1])) return(gateRow("G5", NA, NA_real_, "no power-law columns in the ledger"))
  gateRow("G5", all(pc < p50[1]), max(pc) / p50[1], sprintf("max p_c %.4f, p50 %.4f", max(pc), p50[1]))
}

#' G6: chain ignitions x escape rate x area per escape
#'
#' Each link, the median run's value over observed on the scored years, must be within 0.25-4x; the failing link is named.
#'
#' @param simYears `data.frame(year, run, ignitions, escapes, area50)`. @param obsYears `data.frame(year, ignitions, escapes, area50)`.
#' @return `value` = the link ratio furthest from 1 (on the log scale); `detail` lists the three ratios and names the failing links.
gateG6 <- function(simYears, obsYears, years = unique(simYears$year)) {
  s <- simYears[simYears$year %in% years, ]; o <- obsYears[obsYears$year %in% years, ]
  ## each run's totals over the scored years, then the median run's link
  ign <- perRun(s, "ignitions"); esc <- perRun(s, "escapes"); area <- perRun(s, "area50")
  med <- function(x) { x <- x[is.finite(x)]; if (length(x)) stats::median(x) else NA_real_ }
  ratios <- c(ignitions = med(ign) / sum(o$ignitions),
              escapeRate = med(esc / ign) / (sum(o$escapes) / sum(o$ignitions)),
              areaPerEscape = med(area / esc) / (sum(o$area50) / sum(o$escapes)))
  ratios[!is.finite(ratios)] <- NA_real_
  ok <- inRatio(ratios)
  failing <- names(ratios)[!ok]
  pass <- if (all(is.na(ratios))) NA else length(failing) == 0L
  gateRow("G6", pass, if (all(is.na(ratios))) NA_real_ else ratios[which.max(abs(log(ratios)))],
          paste0(paste0(names(ratios), " ", signif(ratios, 3), collapse = "; "),
                 if (length(failing)) paste0("; failing: ", paste(failing, collapse = ", ")) else ""))
}

#' Run the gates of one ELF
#'
#' @param elf An ELF input list (see [runWave1()]).
#' @return `data.frame(gate, pass, value, detail)` for G1-G6, with an attribute `meaningful`: `FALSE` when G2 or G3 fails
#'   (the skill tests of this ELF are not meaningful).
runGates <- function(elf) {
  f <- elf$fit
  g4 <- gateG4(elf$simSizes, elf$obsSizes)
  g <- rbind(gateG1(f$pipelineObjective, f$objective),
             gateG2(elf$simYears, elf$obsYears),
             gateG3(f$fitAreaHa %||% NA_real_, f$totalAreaHa %||% NA_real_),
             g4[, c("gate", "pass", "value", "detail")],
             gateG5(f$pc, f$p50),
             gateG6(elf$simYears, elf$obsYears))
  g$share100kha <- ifelse(g$gate == "G4", g4$share100kha, NA_real_)
  attr(g, "meaningful") <- !any(g$pass[g$gate %in% c("G2", "G3")] %in% FALSE)
  g
}

`%||%` <- function(a, b) if (is.null(a)) b else a
