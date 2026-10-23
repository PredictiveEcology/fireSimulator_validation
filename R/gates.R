## Gates: pass/fail checks against fixed thresholds, no null. Each returns a one-row data.frame
## `(gate, pass, value, detail)`; `pass` is NA when the inputs needed are missing (not evaluated).

gateRow <- function(gate, pass, value, detail) data.frame(gate = gate, pass = pass, value = value, detail = detail)
inRatio <- function(r, lo = 0.25, hi = 4) is.finite(r) & r >= lo & r <= hi

#' G1: the scoring pipeline reproduces the fit's objective within replicate noise
#'
#' @param pipeline Objective values computed by the scoring pipeline, one per replicate (parameter set or rep).
#' @param ledger The fit's recorded objective (`objFunVal` of the ledger): one value, or one per kept set.
#' @return Passes when the ledger's mean lies within two standard errors of the pipeline mean. `value` is
#'   `mean(pipeline) - mean(ledger)`.
gateG1 <- function(pipeline, ledger) {
  pipeline <- pipeline[is.finite(pipeline)]; ledger <- ledger[is.finite(ledger)]
  if (length(pipeline) < 2L || !length(ledger)) return(gateRow("G1", NA, NA_real_, "pipeline objective not available"))
  d <- mean(pipeline) - mean(ledger); se <- stats::sd(pipeline) / sqrt(length(pipeline))
  gateRow("G1", abs(d) <= 2 * se, d, sprintf("difference %.3g, 2 x replicate SE %.3g", d, 2 * se))
}

#' G2: whole-ELF mean annual area burned within 0.25-4x of the observation
#'
#' @param simYears `data.frame(year, run, area50)`: simulated annual area (ha, fires >= 50 ha) per run.
#' @param obsYears `data.frame(year, area50)`: observed annual area.
#' @param years Years compared (default: the simulated years).
gateG2 <- function(simYears, obsYears, years = unique(simYears$year)) {
  sim <- mean(vapply(years, function(y) mean(simYears$area50[simYears$year == y]), 0))
  obs <- mean(obsYears$area50[obsYears$year %in% years])
  r <- if (obs > 0) sim / obs else if (sim == 0) 1 else Inf
  gateRow("G2", inRatio(r), r, sprintf("simulated %.0f ha/yr, observed %.0f ha/yr", sim, obs))
}

#' G3: the fit's fire set carries at least half of the observed area
#'
#' @param fitAreaHa Area (ha) of the fires the fit used. @param totalAreaHa Area (ha) of all observed fires (NBAC).
gateG3 <- function(fitAreaHa, totalAreaHa) {
  if (!is.finite(fitAreaHa) || !is.finite(totalAreaHa) || totalAreaHa <= 0) return(gateRow("G3", NA, NA_real_, "areas not available"))
  gateRow("G3", fitAreaHa / totalAreaHa >= 0.5, fitAreaHa / totalAreaHa,
          sprintf("fit set %.0f ha of %.0f ha observed", fitAreaHa, totalAreaHa))
}

#' G4: the largest simulated fire is within 10x of the largest observed fire
#'
#' @param simSizes `data.frame(ha)`: simulated fire sizes. @param obsSizes `data.frame(ha)`: observed sizes (whole record).
#' @return `value` = largest simulated / largest observed; the share of simulated area in fires over 100 kha is in
#'   `detail` and in the extra column `share100kha`.
gateG4 <- function(simSizes, obsSizes) {
  r <- max(simSizes$ha, 0) / max(obsSizes$ha, 0)
  share <- if (sum(simSizes$ha) > 0) sum(simSizes$ha[simSizes$ha > 1e5]) / sum(simSizes$ha) else 0
  out <- gateRow("G4", is.finite(r) && r <= 10, r, sprintf("largest simulated %.0f ha, largest observed %.0f ha, %.1f%% of simulated area in fires over 100 kha",
                  max(simSizes$ha, 0), max(obsSizes$ha, 0), 100 * share))
  out$share100kha <- share; out
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
#' Each link, simulated over observed on the scored years, must be within 0.25-4x; the failing link is named.
#'
#' @param simYears `data.frame(year, run, ignitions, escapes, area50)`. @param obsYears `data.frame(year, ignitions, escapes, area50)`.
#' @return `value` = the link ratio furthest from 1 (on the log scale); `detail` lists the three ratios and names the failing links.
gateG6 <- function(simYears, obsYears, years = unique(simYears$year)) {
  s <- simYears[simYears$year %in% years, ]; o <- obsYears[obsYears$year %in% years, ]
  nRun <- length(unique(s$run))
  link <- function(sNum, sDen, oNum, oDen) (sNum / sDen) / (oNum / oDen)
  ratios <- c(ignitions = (sum(s$ignitions) / nRun) / sum(o$ignitions),
              escapeRate = link(sum(s$escapes), sum(s$ignitions), sum(o$escapes), sum(o$ignitions)),
              areaPerEscape = link(sum(s$area50), sum(s$escapes), sum(o$area50), sum(o$escapes)))
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
