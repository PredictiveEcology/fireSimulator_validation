## Bias: the median run over the observation, as numbers, for annual area burned, the 10 largest fires, ignitions and escapes.
## The gates only pass or fail these ratios (0.25-4x, 10x); a bias compounds over a long run, so it is reported as a number.

biasQuantities <- c(area = "area burned", top10 = "10 largest fires", ignitions = "ignitions", escapes = "escapes")
## the unit of the totals' `per1000yr`
biasUnits <- c(area = "ha per 1000 yr", top10 = NA_character_, ignitions = "fires per 1000 yr", escapes = "escapes per 1000 yr")

#' Bias of one ELF: the median run over the observation
#'
#' For each of area burned (`area50`), the median of the 10 largest fires, ignitions and escapes: the ratio of the median over runs of the run's value
#' to the observed value, over the scored years. A run's value is its total over the scored years (area, ignitions, escapes) or the
#' median of its 10 largest fires (as G4, [gateG4()], which this calls). Totals use the statistic of G2 ([totalRatio()]) and the
#' resampling of [bootRuns()].
#'
#' @param elf An ELF input list ([runWave1()]). @param nBoot,seed Bootstrap over runs.
#' @return `data.frame`, one row per quantity: `level` ("ELF"), `ELF`, `quantity`, `observed`, `simMedian` (the median run's value; the same units as `observed`),
#'   `ratio`, `log2ratio`, `lo`, `hi` (bootstrap interval of the median-run ratio), `runLo`, `runHi` (2.5% and 97.5% of the single runs' value over
#'   `observed`), `elfLo`, `elfHi`, `ratioOfSums` (NA here; see [poolBias()]), `obsRank` (share of runs below the observation, 0-1; ties count half),
#'   `nRuns`, `nYears`, `per1000yr` (the implied difference, median run minus observed, accumulated over 1000 years: `(ratio - 1) * observed per year * 1000`;
#'   NA for the largest fires) and `unit` of `per1000yr`.
scoreBias <- function(elf, nBoot = 1000, seed = 1) {
  years <- unique(elf$simYears$year)
  o <- elf$obsYears[elf$obsYears$year %in% years, , drop = FALSE]
  g4 <- gateG4(elf$simSizes, elf$obsSizes, nBoot = nBoot, seed = seed)
  row <- function(q, observed, simMedian, b, i = 1L)
    data.frame(level = "ELF", ELF = elf$elf, quantity = q, observed = observed, simMedian = simMedian, ratio = b$est[i], log2ratio = log2(b$est[i]),
               lo = b$lo[i], hi = b$hi[i], runLo = b$runLo[i], runHi = b$runHi[i], elfLo = NA_real_, elfHi = NA_real_, ratioOfSums = NA_real_,
               obsRank = b$obsRank[i], nRuns = b$nRuns, nYears = length(years),
               per1000yr = if (q == "top10") NA_real_ else (simMedian - observed) / length(years) * 1000, unit = unname(biasUnits[q]))
  tot <- function(q, col) {
    b <- bootRuns(elf$simYears[elf$simYears$year %in% years, , drop = FALSE], totalRatio(col, sum(o[[col]])), nBoot, seed)
    row(q, sum(o[[col]]), b$est[2], b)
  }
  obsTop <- topFires(elf$obsSizes$ha[elf$obsSizes$year %in% years])
  top <- row("top10", obsTop, g4$value * obsTop, list(est = g4$value, lo = g4$lo, hi = g4$hi, runLo = g4$runLo, runHi = g4$runHi,
                                                      obsRank = g4$obsRank, nRuns = g4$nRuns))
  rbind(tot("area", "area50"), top, tot("ignitions", "ignitions"), tot("escapes", "escapes"))
}

#' Bias pooled over ELFs
#'
#' Two pooled ratios, not the same number. `ratio` is the median over ELFs of the ELF ratios (`elfLo`, `elfHi`: their range, lowest and highest ELF);
#' it treats every ELF alike. `ratioOfSums` is the sum over ELFs of the median run's total over the sum of the observed totals (NA for the largest
#' fires, which have no total); large ELFs dominate it. `per1000yr` is the sum over ELFs of the ELFs' differences.
#'
#' @param bias The rows of [scoreBias()] of several ELFs, bound.
#' @return `data.frame` with the columns of [scoreBias()], one row per quantity, `level` "pooled" and `ELF` "all". `observed` and `simMedian` are sums over ELFs (NA
#'   for the largest fires), `nRuns` is NA, `nYears` the sum, `obsRank` the median of the ELF values.
poolBias <- function(bias) {
  do.call(rbind, lapply(unique(bias$quantity), function(q) {
    x <- bias[bias$quantity == q, , drop = FALSE]; tot <- q != "top10"; med <- function(v) stats::median(v[is.finite(v)])
    r <- x$ratio[is.finite(x$ratio)]
    data.frame(level = "pooled", ELF = "all", quantity = q, observed = if (tot) sum(x$observed) else NA_real_, simMedian = if (tot) sum(x$simMedian) else NA_real_,
               ratio = med(x$ratio), log2ratio = med(x$log2ratio), lo = NA_real_, hi = NA_real_, runLo = NA_real_, runHi = NA_real_,
               elfLo = min(r), elfHi = max(r), ratioOfSums = if (tot) sum(x$simMedian) / sum(x$observed) else NA_real_,
               obsRank = med(x$obsRank), nRuns = NA_integer_, nYears = sum(x$nYears), per1000yr = if (tot) sum(x$per1000yr) else NA_real_,
               unit = unname(biasUnits[q]))
  }))
}

#' The bias of every ELF, and pooled
#'
#' @param elfs List of ELF inputs ([runWave1()]). @inheritParams scoreBias
#' @return The rows of [scoreBias()] of every ELF, then those of [poolBias()].
biasTable <- function(elfs, nBoot = 1000, seed = 1) {
  b <- do.call(rbind, lapply(elfs, scoreBias, nBoot = nBoot, seed = seed))
  out <- rbind(b, poolBias(b)); rownames(out) <- NULL
  out
}

## `bias` as rows of the long `results` table of runWave1(): test "B", the quantity in `null`, measure "ratio" and "log2ratio"; `lo` and `hi`
## are the single-run range (ELF rows) or the range over ELFs (pooled rows)
biasResultRows <- function(bias) {
  pooled <- bias$level == "pooled"
  do.call(rbind, lapply(c("ratio", "log2ratio"), function(m) {
    lo <- ifelse(pooled, bias$elfLo, bias$runLo); hi <- ifelse(pooled, bias$elfHi, bias$runHi)
    if (m == "log2ratio") { lo <- log2(lo); hi <- log2(hi) }
    data.frame(test = "B", null = bias$quantity, level = bias$level, ELF = bias$ELF, measure = m, value = bias[[m]], lo = lo, hi = hi,
               nELF = ifelse(pooled, length(unique(bias$ELF[!pooled])), 1L), nUnits = ifelse(pooled, bias$nYears, bias$nRuns), verdict = NA_character_)
  }))
}

## the rows of the page's bias table: one per ELF, then pooled twice (median of the ELF ratios with their range, and the ratio of sums). A ratio
## is `c(ratio, low, high, obsRank)`, low and high the range of single runs (ELF rows) or over ELFs (pooled); `<q>Diff` is the total's per1000yr
biasPageRows <- function(bias) {
  if (is.null(bias)) return(list())
  cell <- function(x, pooled) round(c(x$ratio, if (pooled) x$elfLo else x$runLo, if (pooled) x$elfHi else x$runHi, x$obsRank), 3)
  rowOf <- function(label, x, pooled, sums = FALSE) {
    r <- list(elf = label)
    for (q in names(biasQuantities)) {
      z <- x[x$quantity == q, ]
      r[[q]] <- if (sums) round(c(z$ratioOfSums, NA, NA, NA), 3) else cell(z, pooled)
      if (q != "top10") r[[paste0(q, "Diff")]] <- signif(z$per1000yr, 3)
    }
    r
  }
  c(lapply(unique(bias$ELF[bias$level == "ELF"]), function(id) rowOf(id, bias[bias$level == "ELF" & bias$ELF == id, ], FALSE)),
    list(rowOf("all: median of ELFs", bias[bias$level == "pooled", ], TRUE), rowOf("all: ratio of sums", bias[bias$level == "pooled", ], TRUE, sums = TRUE)))
}
