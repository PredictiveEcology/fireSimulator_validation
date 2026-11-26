## What the bias table and the climate-response table share: the pooling over ELFs, the rows of the long `results` table and the cells of the page.

medFinite <- function(v) stats::median(v[is.finite(v)])

## The range over ELFs: the 2.5% and 97.5% quantiles of the finite values (the same 95% as the range of single runs), NA with none
elfRange <- function(v) quantRange(v[is.finite(v)])
## The 2.5% and 97.5% quantiles (NAs left out; NA, NA with none): the one place of the 95% range
quantRange <- function(v) unname(stats::quantile(v, c(0.025, 0.975), na.rm = TRUE))

## The pooled row of ratios: `x` has one row per ELF with `ratio`, `obsRank` and `nYears`. `ratio`: the median over ELFs of the finite ratios, `elfLo`, `elfHi`: their 2.5-97.5% range (elfRange());
## `nELF`: the ELFs, `nNonFinite`: those left out as not finite; `nYears` summed, `obsRank` the median.
poolRatios <- function(x) {
  r <- elfRange(x$ratio)
  list(ratio = medFinite(x$ratio), elfLo = r[1], elfHi = r[2],
       nELF = nrow(x), nNonFinite = sum(!is.finite(x$ratio)), nYears = sum(x$nYears), obsRank = medFinite(x$obsRank))
}

## Rows of the long `results` table: `value`, `lo`, `hi` for the rows of `x` (`level`, `ELF`, `nELF`, `nYears`, `nRuns`); `nUnits` is the years of a pooled row, else the runs
longRows <- function(x, test, null, measure, value, lo, hi)
  data.frame(test = test, null = null, level = x$level, ELF = x$ELF, measure = measure, value = value, lo = lo, hi = hi,
             nELF = x$nELF, nUnits = ifelse(x$level == "pooled", x$nYears, x$nRuns), verdict = NA_character_)

## A cell of the page's tables: `c(value, low, high, obsRank)`, low and high the range of single runs (ELF rows) or over ELFs (pooled); `value` is the ratio (or, for a difference, `x$ratio`)
rangeCell <- function(x, pooled) ratioCell(x$ratio, if (pooled) x$elfLo else x$runLo, if (pooled) x$elfHi else x$runHi, x$obsRank)
