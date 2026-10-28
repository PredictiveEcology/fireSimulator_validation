## Hindcast validation figures for one simulated series against the observed record. Plain tables in (the ELF input list of
## runWave1(); the wide summary of summariseGates()), ggplot objects out. A run is one `run` id of `simYears` over all its
## years. Observed values are taken in the simulated years only, as the gates do. Statistical choices: medians; densities on
## log10 with one bandwidth (the observed one's x 1.2) so 50 observed fires are not smoothed more than thousands of simulated
## ones, reflected at a hard lower bound; spread of the simulated series = single runs, because the observed record is one
## realisation.

## sizes (inches, dpi) of the figures for Plots(ggsaveArgs = )
panelSize <- list(width = 16, height = 12, dpi = 110, units = "in", bg = "white")
overviewSize <- list(width = 16, height = 9, dpi = 110, units = "in", bg = "white")

## the five ratios of the gates, in the order the figures show them: key (summary column prefix), label
gateRatios <- c(area = "area", top10_median = "median of 10 largest fires", ignitions = "ignitions",
                escapeRate = "escape rate", areaPerEscape = "area per escape")

## colours and line types of the observed and the simulated series; `label` names the simulated one
seriesStyle <- function(label) {
  lv <- c("observed", label)
  list(levels = lv, cols = stats::setNames(c("black", "#0072B2"), lv), lts = stats::setNames(c("solid", "dashed"), lv))
}

## The tables every figure draws on. `grid`: one row per simulated run-year. Areas and sizes are of fires of `minHa` or more.
hindcastSeries <- function(elf, minHa = 50) {
  sy <- elf$simYears; yrs <- sort(unique(sy$year))
  grid <- unique(sy[, c("year", "run")]); grid <- grid[order(grid$run, grid$year), ]
  key <- function(d) paste(d$year, d$run)
  bigSim <- elf$simSizes[elf$simSizes$ha >= minHa, , drop = FALSE]
  a <- tapply(bigSim$ha, key(bigSim), sum)
  grid$area <- as.numeric(a[key(grid)]); grid$area[is.na(grid$area)] <- 0
  oy <- elf$obsYears[elf$obsYears$year %in% yrs, , drop = FALSE]; oy <- oy[order(oy$year), ]
  bigObs <- elf$obsSizes[elf$obsSizes$year %in% yrs & elf$obsSizes$ha >= minHa, , drop = FALSE]
  oa <- tapply(bigObs$ha, factor(bigObs$year, oy$year), sum); oy$area <- as.numeric(oa); oy$area[is.na(oy$area)] <- 0
  merged <- merge(grid, sy[, c("year", "run", "ignitions", "escapes")], by = c("year", "run"))
  list(years = yrs, grid = merged, obsYears = oy, simSizes = bigSim, obsSizes = bigObs)
}

## values of one run (all its years pooled), one list element per run; `v` is aligned with `d` (run-year rows)
byRun <- function(v, run) split(v, factor(run, unique(run)))

## Single-run spread: `vals` = list, one element per run; `stat` maps ONE run's values to a numeric vector. Returns the
## elements x runs matrix.
perRunStat <- function(vals, stat) do.call(cbind, lapply(vals, stat))

## pointwise 2.5-97.5% across runs (NA where fewer than 3 runs have a value)
envelope <- function(m) {
  q <- apply(m, 1, function(z) if (sum(!is.na(z)) < 3) c(NA_real_, NA_real_) else stats::quantile(z, c(0.025, 0.975), names = FALSE, na.rm = TRUE))
  list(lo = q[1, ], hi = q[2, ])
}

## density of log10(v > 0) on [from, to] (256 points), reflected at `lo` when given; area = share of the `ntot` values
kde <- function(v, bw, lo, from, to, ntot) {
  lx <- log10(v[v > 0])
  k <- if (is.null(lo)) stats::density(lx, bw = bw, from = from, to = to, n = 256)
       else stats::density(c(lx, 2 * lo - lx), bw = bw, from = from, to = to, n = 256)
  k$y * (if (is.null(lo)) 1 else 2) * length(lx) / ntot
}

emptyDens <- data.frame(series = character(), ha = numeric(), y = numeric(), lo = numeric(), hi = numeric())

## Density figure of a count or an area: `obs` observed values (zeros included), `runVals` simulated values per run.
dens <- function(obs, runVals, st, xlab, lower = NULL, title, subtitle = NULL, caption = NULL, baseSize = 15) {
  sim <- unlist(runVals, use.names = FALSE); lo <- if (!is.null(lower)) log10(lower)
  bwOf <- function(v) { l <- log10(v[v > 0]); if (length(unique(l)) > 1) stats::bw.nrd0(l) * 1.2 }
  bw <- bwOf(obs); if (is.null(bw)) bw <- bwOf(sim)
  rl <- list(); d <- emptyDens
  one <- function(v, s) {
    lx <- log10(v[v > 0])
    if (is.null(bw) || length(lx) < 2) return(NULL)
    from <- if (is.null(lo)) min(lx) else lo; to <- max(lx)
    if (!to > from) return(NULL)
    x <- seq(from, to, length.out = 256)
    r <- data.frame(series = s, ha = 10^x, y = kde(v, bw, lo, from, to, length(v)), lo = NA_real_, hi = NA_real_)
    if (s != "observed") {
      m <- perRunStat(runVals, function(z) { p <- log10(z[z > 0]); if (length(p) < 2) return(rep(NA_real_, length(x)))
        y <- kde(z, bw, lo, from, to, length(z)); y[x > max(p) | (is.null(lo) & x < min(p))] <- NA_real_; y })
      ee <- envelope(m); r$lo <- ee$lo; r$hi <- ee$hi
      rl[[s]] <<- data.frame(series = s, run = rep(seq_len(ncol(m)), each = nrow(m)), ha = 10^x, y = as.vector(m))
    }
    r
  }
  d <- do.call(rbind, c(list(emptyDens), list(one(obs, "observed"), one(sim, st$levels[2]))))
  d$series <- factor(d$series, st$levels)
  ggplot(d, aes(ha, y, colour = series, linetype = series)) +
    geom_ribbon(data = d[!is.na(d$lo), ], mapping = aes(x = ha, ymin = lo, ymax = hi, fill = series), colour = NA, alpha = 0.2, inherit.aes = FALSE) +
    (if (length(rl)) geom_line(data = do.call(rbind, rl), aes(ha, y, group = run), colour = "grey45", alpha = 0.25, linewidth = 0.3, inherit.aes = FALSE, na.rm = TRUE)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_rug(data = data.frame(ha = obs[obs > 0]), aes(x = ha), inherit.aes = FALSE, sides = "b", colour = "black", alpha = 0.6) +
    scale_x_log10(labels = scales::label_comma()) +
    scale_colour_manual(values = st$cols, drop = FALSE, name = NULL) + scale_fill_manual(values = st$cols, guide = "none", drop = FALSE) +
    scale_linetype_manual(values = st$lts, drop = FALSE, name = NULL) +
    labs(title = title, subtitle = subtitle, x = xlab, y = "density", caption = caption) +
    theme_bw(base_size = baseSize) + theme(legend.position = "bottom")
}

## Bars of the share of values per class, with 95% single-run error bars on the simulated series (each run's own shares).
barsBoot <- function(obs, runVals, st, cutf, title, xlab, ylab, subtitle = NULL, caption = NULL, baseSize = 15) {
  lv <- levels(cutf(1e9)); shr <- function(x) if (length(x)) as.numeric(table(cutf(x)) / length(x)) else rep(NA_real_, length(lv))
  ee <- envelope(perRunStat(runVals, shr))
  b <- rbind(data.frame(series = "observed", cls = factor(lv, lv), share = shr(obs), lo = NA_real_, hi = NA_real_),
             data.frame(series = st$levels[2], cls = factor(lv, lv), share = shr(unlist(runVals, use.names = FALSE)), lo = ee$lo, hi = ee$hi))
  b$series <- factor(b$series, st$levels)
  ggplot(b, aes(cls, share, fill = series)) + geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_errorbar(aes(ymin = lo, ymax = hi, group = series), position = position_dodge(0.8), width = 0.25, colour = "grey25", na.rm = TRUE) +
    scale_fill_manual(values = st$cols, drop = FALSE, name = NULL) + scale_y_continuous(labels = scales::label_percent()) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab, caption = caption) + theme_bw(base_size = baseSize) + theme(legend.position = "bottom")
}

## escape classes
cutEs <- function(x) cut(x, c(-Inf, 0, 1, 2, 4, 9, Inf), labels = c("0", "1", "2", "3-4", "5-9", "10+"))

elfSubtitle <- function(elf, ...) paste(c(paste("ELF", elf$elf), ...), collapse = "; ")
densNote <- "area under a curve = share of non-zero values; one bandwidth for both series; ticks: observed"
barCaption <- "error bars: 95% range across single runs"
runCaption <- "thin lines: single runs; shaded: 95% range across single runs"

#' Area burned by year: observed against the simulated median and 10-90% of runs
#'
#' Fires of `minHa` or more, summed per year (a year without such fires is 0), log y of area + 1. The observed record is
#' taken in the simulated years. The five `plot*` functions each return a ggplot of one ELF.
#'
#' @param elf One ELF input (see [runWave1()]): `elf`, `simYears`, `simSizes`, `obsYears`, `obsSizes`.
#' @param label Name of the simulated series in the legend.
#' @param minHa Smallest fire counted (ha).
#' @param baseSize Base font size of the theme.
#' @return A ggplot.
plotAreaByYear <- function(elf, label = "simulated", minHa = 50, baseSize = 15) {
  s <- hindcastSeries(elf, minHa); st <- seriesStyle(label)
  q <- do.call(rbind, lapply(split(s$grid$area, s$grid$year), function(a)
    data.frame(med = stats::median(a), lo = stats::quantile(a, 0.1, names = FALSE), hi = stats::quantile(a, 0.9, names = FALSE))))
  q$year <- as.integer(rownames(q)); q$series <- factor(label, st$levels)
  o <- data.frame(year = s$obsYears$year, area = s$obsYears$area, series = factor("observed", st$levels))
  ggplot(q, aes(year, med + 1, colour = series)) + geom_pointrange(aes(ymin = lo + 1, ymax = hi + 1), size = 0.5, linewidth = 0.7) +
    geom_point(data = o, aes(year, area + 1), shape = 4, size = 3.5, stroke = 1.3) +
    scale_y_log10(labels = scales::label_comma()) + scale_colour_manual(values = st$cols, drop = FALSE, name = NULL) +
    labs(title = "Area burned by year", subtitle = elfSubtitle(elf, paste0("fires of ", minHa, " ha or more")), x = NULL, y = "area + 1 (ha, log)",
         caption = "x: observed; dots and bars: simulated median and 10-90% of runs") +
    theme_bw(base_size = baseSize) + theme(legend.position = "bottom")
}

#' Density of annual area burned
#'
#' Annual area in fires of `minHa` or more, log10 x, reflected at `minHa`, one bandwidth for both series (the observed one's x
#' 1.2). Thin lines are single runs, the ribbon their pointwise 95% range; the shares of years without such fires are in the subtitle.
#'
#' @inheritParams plotAreaByYear
#' @param showElf Put the ELF id in the subtitle (the panel has it in its title).
plotAnnualAreaDensity <- function(elf, label = "simulated", minHa = 50, baseSize = 15, showElf = TRUE) {
  s <- hindcastSeries(elf, minHa); st <- seriesStyle(label)
  zero <- sprintf("years with none: observed %.0f%%, %s %.0f%%", 100 * mean(s$obsYears$area == 0), label, 100 * mean(s$grid$area == 0))
  dens(s$obsYears$area, byRun(s$grid$area, s$grid$run), st, "area burned per year (ha, log scale)", lower = minHa,
       title = "Annual area burned", subtitle = if (showElf) elfSubtitle(elf, zero) else zero,
       caption = paste0("fires of ", minHa, " ha or more; ", densNote, "; reflected at ", minHa, " ha; ", runCaption), baseSize = baseSize)
}

#' Density of fire sizes
#'
#' Sizes of fires of `minHa` or more, log10 x, reflected at `minHa`, one bandwidth for both series; thin lines are single runs, the
#' ribbon their pointwise 95% range, ticks the observed fires. A run without such fires has no line.
#'
#' @inheritParams plotAreaByYear
plotFireSizes <- function(elf, label = "simulated", minHa = 50, baseSize = 15) {
  s <- hindcastSeries(elf, minHa); st <- seriesStyle(label)
  runs <- split(s$simSizes$ha, factor(s$simSizes$run, unique(s$grid$run)))
  dens(s$obsSizes$ha, runs, st, "fire size (ha, log scale)", lower = minHa, title = "Fire sizes",
       subtitle = elfSubtitle(elf, paste0("fires of ", minHa, " ha or more; ", nrow(s$obsSizes), " observed")),
       caption = paste0("fires of ", minHa, " ha or more; ", densNote, "; reflected at ", minHa, " ha; ", runCaption), baseSize = baseSize)
}

#' Density of ignitions per year
#'
#' Log10 x with one bandwidth for both series; thin lines are single runs, the ribbon their pointwise 95% range.
#'
#' @inheritParams plotAreaByYear
plotIgnitions <- function(elf, label = "simulated", baseSize = 15) {
  s <- hindcastSeries(elf, 0); st <- seriesStyle(label)
  dens(s$obsYears$ignitions, byRun(s$grid$ignitions, s$grid$run), st, "ignitions per year (log scale)",
       title = "Ignitions per year", subtitle = elfSubtitle(elf), caption = paste0(densNote, "; ", runCaption), baseSize = baseSize)
}

#' Share of years by number of escapes
#'
#' An escape is a fire of `minHa` or more (`escapes` of `simYears` and `obsYears`). A small discrete count, so bars of the share
#' of years per class; error bars: 95% range across single runs.
#'
#' @inheritParams plotAreaByYear
plotEscapes <- function(elf, label = "simulated", baseSize = 15) {
  s <- hindcastSeries(elf, 0); st <- seriesStyle(label)
  barsBoot(s$obsYears$escapes, byRun(s$grid$escapes, s$grid$run), st, cutEs, "Escapes per year", "escapes in the year", "share of years",
           subtitle = elfSubtitle(elf), caption = barCaption, baseSize = baseSize)
}

fmtRatio <- function(x) sprintf("%.2f", x)

## "area 0.72 (0.31-1.90), ..." for one row of the summary; the whisker is the single-run 95% range
ratioText <- function(row) {
  txt <- vapply(names(gateRatios), function(k)
    sprintf("%s %s (%s-%s)", gateRatios[[k]], fmtRatio(row[[paste0(k, "_ratio")]]), fmtRatio(row[[paste0(k, "_runLo")]]), fmtRatio(row[[paste0(k, "_runHi")]])), "")
  paste0("simulated / observed, median run (single-run 95% range): ", paste(txt[1:3], collapse = ", "), ",\n", paste(txt[4:5], collapse = ", "))
}

metrics <- c(areaByYear = "Area burned by year", annualArea = "Annual area burned", fireSizes = "Fire sizes",
             ignitions = "Ignitions per year", escapes = "Escapes per year")
## functions, not objects: ggplot2 is attached after R/ is sourced
noLeg <- function() guides(colour = "none", fill = "none", linetype = "none")
panelTheme <- function() theme(plot.title = element_text(size = 22, face = "bold"), plot.subtitle = element_text(size = 14), plot.caption = element_text(size = 12))
shareCaption <- function(minHa) paste0("fires of ", minHa, " ha or more; densities: ", densNote, "; reflected at ", minHa, " ha;\n", runCaption,
                                      "; ", barCaption, "; area by year: x = observed, dots and bars = simulated median and 10-90% of runs")

#' The five plots of one ELF
#'
#' Built once, then laid out by [plotELFPanel()] and [plotMetricPanel()].
#'
#' @inheritParams plotAreaByYear
#' @return Named list of ggplots: `areaByYear`, `annualArea`, `fireSizes`, `ignitions`, `escapes`.
elfPlots <- function(elf, label = "simulated", minHa = 50, baseSize = 15)
  list(areaByYear = plotAreaByYear(elf, label, minHa, baseSize), annualArea = plotAnnualAreaDensity(elf, label, minHa, baseSize, showElf = FALSE),
       fireSizes = plotFireSizes(elf, label, minHa, baseSize), ignitions = plotIgnitions(elf, label, baseSize), escapes = plotEscapes(elf, label, baseSize))

#' One figure per ELF: the five hindcast plots
#'
#' [plotAreaByYear()] full width on top, then [plotAnnualAreaDensity()], [plotFireSizes()], [plotIgnitions()] and [plotEscapes()]
#' in two columns, one shared legend at the bottom and one title, "ELF <id>", with the ratios of the ELF in the subtitle.
#' Made for 16 x 12 inches at 110 dpi (`panelSize`).
#'
#' @inheritParams plotAreaByYear
#' @param summary The wide gate summary ([summariseGates()]); the row of `elf$elf` gives the ratios in the subtitle. Without it there is no subtitle.
#' @param plots [elfPlots()] of `elf`, to lay out plots already built.
#' @return A patchwork object.
plotELFPanel <- function(elf, summary = NULL, label = "simulated", minHa = 50, baseSize = 15, plots = elfPlots(elf, label, minHa, baseSize)) {
  bare <- labs(subtitle = NULL, caption = NULL)
  p <- list(plots$areaByYear + bare + noLeg(), plots$annualArea + labs(caption = NULL), plots$fireSizes + bare + noLeg(),
            plots$ignitions + bare + noLeg(), plots$escapes + bare + noLeg())
  row <- if (!is.null(summary)) summary[as.character(summary$ELF) == as.character(elf$elf), , drop = FALSE]
  patchwork::wrap_plots(p, design = "AA\nBC\nDE", heights = c(1, 1.15, 1.15)) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = paste("ELF", elf$elf), subtitle = if (!is.null(row) && nrow(row)) ratioText(row[1, ]), caption = shareCaption(minHa),
                               theme = panelTheme()) &
    theme(legend.position = "bottom", legend.text = element_text(size = 14))
}

## ELF ids ordered by area ratio, lowest first (the order of plotOverview(), bottom to top)
elfOrder <- function(summary) as.character(summary$ELF)[order(summary$area_ratio, na.last = TRUE)]

metricCols <- 3L
## device size for [plotMetricPanel()]: about 4 inches per row of panes
metricSize <- function(n, ncol = metricCols) c(panelSize[names(panelSize) != "height"], list(height = 4 * ceiling(n / ncol) + 2))

#' One metric, every ELF: a pane per ELF
#'
#' The same plot of each ELF ([elfPlots()]), laid out in `ncol` columns, ordered by area ratio as in [plotOverview()]; pane title =
#' ELF id, one shared legend and caption. Size: [metricSize()].
#'
#' @param plots Named list (names = ELF ids) of [elfPlots()].
#' @param summary The wide gate summary ([summariseGates()]), for the order of the panes.
#' @param metric One of `names(metrics)`: `areaByYear`, `annualArea`, `fireSizes`, `ignitions`, `escapes`.
#' @param ncol Columns of panes.
#' @param minHa Only for the caption.
#' @return A patchwork object.
plotMetricPanel <- function(plots, summary, metric = names(metrics), ncol = metricCols, minHa = 50) {
  metric <- match.arg(metric)
  ids <- intersect(elfOrder(summary), names(plots))
  p <- lapply(seq_along(ids), function(i) plots[[ids[i]]][[metric]] + labs(title = ids[i], subtitle = NULL, caption = NULL) +
                                          (if (i > 1) noLeg()))
  patchwork::wrap_plots(p, ncol = ncol) + patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = paste0(metrics[[metric]], ": every ELF, ordered by area ratio"), caption = shareCaption(minHa), theme = panelTheme()) &
    theme(legend.position = "bottom", legend.text = element_text(size = 14))
}

#' All ELFs in one figure: the five gate ratios
#'
#' One facet per ratio (simulated / observed). Per ELF a point at the median run's ratio and a horizontal whisker over the 95%
#' range of single runs (`runLo`, `runHi`), log x, a reference line at 1, ELFs ordered by area ratio. Ratios outside 0.01-100
#' are drawn at the limit.
#'
#' @param summary Wide table, one row per ELF, as [summariseGates()] returns: `ELF` and, for each of `area`, `top10_median`,
#'   `ignitions`, `escapeRate`, `areaPerEscape`, the columns `<key>_ratio`, `<key>_runLo`, `<key>_runHi`.
#' @param baseSize Base font size of the theme.
#' @return A ggplot.
plotOverview <- function(summary, baseSize = 15) {
  lim <- c(0.01, 100); clamp <- function(x) pmin(pmax(x, lim[1]), lim[2])
  ord <- elfOrder(summary)
  d <- do.call(rbind, lapply(names(gateRatios), function(k) data.frame(
    ELF = factor(as.character(summary$ELF), ord), measure = factor(gateRatios[[k]], gateRatios),
    ratio = clamp(summary[[paste0(k, "_ratio")]]), lo = clamp(summary[[paste0(k, "_runLo")]]), hi = clamp(summary[[paste0(k, "_runHi")]]))))
  ggplot(d, aes(ratio, ELF)) + geom_vline(xintercept = 1, colour = "grey40") +
    geom_pointrange(aes(xmin = lo, xmax = hi), orientation = "y", colour = "#0072B2", size = 0.5, linewidth = 0.7, na.rm = TRUE) +
    facet_wrap(~measure, nrow = 1) +
    scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = scales::label_comma()) +
    labs(title = "Hindcast against observed, all ELFs", subtitle = "simulated / observed: median run (dot) and 95% range of single runs (whisker); 1 = agreement",
         x = "ratio (log scale)", y = "ELF (ordered by area ratio)", caption = "ratios outside 0.01-100 are drawn at the limit") +
    theme_bw(base_size = baseSize) + theme(plot.title = element_text(size = 22, face = "bold"), strip.text = element_text(size = 13), panel.grid.minor = element_blank())
}

#' The gates of [runWave1()] as one row per ELF with the five ratios and their intervals
#'
#' Area (G2) and the median of the 10 largest fires (G4) are the gate `value`; ignitions, escape rate and area per escape are
#' the three links of G6, read from its `detail`. For each ratio `<key>_ratio`, `_lo`, `_hi` (median-run 95%) and `_runLo`, `_runHi`
#' (single runs 95%).
#'
#' @param gates The `gates` table of [runWave1()].
#' @return A data.frame with `ELF` and the columns above.
summariseGates <- function(gates) {
  num <- "(?:NA|NaN|Inf|[0-9.]+(?:e[-+][0-9]+)?)"
  link <- function(detail, k) {
    m <- regmatches(detail, regexec(paste0(k, " (", num, ") \\(median-run 95%: (", num, ")-(", num, "); single runs 95%: (", num, ")-(", num, ")\\)"), detail, perl = TRUE))[[1]]
    stats::setNames(as.numeric(m[2:6]), paste0(k, c("_ratio", "_lo", "_hi", "_runLo", "_runHi")))
  }
  do.call(rbind, lapply(split(gates, factor(gates$ELF, unique(gates$ELF))), function(g) {
    gate <- function(id, k) { r <- g[g$gate == id, ]
      stats::setNames(c(r$value, r$lo, r$hi, r$runLo, r$runHi), paste0(k, c("_ratio", "_lo", "_hi", "_runLo", "_runHi"))) }
    d6 <- g$detail[g$gate == "G6"]
    out <- c(gate("G2", "area"), gate("G4", "top10_median"), link(d6, "ignitions"), link(d6, "escapeRate"), link(d6, "areaPerEscape"))
    cbind(data.frame(ELF = g$ELF[1]), as.data.frame(as.list(out)))
  }))
}
