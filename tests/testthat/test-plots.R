## The figures take the ELF input list of runWave1() and return ggplots. Toy ELF: see helper-source.R.

## a toy ELF with a record year of zero fires and a simulated run with no fire of 50 ha or more
toyPlotElf <- function() {
  e <- toyElf(seed = 3)
  e$obsYears[e$obsYears$year == 1991, c("ignitions", "escapes", "area50")] <- 0
  e$obsSizes <- e$obsSizes[e$obsSizes$year != 1991, ]
  e$simSizes <- e$simSizes[e$simSizes$run != 1, ]
  e$simYears[e$simYears$run == 1, c("escapes", "area50")] <- 0
  e
}
noWarn <- function(p) expect_no_warning(if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else ggplot_build(p))

test_that("every plot function returns a ggplot that builds without warnings, with a zero-fire year and a run without fires", {
  e <- toyPlotElf()
  for (f in list(plotAreaByYear, plotAnnualAreaDensity, plotFireSizes)) {
    p <- f(e, minHa = toyMinHa); expect_s3_class(p, "ggplot"); noWarn(p)
  }
  for (f in list(plotIgnitions, plotEscapes)) { p <- f(e); expect_s3_class(p, "ggplot"); noWarn(p) }
  expect_true(any(e$obsYears$area50 == 0))
  expect_s3_class(plotAreaByYear(e, label = "my fit", minHa = 100), "ggplot")
  expect_true("my fit" %in% levels(plotIgnitions(e, label = "my fit")$data$series))
})

test_that("the annual-area subtitle gives the share of years without fire for both series", {
  p <- plotAnnualAreaDensity(toyPlotElf(), minHa = toyMinHa)
  expect_match(p$labels$subtitle, "years with none: observed [0-9]+%, simulated [0-9]+%")
})

test_that("fire sizes are a density curve, escapes are bars", {
  e <- toyPlotElf()
  expect_s3_class(plotFireSizes(e, minHa = toyMinHa)$layers[[1]]$geom, "GeomRibbon")
  isCol <- function(p) any(vapply(p$layers, function(l) inherits(l$geom, "GeomCol"), TRUE))
  expect_false(isCol(plotFireSizes(e, minHa = toyMinHa))); expect_true(isCol(plotEscapes(e)))
})

test_that("an ELF whose simulation never burns still plots", {
  e <- toyPlotElf(); e$simSizes <- e$simSizes[0, , drop = FALSE]; e$simYears[, c("escapes", "area50")] <- 0
  for (f in list(plotAreaByYear, plotAnnualAreaDensity, plotFireSizes)) noWarn(f(e, minHa = toyMinHa))
  noWarn(plotEscapes(e))
})

test_that("summariseGates gives the five ratios of runWave1()'s gates and plotELFPanel puts them in the subtitle", {
  e <- toyPlotElf()
  g <- do.call(rbind, lapply(list(e, within(e, elf <- "toy2")), function(x) cbind(ELF = x$elf, runGates(x, nBoot = 20))))
  s <- summariseGates(g)
  expect_equal(s$ELF, c("toy", "toy2"))
  expect_equal(s$area_ratio[1], g$value[g$gate == "G2"][1])
  expect_equal(s$top10_median_runHi[1], g$runHi[g$gate == "G4"][1])
  expect_true(all(is.finite(s$escapeRate_ratio)))
  p <- plotELFPanel(e, s, minHa = toyMinHa)
  expect_s3_class(p, "patchwork"); noWarn(p)
  expect_equal(p$patches$annotation$title, "ELF toy")
  expect_match(p$patches$annotation$subtitle, "area [0-9]+\\.[0-9]{2} \\([0-9]+\\.[0-9]{2}-[0-9]+\\.[0-9]{2}\\).*area per escape [0-9]+\\.[0-9]{2}")
  expect_length(p$patches$plots, 4)
  noWarn(plotELFPanel(e, minHa = toyMinHa))
})

test_that("plotOverview orders ELFs by burn rate and builds without warnings", {
  m <- function(r) { x <- list(); for (k in names(gateRatios)) x <- c(x, setNames(list(r, r / 2, r * 2, r / 2, r * 2), paste0(k, c("_ratio", "_lo", "_hi", "_runLo", "_runHi")))); x }
  s <- do.call(rbind, Map(function(id, r) data.frame(ELF = id, m(r)), c("a", "b", "c"), c(2, 0.5, 1)))
  s$burnRate <- c(3, 0.1, 1)
  p <- plotOverview(s)
  expect_s3_class(p, "ggplot"); noWarn(p)
  expect_equal(levels(p$data$ELF), c("b", "c", "a"))
  expect_setequal(levels(p$data$measure), unname(gateRatios))
  s$area_runLo[1] <- 0; s$ignitions_runHi[2] <- Inf        # a zero or infinite whisker end is drawn at the limit
  noWarn(plotOverview(s))
})

test_that("plotOverview draws one line per arm within each ELF when the summary has an arm column", {
  m <- function(r) { x <- list(); for (k in names(gateRatios)) x <- c(x, setNames(list(r, r / 2, r * 2, r / 2, r * 2), paste0(k, c("_ratio", "_lo", "_hi", "_runLo", "_runHi")))); x }
  s <- do.call(rbind, Map(function(id, a, r) data.frame(ELF = id, arm = a, group = substr(a, 1, 1), m(r)),
                          c("a", "a", "b", "b", "b"), c("x1", "y1", "x1", "y1", "y2"), c(2, 0.5, 1, 1.5, 0.8)))
  s$burnRate <- c(3, 3, 0.1, 0.1, 0.1)
  p <- plotOverview(s)
  expect_s3_class(p, "ggplot"); noWarn(p)
  expect_equal(nrow(p$data), 5 * length(gateRatios))                       # every ELF x arm x ratio
  expect_equal(levels(p$data$arm), c("y2", "y1", "x1"))                    # first arm at the top
  expect_s3_class(p$facet, "FacetGrid")
  expect_equal(overviewSizeFor(5)$height, overviewSize$height)
  expect_gt(overviewSizeFor(60)$height, overviewSize$height)
})

test_that("Plots() with .plots = png writes the panel and the overview to disk", {
  skip_if_not_installed("SpaDES.core")
  e <- toyPlotElf(); g <- cbind(ELF = e$elf, runGates(e, nBoot = 20)); s <- summariseGates(g)
  d <- file.path(tempdir(), "valFigs"); dir.create(d)
  SpaDES.core::Plots(data = e, fn = plotELFPanel, summary = s, minHa = toyMinHa, types = "png", path = d, filename = "panel_toy",
                     ggsaveArgs = panelSize, .plotInitialTime = NA)
  SpaDES.core::Plots(data = s, fn = plotOverview, types = "png", path = d, filename = "overview",
                     ggsaveArgs = overviewSize, .plotInitialTime = NA)
  SpaDES.core::Plots(data = list(toy = elfPlots(e, minHa = toyMinHa)), fn = plotMetricPanel, summary = s, metric = "escapes", types = "png", path = d,
                     filename = "metric_escapes", ggsaveArgs = metricSize(1), .plotInitialTime = NA)
  expect_true(all(file.exists(file.path(d, c("panel_toy.png", "overview.png", "metric_escapes.png")))))
  expect_gt(min(file.size(file.path(d, c("panel_toy.png", "overview.png")))), 10000)
})

test_that("plotMetricPanel gives one pane per ELF in area-ratio order, for every metric, from plots built once", {
  e1 <- toyPlotElf(); e2 <- within(toyElf(seed = 5), elf <- "toy2")
  g <- do.call(rbind, lapply(list(e1, e2), function(x) cbind(ELF = x$elf, runGates(x, nBoot = 20))))
  s <- summariseGates(g)
  pl <- list(toy = elfPlots(e1, minHa = toyMinHa), toy2 = elfPlots(e2, minHa = toyMinHa))
  for (m in names(metrics)) {
    p <- plotMetricPanel(pl, s, m)
    expect_s3_class(p, "patchwork"); noWarn(p)
    expect_length(p$patches$plots, 1)
    expect_equal(p$patches$plots[[1]]$labels$title, s$ELF[which.max(s$area_ratio)])
    expect_match(p$patches$annotation$title, metrics[[m]], fixed = TRUE)
  }
  cap <- function(m) plotMetricPanel(pl, s, m)$patches$annotation$caption
  expect_false(grepl("50 ha|error bars|area by year|x: observed", cap("ignitions")))
  expect_match(cap("escapes"), "error bars")
  expect_match(plotELFPanel(e1, s, minHa = toyMinHa, plots = pl$toy)$patches$annotation$caption, "error bars.*x: observed|x: observed.*error bars")
  expect_error(plotMetricPanel(pl, s, "nope"))
  expect_equal(metricSize(14)$height, 4 * 5 + 2)
  f <- tempfile(fileext = ".png")                                    # 43 ELFs: 62 inches, over ggsave's 50-inch limit
  expect_no_error(do.call(ggplot2::ggsave, c(list(f, ggplot2::ggplot()), utils::modifyList(metricSize(43), list(dpi = 5)))))
  expect_equal(plotELFPanel(e1, s, minHa = toyMinHa, plots = pl$toy)$patches$annotation$title, "ELF toy")
})

test_that("elfOrder orders by burn rate, lowest first, an ELF without one at the start (last from the top)", {
  s <- data.frame(ELF = c("a", "b", "c", "d"), area_ratio = c(9, 1, 5, 2), burnRate = c(0.5, NA, 2, 0.1))
  expect_equal(elfOrder(s), c("b", "d", "a", "c"))
  expect_equal(elfOrder(s[, c("ELF", "area_ratio")]), c("a", "b", "c", "d"))        # no burnRate at all: input order
})
