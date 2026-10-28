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
  for (f in list(plotAreaByYear, plotAnnualAreaDensity, plotFireSizes, plotIgnitions, plotEscapes)) {
    p <- f(e); expect_s3_class(p, "ggplot"); noWarn(p)
  }
  expect_true(any(e$obsYears$area50 == 0))
  expect_s3_class(plotAreaByYear(e, label = "my fit", minHa = 100), "ggplot")
  expect_true("my fit" %in% levels(plotIgnitions(e, label = "my fit")$data$series))
})

test_that("the annual-area subtitle gives the share of years without fire for both series", {
  p <- plotAnnualAreaDensity(toyPlotElf())
  expect_match(p$labels$subtitle, "years with none: observed [0-9]+%, simulated [0-9]+%")
})

test_that("fire sizes are a density curve, escapes are bars", {
  e <- toyPlotElf()
  expect_s3_class(plotFireSizes(e)$layers[[1]]$geom, "GeomRibbon")
  isCol <- function(p) any(vapply(p$layers, function(l) inherits(l$geom, "GeomCol"), TRUE))
  expect_false(isCol(plotFireSizes(e))); expect_true(isCol(plotEscapes(e)))
})

test_that("an ELF whose simulation never burns still plots", {
  e <- toyPlotElf(); e$simSizes <- e$simSizes[0, , drop = FALSE]; e$simYears[, c("escapes", "area50")] <- 0
  for (f in list(plotAreaByYear, plotAnnualAreaDensity, plotFireSizes, plotEscapes)) noWarn(f(e))
})

test_that("summariseGates gives the five ratios of runWave1()'s gates and plotELFPanel puts them in the subtitle", {
  e <- toyPlotElf()
  g <- do.call(rbind, lapply(list(e, within(e, elf <- "toy2")), function(x) cbind(ELF = x$elf, runGates(x, nBoot = 20))))
  s <- summariseGates(g)
  expect_equal(s$ELF, c("toy", "toy2"))
  expect_equal(s$area_ratio[1], g$value[g$gate == "G2"][1])
  expect_equal(s$top10_median_runHi[1], g$runHi[g$gate == "G4"][1])
  expect_true(all(is.finite(s$escapeRate_ratio)))
  p <- plotELFPanel(e, s)
  expect_s3_class(p, "patchwork"); noWarn(p)
  expect_equal(p$patches$annotation$title, "ELF toy")
  expect_match(p$patches$annotation$subtitle, "area [0-9]+\\.[0-9]{2} \\([0-9]+\\.[0-9]{2}-[0-9]+\\.[0-9]{2}\\).*area per escape [0-9]+\\.[0-9]{2}")
  expect_length(p$patches$plots, 4)
  noWarn(plotELFPanel(e))
})

test_that("plotOverview orders ELFs by area ratio and builds without warnings", {
  m <- function(r) { x <- list(); for (k in names(gateRatios)) x <- c(x, setNames(list(r, r / 2, r * 2, r / 2, r * 2), paste0(k, c("_ratio", "_lo", "_hi", "_runLo", "_runHi")))); x }
  s <- do.call(rbind, Map(function(id, r) data.frame(ELF = id, m(r)), c("a", "b", "c"), c(2, 0.5, 1)))
  p <- plotOverview(s)
  expect_s3_class(p, "ggplot"); noWarn(p)
  expect_equal(levels(p$data$ELF), c("b", "c", "a"))
  expect_setequal(levels(p$data$measure), unname(gateRatios))
  s$area_runLo[1] <- 0; s$ignitions_runHi[2] <- Inf        # a zero or infinite whisker end is drawn at the limit
  noWarn(plotOverview(s))
})

test_that("Plots() with .plots = png writes the panel and the overview to disk", {
  skip_if_not_installed("SpaDES.core")
  e <- toyPlotElf(); g <- cbind(ELF = e$elf, runGates(e, nBoot = 20)); s <- summariseGates(g)
  d <- file.path(tempdir(), "valFigs"); dir.create(d)
  SpaDES.core::Plots(data = e, fn = plotELFPanel, summary = s, types = "png", path = d, filename = "panel_toy",
                     ggsaveArgs = panelSize, .plotInitialTime = NA)
  SpaDES.core::Plots(data = s, fn = plotOverview, types = "png", path = d, filename = "overview",
                     ggsaveArgs = overviewSize, .plotInitialTime = NA)
  SpaDES.core::Plots(data = list(toy = elfPlots(e)), fn = plotMetricPanel, summary = s, metric = "escapes", types = "png", path = d,
                     filename = "metric_escapes", ggsaveArgs = metricSize(1), .plotInitialTime = NA)
  expect_true(all(file.exists(file.path(d, c("panel_toy.png", "overview.png", "metric_escapes.png")))))
  expect_gt(min(file.size(file.path(d, c("panel_toy.png", "overview.png")))), 10000)
})

test_that("plotMetricPanel gives one pane per ELF in area-ratio order, for every metric, from plots built once", {
  e1 <- toyPlotElf(); e2 <- within(toyElf(seed = 5), elf <- "toy2")
  g <- do.call(rbind, lapply(list(e1, e2), function(x) cbind(ELF = x$elf, runGates(x, nBoot = 20))))
  s <- summariseGates(g)
  pl <- list(toy = elfPlots(e1), toy2 = elfPlots(e2))
  for (m in names(metrics)) {
    p <- plotMetricPanel(pl, s, m)
    expect_s3_class(p, "patchwork"); noWarn(p)
    expect_length(p$patches$plots, 1)
    expect_equal(p$patches$plots[[1]]$labels$title, s$ELF[which.min(s$area_ratio)])
    expect_match(p$patches$annotation$title, metrics[[m]], fixed = TRUE)
  }
  expect_error(plotMetricPanel(pl, s, "nope"))
  expect_equal(metricSize(14)$height, 4 * 5 + 2)
  expect_equal(plotELFPanel(e1, s, plots = pl$toy)$patches$annotation$title, "ELF toy")
})
