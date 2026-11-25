## scoreClimateResponse(): slope, dry10 and wetHalf of the simulated annual area over the observed. Toy ELF: see helper-source.R.

## A toy ELF with a known climate response: observed area = exp(5 + b z) in each scored year; run r is `f(r)` times a response `b * s` (slope factor s)
climToy <- function(b = 0.5, s = 1, f = function(r) 1, nRuns = 5L, years = 1986:2024) {
  z <- as.numeric(scale(sin(seq_along(years) * 1.7) + seq_along(years) / 40))   # distinct values: no ties at the quantiles
  obs <- exp(5 + b * z)
  sim <- do.call(rbind, lapply(seq_len(nRuns), function(r) data.frame(year = years, run = r, ignitions = 1, escapes = 1, area50 = f(r) * exp(5 + b * s * z))))
  list(elf = "clim", simYears = sim, simSizes = data.frame(year = years, run = 1L, ha = 100),
       obsYears = data.frame(year = years, ignitions = 1, escapes = 1, area50 = obs), obsSizes = data.frame(year = years, ha = 100),
       covariate = data.frame(year = years, value = 100 + 20 * z), z = z, obs = obs)
}

test_that("a simulator with the observed response has ratios of 1; a doubled slope has a slope ratio of 2", {
  e <- climToy()
  cr <- scoreClimateResponse(e)
  expect_equal(cr$metric, c("slope", "dry10", "wetHalf")); expect_equal(cr$ratio, c(1, 1, 1)); expect_equal(cr$runLo, c(1, 1, 1)); expect_equal(cr$runHi, c(1, 1, 1))
  expect_equal(cr$observed[1], 0.5); expect_equal(unique(cr$nRuns), 5L); expect_equal(unique(cr$nYears), 39L)
  cr2 <- scoreClimateResponse(climToy(s = 2))
  expect_equal(cr2$ratio[1], 2)                                                # the covariate is standardised: slope = b
  expect_gt(cr2$ratio[2], 1); expect_lt(cr2$ratio[3], 1)                      # steeper: more area in the driest years, less in the wetter half
})

test_that("dry10 is the area in the driest 10% of years and wetHalf that in the wetter half, simulated over observed", {
  e <- climToy(f = function(r) 2)
  cr <- scoreClimateResponse(e)
  expect_equal(cr$ratio, c(1, 2, 2))                                           # a constant factor leaves the slope, doubles the areas
  z <- e$z; dry <- z >= quantile(z, 0.9); wet <- z <= median(z)
  expect_equal(sum(dry), 4L); expect_equal(sum(wet), 20L)
  expect_equal(cr$observed[2], sum(e$obs[dry])); expect_equal(cr$observed[3], sum(e$obs[wet]))
  expect_equal(cr$simMedian[2], 2 * sum(e$obs[dry]))
})

test_that("the ratio is the median of the single runs and the range their 2.5-97.5% quantiles", {
  e <- climToy(nRuns = 9L, f = function(r) r)                                  # run r has r times the observed area
  cr <- scoreClimateResponse(e)
  r <- 1:9
  expect_equal(cr$ratio[2:3], c(median(r), median(r))); expect_equal(cr$runLo[2], unname(quantile(r, 0.025))); expect_equal(cr$runHi[3], unname(quantile(r, 0.975)))
  expect_equal(cr$ratio[1], 1)
})

test_that("no covariate: no result; years without fire are left out of a slope, which is NA with fewer than 3", {
  e <- climToy(); e$covariate <- NULL
  expect_null(scoreClimateResponse(e)); expect_null(climateTable(list(e)))
  e <- climToy(); e$simYears$area50[e$simYears$run == 1 & e$simYears$year > 1987] <- 0
  cr <- scoreClimateResponse(e); expect_equal(cr$nRuns[1], 5L); expect_equal(cr$ratio[1], 1)   # the median of the other runs; the NA run is left out
  ## fewer than 3 years with fire: the run's slope is NA, and so is the ELF's when every run has so few
  e$simYears$area50[e$simYears$year > 1987] <- 0
  cr <- scoreClimateResponse(e); expect_true(is.na(cr$simMedian[1]) && is.na(cr$ratio[1])); expect_true(all(is.finite(cr$ratio[2:3])))
})

test_that("an observed slope of 0 gives an NA slope ratio, the area ratios stay", {
  e <- climToy(); e$obsYears$area50 <- 100                                     # constant observed area: slope exactly 0
  cr <- scoreClimateResponse(e)
  expect_equal(cr$observed[1], 0, tolerance = 1e-12); expect_true(is.na(cr$ratio[1])); expect_true(all(is.finite(cr$ratio[2:3])))
})

test_that("a run missing scored years does not error: the missing years count as no area", {
  e <- climToy(); miss <- e$simYears$run == 2 & e$simYears$year %in% c(1990, 2000)
  e$simYears <- e$simYears[!miss, ]
  cr <- scoreClimateResponse(e)
  expect_equal(cr$nRuns[1], 5L); expect_true(all(is.finite(cr$ratio)))
  expect_lte(cr$runLo[2], 1); expect_equal(cr$runHi[2], 1)                      # run 2 has less area in the driest and wetter years it lacks
  expect_true(all(c("obsRank", "nNonFinite") %in% names(cr)))
})

test_that("the single-run summary is bootRuns()'s: the range and the observed rank", {
  single <- rbind(a = c(0.5, 1, 2, 4), b = c(1, 1, 1, 1))
  rs <- runSummary(single)
  expect_equal(rs$runLo, unname(c(quantile(c(0.5, 1, 2, 4), 0.025), 1))); expect_equal(unname(rs$obsRank), c(0.375, 0.5))
  cr <- scoreClimateResponse(climToy(nRuns = 9L, f = function(r) r)); expect_equal(cr$obsRank[2], 1 / 9 * (0 + 0.5))   # one run (r = 1) at the observed, none below
})

test_that("pooled: the median of the ELF ratios with their range; the rows are in runWave1()'s results as test C", {
  a <- climToy(s = 2); a$elf <- "a"; b <- climToy(s = 3); b$elf <- "b"; c3 <- climToy(s = 4); c3$elf <- "c"
  ct <- climateTable(list(a, b, c3))
  p <- ct[ct$level == "pooled" & ct$metric == "slope", ]
  expect_equal(c(p$ratio, p$elfLo, p$elfHi, p$nELF), c(3, 2, 4, 3))
  rr <- climateResultRows(ct)
  expect_equal(unique(rr$test), "C"); expect_equal(unique(rr$null), c("slope", "dry10", "wetHalf")); expect_equal(nrow(rr), 12L)
  expect_equal(rr$value[rr$level == "pooled" & rr$null == "slope"], 3); expect_equal(rr$hi[rr$level == "pooled" & rr$null == "slope"], 4)
  pg <- climatePageRows(ct)
  expect_equal(pg[[1]]$elf, "all: median of ELFs"); expect_equal(pg[[1]]$slope[1:3], c(3, 2, 4)); expect_equal(pg[[2]]$obsSlope, 0.5)
})

test_that("runWave1() returns the climate table and its rows in results; the page has the table", {
  skip_if_not_installed("jsonlite")
  e <- toyElf(seed = 3, name = "a-x"); e2 <- toyElf(seed = 5, name = "b-y")
  r <- runWave1(list(e, e2), minHa = toyMinHa, nBoot = 10, nNull = 5, nPerm = 20)
  expect_equal(sort(unique(r$climate$ELF)), c("a-x", "all", "b-y"))
  expect_equal(sum(r$results$test == "C"), 9L)
  g <- do.call(rbind, lapply(list(e, e2), function(x) cbind(ELF = x$elf, runGates(x, nBoot = 20))))
  f <- tempfile(fileext = ".html")
  buildValidationPage(list(e, e2), g, f, testthat::test_path("..", "..", "pageTemplate.html"), minHa = toyMinHa, dpi = 30, bias = r$bias, climate = r$climate)
  h <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_false(grepl("/*CLIMATE*/", h, fixed = TRUE)); expect_match(h, "const CLIMATE = [{", fixed = TRUE); expect_match(h, "all: median of ELFs", fixed = TRUE)
})
