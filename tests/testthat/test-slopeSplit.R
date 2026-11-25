## scoreSlopeSplit(): the climate slope split into the slopes of log ignitions, log escape rate and log mean escaped-fire size. Toy ELF: climToy() of helper-source.R.

## Observed: 100 ignitions every year, escapes = 20 exp(0.3 z), 200 ha per escape. A run has the slopes `sl = c(ignitions, escape rate, size)` (log units per sd of z).
splitToy <- function(sl = c(0, 0.3, 0), obsSl = c(0, 0.3, 0), nRuns = 3L) {
  e <- climToy(nRuns = nRuns); z <- e$z
  rec <- function(s) { ign <- 100 * exp(s[1] * z); esc <- ign * 0.2 * exp(s[2] * z); data.frame(ignitions = ign, escapes = esc, area50 = esc * 200 * exp(s[3] * z)) }
  e$obsYears[, c("ignitions", "escapes", "area50")] <- rec(obsSl)
  e$simYears[, c("ignitions", "escapes", "area50")] <- do.call(rbind, rep(list(rec(sl)), nRuns))
  e
}
parts <- c("slopeIgnitions", "slopeEscapeRate", "slopeSize")

test_that("a toy where only the escape rate responds to climate: the escape component carries the whole difference, the others are 0", {
  ss <- scoreSlopeSplit(splitToy(sl = c(0, 0.1, 0), obsSl = c(0, 0.3, 0)))
  expect_equal(ss$metric, parts)
  expect_equal(ss$observed, c(0, 0.3, 0)); expect_equal(ss$simMedian, c(0, 0.1, 0)); expect_equal(ss$ratio, c(0, -0.2, 0))
  expect_equal(c(ss$runLo[2], ss$runHi[2]), c(-0.2, -0.2)); expect_equal(ss$obsRank[2], 1)   # every run's escape slope is below the observed
  expect_equal(unique(ss$nRuns), 3L); expect_equal(unique(ss$nYears), 39L)
})

test_that("the three parts add up to the total slope, observed, simulated and as differences", {
  e <- splitToy(sl = c(0.2, -0.1, 0.15), obsSl = c(0.05, 0.3, 0.1), nRuns = 1L)
  ss <- scoreSlopeSplit(e); cr <- scoreClimateResponse(e)
  expect_equal(sum(ss$observed), cr$observed[1]); expect_equal(sum(ss$simMedian), cr$simMedian[1]); expect_equal(sum(ss$ratio), cr$ratio[1])
  expect_equal(ss$ratio, c(0.15, -0.4, 0.05))
})

test_that("years without an escape are left out of all three parts, and the parts still add up to the total slope on the same years", {
  e <- splitToy(sl = c(0.2, -0.1, 0.15), obsSl = c(0.05, 0.3, 0.1), nRuns = 1L)
  e$simYears$ignitions <- e$simYears$ignitions * exp(0.3 * sin(3 * seq_len(nrow(e$simYears))))   # noise, so that dropping a year changes the slope
  before <- scoreSlopeSplit(e)$simMedian
  gone <- e$simYears$year %in% c(1990, 2000, 2010)
  e$simYears$escapes[gone] <- 0; e$simYears$area50[gone] <- 0                  # the run has no escape (and no area) in these years
  ss <- scoreSlopeSplit(e); cr <- scoreClimateResponse(e)
  expect_equal(sum(ss$simMedian), cr$simMedian[1]); expect_equal(sum(ss$ratio), cr$ratio[1])
  expect_true(all(is.finite(ss$ratio)))
  ## and not as if the years had been present
  expect_false(isTRUE(all.equal(ss$simMedian, before)))
  e$simYears$escapes[e$simYears$year > 1987] <- 0; e$simYears$area50[e$simYears$year > 1987] <- 0   # fewer than 3 years with an escape: NA, not an error
  expect_true(all(is.na(scoreSlopeSplit(e)$simMedian)))
})

test_that("pooled: the median of the ELF differences with their 2.5-97.5% range over ELFs; the rows are in results as test C with measure difference", {
  a <- splitToy(sl = c(0, 0.4, 0)); a$elf <- "a"; b <- splitToy(sl = c(0, 0.5, 0)); b$elf <- "b"; c3 <- splitToy(sl = c(0, 0.6, 0)); c3$elf <- "c"   # escape differences 0.1, 0.2, 0.3
  st <- slopeSplitTable(list(a, b, c3))
  p <- st[st$level == "pooled" & st$metric == "slopeEscapeRate", ]
  expect_equal(c(p$ratio, p$elfLo, p$elfHi, p$nELF), c(0.2, unname(quantile(c(0.1, 0.2, 0.3), c(0.025, 0.975))), 3))
  expect_equal(st$ratio[st$level == "pooled" & st$metric != "slopeEscapeRate"], c(0, 0))
  rr <- slopeSplitResultRows(st)
  expect_equal(unique(rr$test), "C"); expect_equal(unique(rr$null), parts); expect_equal(unique(rr$measure), "difference"); expect_equal(nrow(rr), 12L)
  expect_equal(rr$value[rr$level == "pooled" & rr$null == "slopeEscapeRate"], 0.2); expect_equal(rr$hi[rr$level == "pooled" & rr$null == "slopeEscapeRate"], p$elfHi)
  expect_null(slopeSplitTable(list({ d <- a; d$covariate <- NULL; d }))); expect_null(slopeSplitResultRows(NULL))
  pg <- slopeSplitPageRows(st)
  expect_equal(pg[[1]]$elf, "all: median of ELFs"); expect_equal(pg[[1]]$slopeEscapeRate[1], 0.2); expect_true(all(is.na(pg[[1]]$simslopeEscapeRate)))
  expect_equal(pg[[2]]$obsslopeEscapeRate, 0.3); expect_equal(pg[[2]]$simslopeEscapeRate[1], 0.4); expect_equal(pg[[2]]$slopeEscapeRate[1], 0.1)
})

test_that("runWave1() returns the split and its rows in results; the page has the table with a note that the parts add up", {
  skip_if_not_installed("jsonlite")
  e <- toyElf(seed = 3, name = "a-x"); e2 <- toyElf(seed = 5, name = "b-y")
  r <- runWave1(list(e, e2), minHa = toyMinHa, nBoot = 10, nNull = 5, nPerm = 20)
  expect_equal(sort(unique(r$slopeSplit$ELF)), c("a-x", "all", "b-y"))
  expect_equal(sum(r$results$test == "C" & r$results$null %in% parts), 9L); expect_equal(unique(r$results$measure[r$results$null %in% parts]), "difference")
  d <- r$slopeSplit[r$slopeSplit$ELF == "a-x", ]; tot <- r$climate[r$climate$ELF == "a-x" & r$climate$metric == "slope", ]
  expect_equal(sum(d$observed), tot$observed, tolerance = 1e-8)                # toy years with area have ignitions and escapes
  g <- do.call(rbind, lapply(list(e, e2), function(x) cbind(ELF = x$elf, runGates(x, nBoot = 20))))
  f <- tempfile(fileext = ".html")
  buildValidationPage(list(e, e2), g, f, testthat::test_path("..", "..", "pageTemplate.html"), minHa = toyMinHa, dpi = 30, bias = r$bias, climate = r$climate, slopeSplit = r$slopeSplit)
  h <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_false(grepl("/*SLOPESPLIT*/", h, fixed = TRUE)); expect_match(h, "const SLOPESPLIT = [{", fixed = TRUE)
  expect_match(h, "id=\"climateSplitRows\"", fixed = TRUE); expect_match(h, "data-s=\"slopeEscapeRate\"", fixed = TRUE)
  expect_match(h, "add up to the total slope difference", fixed = TRUE)
  expect_lt(regexpr("id=\"climateRows\"", h), regexpr("id=\"climateSplitRows\"", h))
  expect_lt(regexpr("data-s=\"slopeIgnitions\"", h), regexpr("data-s=\"slopeEscapeRate\"", h)); expect_lt(regexpr("data-s=\"slopeEscapeRate\"", h), regexpr("data-s=\"slopeSize\"", h))
})
