skillOf <- function(d, null) { x <- d[d$null == null, ]; skillScore(x$lossModel, x$lossNull) }

test_that("S1-S3: a perfect model gets skill near 1 against N1 and N2; a climatological model does not beat N1", {
  p <- toyElf(perfect = TRUE); c0 <- toyElf(perfect = FALSE)
  for (f in list(scoreS1, scoreS2, scoreS3)) {
    s <- f(p)$scores
    expect_gt(skillOf(s, "N1"), 0.95); expect_gt(skillOf(s, "N2"), 0.9)
    expect_equal(sort(unique(s$null)), c("N0", "N1", "N2"))
  }
  expect_lt(skillOf(scoreS1(c0)$scores, "N1"), 0.2)
  expect_lt(skillOf(scoreS3(c0)$scores, "N1"), 0.2)
})

test_that("S1 calibration: coverage and log score columns; a perfect model's intervals contain the observation", {
  r <- scoreS1(toyElf(perfect = TRUE))
  expect_true(all(c("cover50", "cover90", "logScore") %in% names(r$cal)))
  expect_true(all(r$cal$cover50)); expect_true(all(is.finite(r$cal$logScore)))
})

test_that("S2 drops years without observed ignitions", {
  e <- toyElf(perfect = TRUE); e$obsYears$ignitions[e$obsYears$year == 1991] <- 0; e$obsYears$escapes[e$obsYears$year == 1991] <- 0
  expect_false(1991L %in% scoreS2(e)$scores$year)
})

test_that("S3 uses fires of 50 ha or more only: it reads area50, not total area", {
  e <- toyElf(perfect = TRUE); s0 <- scoreS3(e)$scores
  e$simYears$ignitions <- e$simYears$ignitions * 100       # ignitions do not enter S3
  expect_equal(scoreS3(e)$scores, s0)
})

test_that("S4: perfect sizes beat the pooled-sizes null; sizes far too big are worse than N1; the largest-fire CRPS is 0 for a perfect model", {
  skip_if_not_installed("kSamples")
  p <- toyElf(perfect = TRUE)
  s <- scoreS4(p, nNull = 49)
  ad <- s$ad$scores; lg <- s$largest$scores
  expect_gt(skillOf(ad, "N1"), 0.5)
  expect_equal(unique(lg$lossModel), 0)
  expect_gt(skillOf(lg, "N1"), 0.95)
  b <- p; b$simSizes$ha <- b$simSizes$ha * 50
  expect_lt(skillOf(scoreS4(b, nNull = 49)$ad$scores, "N1"), 0)
})

test_that("scoreSizesVsNull keeps the scored year in the pool when leaveOut = FALSE", {
  skip_if_not_installed("kSamples")
  hist <- data.frame(year = c(rep(1:3, each = 5)), ha = c(rep(60, 5), rep(100, 5), rep(5e4, 5)))
  sim <- data.frame(year = 3, ha = rep(5e4, 5)); obs <- hist[hist$year == 3, ]
  a <- scoreSizesVsNull(sim, obs, hist, nNull = 99, nulls = "history", leaveOut = TRUE)$byYear
  b <- scoreSizesVsNull(sim, obs, hist, nNull = 99, nulls = "history", leaveOut = FALSE)$byYear
  expect_gt(a$adNullMedian, b$adNullMedian)           # the year's own sizes help the in-sample null
  expect_equal(a$adModel, b$adModel)
})

test_that("S5: a model that ranks the years right has rho near 1, and the difference to an index is reported", {
  p <- toyElf(perfect = TRUE); r <- scoreS5(p)
  expect_gt(r$rho, 0.99); expect_equal(r$rhoDiff, r$rho - r$rhoIndex)
  expect_equal(nrow(attr(r, "years")), 8)
})

test_that("the runner returns the tidy table, gates, per-ELF tables and the flagged ELFs", {
  skip_if_not_installed("kSamples")
  elfs <- list(toyElf("a", perfect = TRUE, seed = 1), toyElf("b", perfect = TRUE, seed = 2), toyElf("c", perfect = FALSE, seed = 3))
  elfs[[3]]$fit$fitAreaHa <- 0                                  # fails G3
  out <- runWave1(elfs, nBoot = 49, nNull = 29, nPerm = 99)
  expect_named(out, c("results", "gates", "perELF", "flagged"))
  expect_equal(out$flagged, "c")
  expect_equal(nrow(out$gates), 18)
  r <- out$results
  expect_named(r, c("test", "null", "level", "ELF", "measure", "value", "lo", "hi", "nELF", "nUnits", "verdict"))
  pooled <- r[r$level == "pooled" & r$measure == "skill", ]
  expect_setequal(unique(pooled$test), c("S1", "S2", "S3", "S4ad", "S4largest"))
  expect_true(all(pooled$nELF == 3L | pooled$nELF == 2L | pooled$nELF == 1L))
  expect_true(all(is.na(r$verdict[r$level == "ELF" | r$null %in% "N0"])))           # per-ELF numbers are reported, not judged
  s1 <- pooled[pooled$test == "S1" & pooled$null == "N1", ]; expect_true(s1$lo < s1$value && s1$value < s1$hi + 1e-12)
  expect_true(all(c("cover50", "cover90", "logScore") %in% r$measure))
  expect_true(all(c("rhoModel", "rhoIndex", "rhoDiff") %in% r$measure[r$test == "S5"]))
  expect_named(out$perELF, c("S1", "S2", "S3", "S4ad", "S4largest", "S5", "calibration"))
})
