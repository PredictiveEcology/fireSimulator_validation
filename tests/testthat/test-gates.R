test_that("G1: within replicate noise passes, far outside fails, missing is NA", {
  expect_true(gateG1(c(100, 101, 99, 100.5), 100.2)$pass)
  expect_false(gateG1(c(100, 101, 99, 100.5), 130)$pass)
  expect_true(is.na(gateG1(NULL, 100)$pass))
})

test_that("G2: ratio of mean annual area to the observation, 0.25-4x", {
  sim <- data.frame(year = rep(1:2, each = 3), run = rep(1:3, 2), area50 = c(100, 200, 300, 100, 200, 300))
  obs <- data.frame(year = 1:2, area50 = c(200, 200))
  expect_true(gateG2(sim, obs)$pass); expect_equal(gateG2(sim, obs)$value, 1)
  expect_false(gateG2(sim, data.frame(year = 1:2, area50 = c(1e5, 1e5)))$pass)
  expect_false(gateG2(sim, data.frame(year = 1:2, area50 = c(10, 10)))$pass)
  expect_equal(gateG2(sim, data.frame(year = 1:2, area50 = c(800, 800)))$value, 0.25)
})

test_that("G3: the fit set must carry half of the area", {
  expect_true(gateG3(60, 100)$pass); expect_false(gateG3(40, 100)$pass); expect_true(is.na(gateG3(NA, 100)$pass))
})

test_that("G4: median of the 10 largest simulated fires within 10x of the observed; share over 100 kha reported", {
  g <- gateG4(data.frame(ha = c(100, 5e5, 1e5)), data.frame(ha = c(200, 60000)))
  expect_equal(g$value, 1e5 / 30100); expect_true(g$pass)   # medians of all (fewer than 10) fires
  expect_equal(g$share100kha, 5e5 / (5e5 + 1e5 + 100))
  expect_match(g$detail, "largest fire of the median run 500000 ha, largest observed 60000 ha")
  expect_false(gateG4(data.frame(ha = 7e5), data.frame(ha = 6e4))$pass)
  ## one huge observed fire no longer sets the target: the 10 largest do
  obs <- data.frame(ha = c(rep(1000, 9), 1e6)); sim <- data.frame(ha = rep(1000, 10))
  expect_equal(gateG4(sim, obs)$value, 1)
})

test_that("topFires is the median of the k largest", {
  expect_equal(topFires(c(1:20), 10), 15.5); expect_equal(topFires(c(5, 1), 10), 3); expect_equal(topFires(numeric(0)), 0)
})

test_that("G5: every kept p_c below p50", {
  expect_true(gateG5(c(0.2, 0.24), 0.25)$pass); expect_false(gateG5(c(0.2, 0.26), 0.25)$pass); expect_true(is.na(gateG5(0.2, NULL)$pass))
})

test_that("G6: each link within 0.25-4x and the failing link is named", {
  sim <- data.frame(year = rep(1:2, each = 2), run = rep(1:2, 2), ignitions = 100, escapes = 10, area50 = 1000)
  obs <- data.frame(year = 1:2, ignitions = c(100, 100), escapes = c(10, 10), area50 = c(1000, 1000))
  expect_true(gateG6(sim, obs)$pass)
  obs$escapes <- c(80, 80)                       # observed escape rate 8x the simulated one; area per escape is 1/8x
  g <- gateG6(sim, obs)
  expect_false(g$pass); expect_match(g$detail, "failing: escapeRate, areaPerEscape")
  obs2 <- data.frame(year = 1:2, ignitions = c(1000, 1000), escapes = c(100, 100), area50 = c(10000, 10000))
  expect_match(gateG6(sim, obs2)$detail, "failing: ignitions")
})

test_that("runGates flags an ELF whose G2 or G3 fails", {
  e <- toyElf("ok", perfect = TRUE)
  expect_true(attr(runGates(e), "meaningful"))
  e$fit$fitAreaHa <- 0.1 * e$fit$totalAreaHa
  g <- runGates(e); expect_false(attr(g, "meaningful")); expect_false(g$pass[g$gate == "G3"])
  expect_equal(g$gate, paste0("G", 1:6))
})

test_that("gates summarise runs by the median: one runaway run does not move them", {
  sim <- data.frame(year = rep(1:2, each = 3), run = rep(1:3, 2), area50 = c(100, 200, 300, 100, 200, 300))
  obs <- data.frame(year = 1:2, area50 = c(200, 200))
  runaway <- sim; runaway$area50[runaway$run == 3] <- 1e7
  expect_equal(gateG2(runaway, obs)$value, gateG2(sim, obs)$value)
  sizes <- data.frame(year = 1, run = c(1, 2, 3), ha = c(1000, 2000, 1e7))
  expect_equal(gateG4(sizes, data.frame(year = 1, ha = 2000))$value, 1)  # median run's top-10 median = 2000
  s6 <- data.frame(year = 1, run = 1:3, ignitions = 100, escapes = 10, area50 = c(1000, 1000, 1e7))
  o6 <- data.frame(year = 1, ignitions = 100, escapes = 10, area50 = 1000)
  expect_match(gateG6(s6, o6)$detail, "areaPerEscape 1 \\(")
  expect_equal(gateG1(c(1, 2, 3, 1000), 2.5)$value, 0)                    # median 2.5, not mean 251.5
})

## run-to-run uncertainty: gate values carry a bootstrap interval over runs (resampled within each year)
synthRuns <- function(nRuns, seed = 3, years = 1:6) {
  set.seed(seed)
  sy <- do.call(rbind, lapply(years, function(y) data.frame(year = y, run = seq_len(nRuns), ignitions = rpois(nRuns, 100),
          escapes = rpois(nRuns, 10), area50 = rgamma(nRuns, 2, 1 / 500))))
  ss <- do.call(rbind, lapply(seq_len(nrow(sy)), function(i) data.frame(year = sy$year[i], run = sy$run[i], ha = 50 * exp(rexp(5)))))
  list(simYears = sy, simSizes = ss,
       obsYears = data.frame(year = years, ignitions = 100, escapes = 10, area50 = 1000),
       obsSizes = data.frame(year = rep(years, each = 5), ha = 50 * exp(rexp(5 * length(years)))))
}

test_that("gates G2, G4 and G6 report lo and hi beside value, and the detail shows the interval", {
  d <- synthRuns(20)
  for (g in list(gateG2(d$simYears, d$obsYears), gateG4(d$simSizes, d$obsSizes), gateG6(d$simYears, d$obsYears))) {
    expect_true(all(c("lo", "hi") %in% names(g)))
    expect_true(g$lo <= g$value && g$value <= g$hi)
    expect_match(g$detail, "\\(\\d")
  }
  expect_equal(gateG2(d$simYears, d$obsYears, nBoot = 200, seed = 5), gateG2(d$simYears, d$obsYears, nBoot = 200, seed = 5))
  expect_true(is.na(gateG1(c(100, 101, 99), 100)$lo))
  g <- runGates(toyElf("t")); expect_true(all(c("lo", "hi") %in% names(g)))
  expect_true(all(is.finite(g$lo[g$gate %in% c("G2", "G4", "G6")])))
})

test_that("the interval narrows as runs increase, and the pass rule does not use it", {
  w <- function(n) { d <- synthRuns(n); g <- gateG2(d$simYears, d$obsYears); g$hi - g$lo }
  expect_lt(w(200), w(10))
  one <- synthRuns(1); g <- gateG2(one$simYears, one$obsYears)   # one run: nothing to resample
  expect_equal(g$lo, g$value); expect_equal(g$hi, g$value)
  d <- synthRuns(20); expect_equal(gateG2(d$simYears, d$obsYears, nBoot = 10)$pass, gateG2(d$simYears, d$obsYears, nBoot = 500)$pass)
})

test_that("runWave1 gates table has lo and hi", {
  r <- runWave1(list(toyElf("a")), nBoot = 20, nNull = 5, nPerm = 99)
  expect_true(all(c("lo", "hi") %in% names(r$gates)))
})
