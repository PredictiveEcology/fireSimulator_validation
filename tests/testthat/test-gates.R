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
    expect_match(g$detail, "\\(median-run 95%: \\d")
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

## spread of single runs: each gate's ratio computed for each run on its own, 2.5% and 97.5% across runs
test_that("gates G2, G4 and G6 report runLo and runHi, the spread of single runs, beside lo and hi", {
  d <- synthRuns(30)
  for (g in list(gateG2(d$simYears, d$obsYears), gateG4(d$simSizes, d$obsSizes), gateG6(d$simYears, d$obsYears))) {
    expect_true(all(c("runLo", "runHi") %in% names(g)))
    expect_true(g$runLo <= g$lo && g$hi <= g$runHi)          # single runs vary more than the median run
    expect_match(g$detail, "median-run 95%: [0-9.e+-]+-[0-9.e+-]+; single runs 95%: [0-9.e+-]+-[0-9.e+-]+")
  }
  expect_true(is.na(gateG1(c(100, 101, 99), 100)$runLo))
  g <- runGates(toyElf("t")); expect_true(all(c("runLo", "runHi") %in% names(g)))
  expect_true(all(is.finite(g$runLo[g$gate %in% c("G2", "G4", "G6")])))
  expect_true(all(c("runLo", "runHi") %in% names(runWave1(list(toyElf("a")), nBoot = 20, nNull = 5, nPerm = 99)$gates)))
})

test_that("the single-run range is each run's own ratio and widens with between-run variance", {
  sim <- data.frame(year = rep(1:2, each = 4), run = rep(1:4, 2), area50 = rep(c(100, 200, 300, 400), 2))
  obs <- data.frame(year = 1:2, area50 = c(250, 250))
  g <- gateG2(sim, obs)                                        # run totals 200, 400, 600, 800 over observed 500
  expect_equal(g$runLo, unname(quantile(c(200, 400, 600, 800) / 500, 0.025)))
  expect_equal(g$runHi, unname(quantile(c(200, 400, 600, 800) / 500, 0.975)))
  spread <- function(sd) {
    set.seed(2); s <- do.call(rbind, lapply(1:6, function(y) data.frame(year = y, run = 1:50, area50 = 500 * exp(rnorm(50, 0, sd)))))
    g <- gateG2(s, data.frame(year = 1:6, area50 = 500)); c(run = g$runHi - g$runLo, med = g$hi - g$lo)
  }
  expect_gt(spread(1)[["run"]], spread(0.2)[["run"]])
  expect_gt(spread(1)[["run"]], spread(1)[["med"]])
  one <- synthRuns(1); g <- gateG2(one$simYears, one$obsYears); expect_equal(g$runLo, g$value); expect_equal(g$runHi, g$value)
})

## where the observed record ranks among the single runs (a ratio sim/obs, so observed = 1)
test_that("bootRuns obsRank is 0 when every run is above the observed, 1 when all are below, 0.5 with ties at 1", {
  mk <- function(v) data.frame(year = 1, run = seq_along(v), x = v)
  stat <- function(s) c(ratio = stats::median(s$x))
  rk <- function(v) bootRuns(mk(v), stat, nBoot = 10)
  expect_equal(rk(c(2, 3, 4))$obsRank, c(ratio = 0))
  expect_equal(rk(c(0.2, 0.3, 0.4))$obsRank, c(ratio = 1))
  expect_equal(rk(c(1, 1, 1, 1))$obsRank, c(ratio = 0.5))
  expect_equal(rk(c(0.5, 1, 2, 3))$obsRank, c(ratio = 0.375))   # (1 + 0.5) / 4
  b <- rk(c(0.5, 1, 2, 3))
  expect_equal(b$nRuns, 4L); expect_equal(dim(b$single), c(1L, 4L)); expect_equal(as.vector(b$single), c(0.5, 1, 2, 3))
})

test_that("G2 and G4 rows carry obsRank and nRuns; G6 and the other gates do not", {
  d <- synthRuns(30)
  g2 <- gateG2(d$simYears, d$obsYears); g4 <- gateG4(d$simSizes, d$obsSizes)
  for (g in list(g2, g4)) { expect_true(g$obsRank >= 0 && g$obsRank <= 1); expect_equal(g$nRuns, 30L) }
  expect_equal(g2$obsRank, bootRuns(d$simYears, function(s) stats::median(perRun(s, "area50")) / sum(d$obsYears$area50), 10)$obsRank[[1]])
  expect_true(all(is.na(c(gateG1(c(100, 101, 99), 100)$obsRank, gateG3(5, 10)$obsRank, gateG5(0.1, 0.2)$obsRank))))
  all <- runGates(toyElf("t"))
  expect_true(all(c("obsRank", "nRuns") %in% names(all)))
  expect_true(all(is.finite(all$obsRank[all$gate %in% c("G2", "G4")])))
})

test_that("detail carries the observed rank, after the unchanged interval text, for G2, G4 and every G6 link", {
  d <- synthRuns(30)
  old <- "ratio [0-9.e+-]+ \\(median-run 95%: [0-9.e+-]+-[0-9.e+-]+; single runs 95%: [0-9.e+-]+-[0-9.e+-]+\\)"
  for (g in list(gateG2(d$simYears, d$obsYears), gateG4(d$simSizes, d$obsSizes))) {
    expect_match(g$detail, old)
    expect_match(g$detail, "\\); observed above [0-9]+% of 30 runs")
  }
  g6 <- gateG6(d$simYears, d$obsYears)$detail
  for (k in c("ignitions", "escapeRate", "areaPerEscape"))
    expect_match(g6, paste0(k, " [0-9.e+-]+ \\(median-run 95%: [0-9.e+-]+-[0-9.e+-]+; single runs 95%: [0-9.e+-]+-[0-9.e+-]+\\); observed above [0-9]+% of 30 runs"))
  expect_equal(lengths(regmatches(g6, gregexpr("observed above", g6))), 3L)
})
