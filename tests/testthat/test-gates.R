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

test_that("G4: largest simulated fire within 10x of the largest observed; share in fires over 100 kha reported", {
  g <- gateG4(data.frame(ha = c(100, 5e5, 1e5)), data.frame(ha = c(200, 60000)))
  expect_equal(g$value, 5e5 / 6e4); expect_true(g$pass)
  expect_equal(g$share100kha, 5e5 / (5e5 + 1e5 + 100))
  expect_false(gateG4(data.frame(ha = 7e5), data.frame(ha = 6e4))$pass)
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
  expect_equal(gateG4(sizes, data.frame(year = 1, ha = 2000))$value, 1)  # median run's largest = 2000
  s6 <- data.frame(year = 1, run = 1:3, ignitions = 100, escapes = 10, area50 = c(1000, 1000, 1e7))
  o6 <- data.frame(year = 1, ignitions = 100, escapes = 10, area50 = 1000)
  expect_match(gateG6(s6, o6)$detail, "areaPerEscape 1(;|$)")
  expect_equal(gateG1(c(1, 2, 3, 1000), 2.5)$value, 0)                    # median 2.5, not mean 251.5
})
