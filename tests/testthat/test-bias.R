## scoreBias(): the median run over the observation, as numbers. Toy ELF: see helper-source.R.

## toyElf() with every simulated run scaled: area, ignitions, escapes by `f` (sizes too, so the 10 largest scale with it)
scaleSim <- function(e, f) {
  e$simYears$area50 <- e$simYears$area50 * f; e$simYears$ignitions <- e$simYears$ignitions * f; e$simYears$escapes <- e$simYears$escapes * f
  e$simSizes$ha <- e$simSizes$ha * f
  e
}
bq <- function(b, q) b[b$quantity == q, ]

test_that("a simulator with twice the observed area has a ratio of 2 and the observation below every run; a perfect one has 1", {
  e <- toyElf(seed = 4, perfect = TRUE)
  b1 <- scoreBias(e, nBoot = 20)
  expect_equal(b1$quantity, c("area", "top10", "ignitions", "escapes"))
  expect_equal(b1$ratio, rep(1, 4)); expect_equal(b1$log2ratio, rep(0, 4)); expect_equal(b1$per1000yr[b1$quantity != "top10"], rep(0, 3))
  e2 <- scaleSim(e, 2); b2 <- scoreBias(e2, nBoot = 20)
  expect_equal(b2$ratio, rep(2, 4)); expect_equal(b2$log2ratio, rep(1, 4))
  expect_equal(b2$obsRank, rep(0, 4)); expect_equal(b2$runLo, rep(2, 4)); expect_equal(b2$runHi, rep(2, 4)); expect_equal(unique(b2$nRuns), 20L)
  ## per1000yr: (ratio - 1) x the observed annual mean x 1000, in the units of the quantity
  yrs <- unique(e$simYears$year); o <- e$obsYears[e$obsYears$year %in% yrs, ]
  expect_equal(bq(b2, "area")$per1000yr, mean(o$area50) * 1000)
  expect_equal(bq(b2, "ignitions")$per1000yr, mean(o$ignitions) * 1000); expect_equal(bq(b2, "escapes")$per1000yr, mean(o$escapes) * 1000)
  expect_equal(bq(b2, "area")$unit, "ha per 1000 yr"); expect_true(is.na(bq(b2, "top10")$per1000yr))
})

test_that("a halved simulator has a ratio of 0.5 and the observation above every run; ignitions and escapes are separate ratios", {
  e <- toyElf(seed = 4, perfect = TRUE)
  e$simYears$ignitions <- e$simYears$ignitions * 3                       # only ignitions off
  b <- scoreBias(e, nBoot = 20)
  expect_equal(bq(b, "ignitions")$ratio, 3); expect_equal(bq(b, "escapes")$ratio, 1); expect_equal(bq(b, "area")$ratio, 1)
  b <- scoreBias(scaleSim(toyElf(seed = 4, perfect = TRUE), 0.5), nBoot = 20)
  expect_equal(b$ratio, rep(0.5, 4)); expect_equal(b$obsRank, rep(1, 4)); expect_true(all(b$per1000yr[b$quantity != "top10"] < 0))
})

test_that("the ratio of the 10 largest fires is G4's value, the area ratio G2's, and the intervals are G4's and G2's", {
  e <- toyElf(seed = 6)
  b <- scoreBias(e, nBoot = 30, seed = 2); g <- runGates(e, nBoot = 30, seed = 2)
  for (cmp in list(c("top10", "G4"), c("area", "G2"))) {
    x <- bq(b, cmp[1]); y <- g[g$gate == cmp[2], ]
    expect_equal(x$ratio, y$value); expect_equal(c(x$lo, x$hi, x$runLo, x$runHi, x$obsRank), c(y$lo, y$hi, y$runLo, y$runHi, y$obsRank))
  }
})

test_that("pooled: the median of the ELF ratios with their range, and the ratio of sums, which differ when ELFs differ in size", {
  big <- toyElf(seed = 4, perfect = TRUE, name = "big"); big$obsYears$area50 <- big$obsYears$area50 * 10     # observed area 10x, simulated as before
  big$simYears$area50 <- big$simYears$area50 * 10
  a <- scaleSim(toyElf(seed = 5, perfect = TRUE, name = "a"), 2); c3 <- scaleSim(toyElf(seed = 7, perfect = TRUE, name = "c"), 4)
  tb <- biasTable(list(big, a, c3), nBoot = 20)
  expect_equal(unique(tb$level), c("ELF", "pooled")); expect_equal(sum(tb$level == "pooled"), 4L)
  p <- bq(tb[tb$level == "pooled", ], "area")
  expect_equal(p$ratio, 2); expect_equal(c(p$elfLo, p$elfHi), c(1, 4))           # median of 1, 2, 4
  e <- bq(tb[tb$level == "ELF", ], "area")
  expect_equal(p$ratioOfSums, sum(e$simMedian) / sum(e$observed)); expect_false(isTRUE(all.equal(p$ratioOfSums, p$ratio)))
  expect_equal(p$per1000yr, sum(e$per1000yr)); expect_true(is.na(bq(tb[tb$level == "pooled", ], "top10")$ratioOfSums))
})

test_that("runWave1 returns the bias table and its rows in the long results table", {
  skip_if_not_installed("fireSenseUtils")
  e <- scaleSim(toyElf(seed = 3), 1)
  r <- runWave1(list(e), minHa = toyMinHa, nBoot = 20, nNull = 19, nPerm = 99)
  expect_equal(r$bias, biasTable(list(e), nBoot = 20))
  bb <- r$results[r$results$test == "B", ]
  expect_setequal(bb$measure, c("ratio", "log2ratio")); expect_setequal(bb$null, names(biasQuantities)); expect_setequal(bb$level, c("ELF", "pooled"))
  x <- bb[bb$level == "ELF" & bb$measure == "ratio" & bb$null == "area", ]
  expect_equal(x$value, bq(r$bias, "area")$ratio[1]); a <- bq(r$bias, "area")[1, ]; expect_equal(c(x$lo, x$hi), c(a$runLo, a$runHi))
  expect_equal(r$results$value[r$results$test == "B" & r$results$measure == "log2ratio"], log2(r$results$value[r$results$test == "B" & r$results$measure == "ratio"]))
  expect_named(r$results, c("test", "null", "level", "ELF", "measure", "value", "lo", "hi", "nELF", "nUnits", "verdict"))
})

test_that("the page carries a sortable bias table: every cell has data-v, the pooled rows are there, and without bias there is no table", {
  skip_if_not_installed("jsonlite")
  e1 <- toyElf(seed = 3, name = "a-x"); e2 <- scaleSim(toyElf(seed = 5, name = "b-y"), 2)
  g <- do.call(rbind, lapply(list(e1, e2), function(x) cbind(ELF = x$elf, runGates(x, nBoot = 20))))
  tb <- biasTable(list(e1, e2), nBoot = 20)
  f <- tempfile(fileext = ".html")
  buildValidationPage(list(e1, e2), g, f, testthat::test_path("..", "..", "pageTemplate.html"), minHa = toyMinHa, dpi = 30, bias = tb)
  h <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_false(grepl("/\\*BIAS\\*/\\[\\]", h)); expect_match(h, "all: median of ELFs", fixed = TRUE); expect_match(h, "all: ratio of sums", fixed = TRUE)
  expect_match(h, "data-v=", fixed = TRUE); expect_match(h, "data-b=\"areaDiff\"", fixed = TRUE)
  rows <- biasPageRows(tb)
  expect_equal(length(rows), 4L); expect_equal(rows[[2]]$area[1], round(bq(tb, "area")$ratio[2], 3))
  buildValidationPage(list(e1, e2), g, f, testthat::test_path("..", "..", "pageTemplate.html"), minHa = toyMinHa, dpi = 30)
  expect_match(paste(readLines(f, warn = FALSE), collapse = "\n"), "const BIAS = [];", fixed = TRUE)
})

test_that("an ELF with no observed fire gives a ratio of 1 (none simulated) or Inf (some simulated), never NaN, for the 10 largest as for the totals", {
  e <- toyElf(seed = 4, perfect = TRUE); e$obsSizes <- e$obsSizes[integer(0), , drop = FALSE]
  e$obsYears$area50 <- 0; e$obsYears$escapes <- 0
  b <- scoreBias(e, nBoot = 20)
  expect_false(anyNA(b$ratio)); expect_equal(bq(b, "top10")$ratio, Inf); expect_equal(bq(b, "area")$ratio, Inf)
  e$simSizes <- e$simSizes[integer(0), , drop = FALSE]; e$simYears$area50 <- 0; e$simYears$escapes <- 0
  b <- scoreBias(e, nBoot = 20)
  expect_equal(bq(b, "top10")$ratio, 1); expect_equal(bq(b, "area")$ratio, 1); expect_equal(bq(b, "escapes")$ratio, 1)
})

test_that("pooling reports how many ELFs are in the median and how many were left out for a ratio that is not finite", {
  ok <- scoreBias(scaleSim(toyElf(seed = 5, perfect = TRUE, name = "a"), 2), nBoot = 20)
  z <- toyElf(seed = 4, perfect = TRUE, name = "z"); z$obsSizes <- z$obsSizes[integer(0), , drop = FALSE]
  bad <- scoreBias(z, nBoot = 20)
  p <- poolBias(rbind(ok, bad)); a <- bq(p, "top10")
  expect_equal(a$nELF, 2L); expect_equal(a$nNonFinite, 1L); expect_equal(a$ratio, bq(ok, "top10")$ratio)
  expect_equal(poolBias(ok)$nNonFinite, rep(0L, 4)); expect_equal(unique(poolBias(ok)$nELF), 1L)
})
