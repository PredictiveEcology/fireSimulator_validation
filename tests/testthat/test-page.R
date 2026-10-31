## runTablesToELF() and buildValidationPage(). Toy ELF: see helper-source.R.

test_that("runTablesToELF rebuilds the run-year table from fire sizes and ignitions, zero-filling run-years without fires", {
  e <- toyElf(seed = 2)
  ign <- e$simYears[, c("year", "run", "ignitions")]
  sz <- rbind(e$simSizes, data.frame(year = 1986L, run = 1L, ha = 10))                  # a fire under minHa is dropped
  r <- runTablesToELF("rt", sz, ign, e$obsYears, e$obsSizes, years = sort(unique(ign$year)), runs = 20L)
  key <- paste(r$simYears$year, r$simYears$run); ref <- paste(e$simYears$year, e$simYears$run)
  expect_equal(r$simYears[match(ref, key), c("ignitions", "escapes")], e$simYears[, c("ignitions", "escapes")], ignore_attr = TRUE)
  expect_equal(r$simYears$area50[match(ref, key)], e$simYears$area50)
  expect_true(all(r$simSizes$ha >= 50))
  z <- runTablesToELF("rt", sz[0, ], ign, e$obsYears, e$obsSizes)
  expect_true(all(z$simYears$escapes == 0 & z$simYears$area50 == 0))
})

test_that("buildValidationPage writes one file with the rows, the ranks and the images inside", {
  skip_if_not_installed("jsonlite")
  e1 <- toyElf(seed = 3, name = "a-x"); e2 <- toyElf(seed = 5, name = "b-y")
  g <- do.call(rbind, lapply(list(e1, e2), function(x) cbind(ELF = x$elf, runGates(x, nBoot = 20))))
  f <- tempfile(fileext = ".html")
  buildValidationPage(list(e1, e2), g, f, testthat::test_path("..", "..", "pageTemplate.html"), title = "Toy validation", subtitle = "toy run", dpi = 30)
  h <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_false(grepl("__[A-Z]+__|/\\*ROWS\\*/\\[\\]|/\\*IMG\\*/\\{\\}", h))
  expect_match(h, "<title>Toy validation</title>", fixed = TRUE); expect_match(h, "<h1>Toy validation</h1>", fixed = TRUE)
  expect_match(h, "toy run", fixed = TRUE); expect_match(h, "data:image/png;base64,", fixed = TRUE)
  expect_match(h, "\"panel_a-x\"", fixed = TRUE); expect_match(h, "\"metric_fireSizes\"", fixed = TRUE)
  rk <- observedRanks(e1)
  expect_true(all(rk >= 0 & rk <= 1)); expect_named(rk, c("area", "top10", "ign", "esc", "ape"))
})
