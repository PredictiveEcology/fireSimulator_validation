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

test_that("observedRanks compares the runs with the observed record of the scored years only", {
  rk <- observedRanks(toyElf(perfect = TRUE))                        # every run is the observed history: ties, rank 0.5
  expect_equal(unname(rk), rep(0.5, 5))
})

test_that("the page carries each ELF's burn rate, sorts by it by default and puts an ELF without one last", {
  skip_if_not_installed("jsonlite")
  e1 <- toyElf(seed = 3, name = "a-x"); e2 <- toyElf(seed = 5, name = "b-y"); e3 <- toyElf(seed = 7, name = "c-z")
  e1$burnRate <- 0.2; e2$burnRate <- 1.5                                         # e3 has none
  g <- do.call(rbind, lapply(list(e1, e2, e3), function(x) cbind(ELF = x$elf, runGates(x, nBoot = 20))))
  f <- tempfile(fileext = ".html")
  buildValidationPage(list(e1, e2, e3), g, f, testthat::test_path("..", "..", "pageTemplate.html"), dpi = 30)
  h <- paste(readLines(f, warn = FALSE), collapse = "\n")
  rows <- jsonlite::fromJSON(sub(";.*", "", sub(".*const ROWS = ", "", h)), simplifyVector = FALSE)
  expect_equal(lapply(rows, `[[`, "burnRate"), list(0.2, 1.5, NULL), ignore_attr = TRUE)
  expect_match(h, 'let sortK = "burnRate", sortD = -1', fixed = TRUE)
  expect_match(h, 'data-k="burnRate">Burned %/yr', fixed = TRUE)
})

test_that("yearCoverage is the share of scored years with observed inside the runs' range", {
  expect_equal(yearCoverage(toyElf(perfect = TRUE)), c(area = 1, ignitions = 1, escapes = 1))   # every run is the observed year
  e <- toyElf(seed = 4)
  e$simYears[c("ignitions", "escapes", "area50")] <- 0                                          # every run: no fire
  e$obsYears$ignitions[e$obsYears$year == 1986L] <- 0                                            # one of the 8 scored years observed 0
  cv <- yearCoverage(e)
  expect_equal(cv[["ignitions"]], mean(e$obsYears$ignitions[e$obsYears$year %in% scoredYears(e)] == 0))
  expect_equal(cv[["ignitions"]], 1 / 8)
})

test_that("buildValidationPage draws a coverage map per per-year metric with the caller's map function, none without one", {
  skip_if_not_installed("jsonlite")
  e1 <- toyElf(seed = 3, name = "a-x"); e2 <- toyElf(perfect = TRUE, name = "b-y")
  g <- do.call(rbind, lapply(list(e1, e2), function(x) cbind(ELF = x$elf, runGates(x, nBoot = 20))))
  got <- list()
  map <- function(values, title, file) {
    got[[length(got) + 1L]] <<- list(values = values, title = title)
    grDevices::png(file, width = 400, height = 320); graphics::plot.new(); grDevices::dev.off()
  }
  f <- tempfile(fileext = ".html")
  buildValidationPage(list(e1, e2), g, f, testthat::test_path("..", "..", "pageTemplate.html"), dpi = 30, map = map)
  h <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_length(got, length(coverageOfMetric))
  expect_equal(got[[1]]$values, c("a-x" = yearCoverage(e1)[["area"]], "b-y" = 1))
  expect_match(got[[4]]$title, "^Escapes per year: ")
  for (m in names(coverageOfMetric)) expect_match(h, paste0("\"map_", m, "\":\"data:image/png"), fixed = TRUE)
  expect_false(grepl("\"map_fireSizes\":", h, fixed = TRUE))
  f2 <- tempfile(fileext = ".html")
  buildValidationPage(list(e1, e2), g, f2, testthat::test_path("..", "..", "pageTemplate.html"), dpi = 30)
  expect_false(grepl("\"map_[A-Za-z]+\":", paste(readLines(f2, warn = FALSE), collapse = "\n")))
})
