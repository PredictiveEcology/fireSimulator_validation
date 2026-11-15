## Shared climate folder, and scoring per parameter set. Synthetic hindcast folders, no fit inputs.

## a CMD raster stack of `vals` (one layer per year of `years`) written to `dir`
writeToyClimate <- function(dir, elf, years, vals, variable = "CMD") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  g <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = "EPSG:3978", vals = 0)
  terra::writeRaster(terra::rast(lapply(vals, function(v) terra::setValues(g, v))) |> stats::setNames(paste0(variable, "_", years)),
                     file.path(dir, paste0(variable, "_historical_", elf, ".tif")))
  g
}

test_that("readHindcastCovariate finds the climate in a year folder (old layout)", {
  root <- withr::local_tempdir()
  g <- writeToyClimate(file.path(root, "9.9", 1986, 1, "climate"), "9.9", 1986:1987, c(5, 7))
  r <- readHindcastCovariate(root, "9.9", terra::setValues(g, 1))
  expect_equal(r$year, 1986:1987); expect_equal(r$value, c(5, 7))
})

test_that("readHindcastCovariate finds the climate in one shared folder per ELF (new layout), which wins over a year folder", {
  root <- withr::local_tempdir()
  g <- writeToyClimate(file.path(root, "9.9", "climate"), "9.9", 1986:1987, c(5, 7))
  dir.create(file.path(root, "9.9", 1986, 1), recursive = TRUE)          # a year folder with no climate copy
  r <- readHindcastCovariate(root, "9.9", terra::setValues(g, 1))
  expect_equal(r$value, c(5, 7))
  writeToyClimate(file.path(root, "9.9", 1986, 1, "climate"), "9.9", 1986:1987, c(50, 70))
  expect_equal(readHindcastCovariate(root, "9.9", terra::setValues(g, 1))$value, c(5, 7))
})

## One year, two parameter sets in one folder (paramSet column), 3 reps each; set 2 has twice the area of set 1
mkTwoSets <- function(root, elf, years = c(1986L, 1991L)) {
  for (y in years) {
    d <- file.path(root, elf, y, 1); dir.create(d, recursive = TRUE)
    reps <- data.frame(paramSet = rep(1:2, each = 3), rep = rep(1:3, 2), year = y, time = y, ignitions = 10, escapes = 3)
    write.csv(reps, file.path(d, "hindcast_reps.csv"), row.names = FALSE)
    z <- data.frame(paramSet = reps$paramSet, rep = reps$rep, fire_id = 1, ha = 100 * reps$paramSet * (1 + reps$rep))
    write.csv(z, file.path(d, "hindcast_sizes.csv"), row.names = FALSE)
  }
}

test_that("readHindcastSim keeps the paramSet column: a run is paramSet x rep", {
  root <- withr::local_tempdir(); mkTwoSets(root, "9.9", 1986L)
  r <- readHindcastSim(root, "9.9", 1986L, sets = 1L, minHa = toyMinHa)
  expect_equal(sort(unique(r$simYears$run)), c(1001, 1002, 1003, 2001, 2002, 2003))
  expect_equal(r$simYears$paramSet, rep(1:2, each = 3)); expect_equal(r$simSizes$paramSet, rep(1:2, each = 3))
})

test_that("scoreBySet: one row per paramSet x ELF, medians with single-run ranges; the set with more area has the bigger ratio", {
  e <- toyElf(nRuns = 3L, scored = c(1986L, 1991L))
  root <- withr::local_tempdir(); mkTwoSets(root, e$elf)
  s <- readHindcastSim(root, e$elf, c(1986L, 1991L), sets = 1L, minHa = toyMinHa)
  e$simYears <- s$simYears; e$simSizes <- s$simSizes
  out <- scoreBySet(e, nBoot = 20)
  expect_equal(nrow(out), 2L); expect_equal(out$ELF, rep(e$elf, 2)); expect_equal(out$paramSet, 1:2)
  expect_true(all(c("area_ratio", "area_runLo", "area_runHi", "top10_median_ratio", "ignitions_ratio", "escapeRate_ratio",
                    "areaPerEscape_ratio", "annualArea_skillN1", "annualArea_skillN2", "G2", "G4", "G6", "nRuns") %in% names(out)))
  expect_gt(out$area_ratio[2], out$area_ratio[1])
  expect_equal(out$area_ratio[1], gateG2(e$simYears[e$simYears$paramSet == 1, ], e$obsYears, nBoot = 20)$value)   # the module's own gate, on that set's runs
  expect_equal(out$nRuns, c(3L, 3L))
})
