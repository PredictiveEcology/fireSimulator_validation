mkSet <- function(root, elf, y, s, nRep = 2) {
  d <- file.path(root, elf, y, s); dir.create(d, recursive = TRUE)
  write.csv(data.frame(rep = 1:nRep, year = y, time = y, ignitions = c(10, 20)[1:nRep], escapes = c(7, 8)[1:nRep], areaBurnedHa = c(300, 600)[1:nRep], nFires50 = c(1, 2)[1:nRep]), file.path(d, "hindcast_reps.csv"), row.names = FALSE)
  write.csv(data.frame(rep = c(1, 2, 2), fire_id = 1:3, ha = c(300, 100, 500)), file.path(d, "hindcast_sizes.csv"), row.names = FALSE)
  write.csv(data.frame(year = y, ignitions = 5, escapes = 3, areaBurnedHa = 100, nFires50 = 1), file.path(d, "hindcast_observed.csv"), row.names = FALSE)
  write.csv(data.frame(cause = "natural", ha = 100), file.path(d, "hindcast_observedSizes.csv"), row.names = FALSE)
}

test_that("hindcastFinishedYears needs every file of every set", {
  root <- withr::local_tempdir()
  for (s in 1:2) mkSet(root, "9.9", 1986, s)
  for (s in 1:2) mkSet(root, "9.9", 1991, s)
  file.remove(file.path(root, "9.9", 1991, 2, "hindcast_sizes.csv"))
  expect_equal(hindcastFinishedYears(root, "9.9", c(1986L, 1991L, 1996L), sets = 1:2), 1986L)
})

test_that("readHindcastSim: runs are set x rep, area50 sums the fires of 50 ha or more", {
  root <- withr::local_tempdir()
  for (s in 1:2) mkSet(root, "9.9", 1986, s)
  r <- readHindcastSim(root, "9.9", 1986L, sets = 1:2)
  expect_equal(nrow(r$simYears), 4); expect_equal(sort(unique(r$simYears$escapes)), c(1, 2)); expect_equal(sort(r$simYears$run), c(1001, 1002, 2001, 2002))
  expect_equal(r$simYears$area50[r$simYears$run == 1002], 500 + 100)
  expect_equal(r$simYears$area50[r$simYears$run == 1001], 300)
  expect_equal(nrow(r$simSizes), 6)
})

test_that("readHindcastObserved: scored years come from the hindcast's own files, area50 from natural fires of 50 ha or more", {
  root <- withr::local_tempdir(); mkSet(root, "9.9", 1986, 1)
  rec <- list(counts = data.frame(year = 1986:1987, ignitions = c(99, 7), escapes = c(9, 9), areaBurnedHa = 0, nFires50 = c(9, 2)),
              sizes = data.frame(year = c(1986, 1987, 1987, 1987, 1988), cause = c("natural", "natural", "natural", "allCause", NA), ha = c(5000, 80, 10, 400, NA)))
  o <- readHindcastObserved(rec, root, "9.9", 1986L)
  expect_equal(o$obsYears$ignitions, c(5, 7)); expect_equal(o$obsYears$area50, c(100, 80))
  expect_equal(o$obsYears$escapes, c(1, 2))                  # fires of 50 ha or more, not the ignition-point escapes
  expect_equal(o$nbacAreaHa, 400); expect_equal(o$fitAreaHa, 400)
})

test_that("readLedgerRow reads p_c, p50 and the objective", {
  f <- withr::local_tempfile(fileext = ".rds")
  sp <- data.frame(polygonID = c("1.1", "9.9")); sp$params <- list(data.frame(maxAsymptote = 0.1), data.frame(maxAsymptote = c(0.2, 0.21)))
  sp$objFunVal <- list(1, c(100, 101)); sp$powerLawPc50_spread <- c(0.3, 0.25)
  saveRDS(sp, f)
  r <- readLedgerRow(f, "9.9")
  expect_equal(r$pc, c(0.2, 0.21)); expect_equal(r$p50, 0.25); expect_equal(r$objective, c(100, 101))
})

test_that("the module file hands its input to the runner and nothing else", {
  exprs <- parse(testthat::test_path("..", "..", "fireSimulator_validation.R"), keep.source = FALSE)
  txt <- paste(deparse(exprs[[length(exprs)]]), collapse = " ")
  expect_match(txt, "runWave1\\(sim\\$validationELFs")
  expect_false(grepl("scheduleEvent", paste(vapply(exprs, function(e) paste(deparse(e), collapse = " "), ""), collapse = " ")))
})

test_that("observedRecord counts observed area on the burnable core, without a flammableRTM warning, and rebuilds an old cache", {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = "EPSG:3978")
  core <- terra::setValues(r, ifelse(terra::colFromCell(r, seq_len(terra::ncell(r))) == 1, 1, NA))   # burnable cells only
  pts <- terra::vect(cbind(500, 5500), crs = "EPSG:3978"); pts$YEAR <- 2000L; pts$CAUSE <- "N"; pts$SIZE_HA <- 10
  poly <- terra::vect("POLYGON ((500 500, 1500 500, 1500 1500, 500 1500, 500 500))", crs = "EPSG:3978")
  poly$YEAR <- 2000L; poly$FIRECAUS <- "Natural"                                                   # 100 ha, 50 ha in column 1
  testthat::local_mocked_bindings(load_nfdb_points = function(...) pts, load_nbac_polys = function(...) poly,
                                  .package = "fireregimetools")
  cacheFile <- file.path(withr::local_tempdir(), "9.9.rds")
  saveRDS("a cache of the whole-polygon count", cacheFile)
  expect_no_warning(rec <- observedRecord(core, "nfdb", "nbac", years = 2000L, cacheFile = cacheFile))
  expect_equal(rec$counts$areaBurnedHaAllCause, 50, tolerance = 0.03)
  expect_true(file.exists(sub("\\.rds$", "_burnable.rds", cacheFile)))
  expect_identical(readRDS(cacheFile), "a cache of the whole-polygon count")                       # the old file is left alone
})
