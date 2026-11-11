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
  r <- readHindcastSim(root, "9.9", 1986L, sets = 1:2, minHa = toyMinHa)
  expect_equal(nrow(r$simYears), 4); expect_equal(sort(unique(r$simYears$escapes)), c(1, 2)); expect_equal(sort(r$simYears$run), c(1001, 1002, 2001, 2002))
  expect_equal(r$simYears$area50[r$simYears$run == 1002], 500 + 100)
  expect_equal(r$simYears$area50[r$simYears$run == 1001], 300)
  expect_equal(nrow(r$simSizes), 6)
})

test_that("readHindcastObserved: scored years come from the hindcast's own files; escapes are the escape set's, not capped by ignitions", {
  root <- withr::local_tempdir(); mkSet(root, "9.9", 1986, 1)
  rec <- list(counts = data.frame(year = 1986:1987, ignitions = c(99, 7), escapes = c(9, 9), areaBurnedHa = 0, nFires50 = c(9, 2)),
              sizes = data.frame(year = c(1986, 1987, 1987, 1987, 1988), cause = c("natural", "natural", "natural", "allCause", NA), ha = c(5000, 80, 10, 400, NA)),
              escapeSizeHa = toyMinHa, fireCauses = c("natural", "unknown"))
  o <- readHindcastObserved(rec, root, "9.9", 1986L, minHa = toyMinHa)
  expect_equal(o$obsYears$ignitions, c(5, 7)); expect_equal(o$obsYears$area50, c(100, 80))
  expect_equal(o$obsYears$escapes, c(3, 9))                  # the files' and the record's escapes; 1987 has more escapes than ignitions, kept
  expect_equal(o$nbacAreaHa, 400); expect_equal(o$fitAreaHa, 400)
  expect_equal(o$escapeSizeHa, toyMinHa); expect_equal(o$fireCauses, c("natural", "unknown"))   # the record's definition travels with it
})

test_that("minHa is honoured: the simulated and observed escapes, area and sizes follow it, not the escape size of the record", {
  root <- withr::local_tempdir()
  for (s in 1:2) mkSet(root, "9.9", 1986, s)
  r <- readHindcastSim(root, "9.9", 1986L, sets = 1:2, minHa = 200)       # fires of 300, 100, 500 ha
  expect_equal(as.vector(r$simYears$escapes[order(r$simYears$run)]), c(1, 1, 1, 1))
  expect_equal(r$simYears$area50[r$simYears$run == 1002], 500); expect_equal(r$simYears$area50[r$simYears$run == 1001], 300)
  rec <- list(counts = data.frame(year = 1986:1987, ignitions = c(99, 7), escapes = c(9, 9)),
              sizes = data.frame(year = c(1987, 1987, 1987), cause = c("natural", "natural", "allCause"), ha = c(80, 300, 400)),
              escapeSizeHa = toyMinHa, fireCauses = "natural")
  o <- readHindcastObserved(rec, root, "9.9", 1986L, minHa = 200)
  expect_equal(o$obsYears$escapes, c(0, 1)); expect_equal(o$obsYears$area50, c(0, 300)); expect_equal(o$obsSizes$ha, 300)
  expect_equal(o$fitAreaHa, 400)
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
  pts <- terra::vect(cbind(500, 500), crs = "EPSG:3978"); pts$YEAR <- 2000L; pts$CAUSE <- "N"; pts$SIZE_HA <- 10
  poly <- terra::vect("POLYGON ((0 0, 2000 0, 2000 1000, 0 1000, 0 0))", crs = "EPSG:3978")
  poly$YEAR <- 2000L; poly$FIRECAUS <- "Natural"                                                   # 200 ha, 100 ha in column 1
  testthat::local_mocked_bindings(load_nfdb_points = function(...) pts, load_nbac_polys = function(...) poly,
                                  .package = "fireregimetools")
  cacheFile <- file.path(withr::local_tempdir(), "9.9.rds")
  saveRDS("a cache of the whole-polygon count", cacheFile)
  expect_no_warning(rec <- observedRecord(core, core, "nfdb", "nbac", years = 2000L, cacheFile = cacheFile))
  expect_equal(rec$counts$areaBurnedHaAllCause, 100, tolerance = 0.03)
  expect_true(file.exists(sub("\\.rds$", "_coreStart_natural-unknown_50ha.rds", cacheFile)))
  expect_identical(readRDS(cacheFile), "a cache of the whole-polygon count")                       # the old file is left alone
})

test_that("readHindcastSim: a set-year with no simulated fires (header-only sizes file) reads as zero area", {
  root <- withr::local_tempdir()
  for (s in 1:2) mkSet(root, "9.9", 1986, s)
  write.csv(data.frame(rep = integer(0), fire_id = integer(0), ha = numeric(0)), file.path(root, "9.9", 1986, 2, "hindcast_sizes.csv"), row.names = FALSE)
  r <- readHindcastSim(root, "9.9", 1986L, sets = 1:2, minHa = toyMinHa)
  expect_equal(nrow(r$simSizes), 3)
  expect_equal(as.vector(r$simYears$area50[r$simYears$run %in% c(2001, 2002)]), c(0, 0))
  expect_equal(r$simYears$area50[r$simYears$run == 1002], 600)
})

test_that("observedRecord: the default fireCauses is natural + unknown, as in the fit; natural alone leaves the unknown-cause fire out, in its own cache", {
  skip_if_not("fireCauseCodes" %in% getNamespaceExports("fireSenseUtils"))
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = "EPSG:3978")
  core <- terra::setValues(r, ifelse(terra::colFromCell(r, seq_len(terra::ncell(r))) <= 5, 1, NA))
  pts <- terra::vect(cbind(c(500, 2500), 500), crs = "EPSG:3978"); pts$YEAR <- 2000L; pts$CAUSE <- c("N", "U"); pts$SIZE_HA <- 10
  poly <- terra::vect(c("POLYGON ((0 0, 1000 0, 1000 1000, 0 1000, 0 0))",
                        "POLYGON ((2000 0, 3000 0, 3000 1000, 2000 1000, 2000 0))"), crs = "EPSG:3978")
  poly$YEAR <- 2000L; poly$FIRECAUS <- c("Natural", "Undetermined")
  testthat::local_mocked_bindings(load_nfdb_points = function(...) pts, load_nbac_polys = function(...) poly,
                                  .package = "fireregimetools")
  cacheFile <- file.path(withr::local_tempdir(), "9.9.rds")
  dflt <- observedRecord(core, core, "nfdb", "nbac", years = 2000L, cacheFile = cacheFile)
  expect_equal(dflt$counts$ignitions, 2)
  expect_equal(dflt$counts$areaBurnedHa, 200, tolerance = 0.03)
  nat <- observedRecord(core, core, "nfdb", "nbac", years = 2000L, cacheFile = cacheFile, fireCauses = "natural")
  expect_equal(nat$counts$ignitions, 1)
  expect_equal(nat$counts$areaBurnedHa, 100, tolerance = 0.03)
  expect_true(file.exists(sub("\\.rds$", "_coreStart_natural-unknown_50ha.rds", cacheFile)))
  expect_true(file.exists(sub("\\.rds$", "_coreStart_natural_50ha.rds", cacheFile)))
  expect_identical(formals(readFireSenseHindcast)$fireCauses, quote(fireSenseUtils::fireSenseFireCauses))
  expect_identical(formals(observedRecord)$fireCauses, quote(fireSenseUtils::fireSenseFireCauses))
  expect_identical(formals(observedRecord)$escapeSizeHa, quote(fireSenseUtils::fireSenseEscapeSizeHa))
})

test_that("hindcastCore finds the core raster with or without the NRV folders", {
  g <- terra::rast(nrows = 2, ncols = 2, vals = 1:4)
  for (sub in c("", "NRV/370/1")) {
    root <- withr::local_tempdir()
    dir.create(file.path(root, "hc", "9.9", "1986", "1"), recursive = TRUE)
    terra::writeRaster(g, file.path(root, "hc", "9.9", "1986", "1", "hindcast_burnMap.tif"))
    dir.create(file.path(root, "in", "9.9", sub), recursive = TRUE)
    terra::writeRaster(terra::rast(nrows = 2, ncols = 2, vals = c(2, 2, 1, 2)), file.path(root, "in", "9.9", sub, "rtml_9.9.tif"))
    core <- hindcastCore(file.path(root, "hc"), "9.9", file.path(root, "in"))
    expect_equal(as.vector(terra::values(core)), c(1, 1, NA, 1))
  }
})

test_that("burnRate is the mean annual observed area over the burnable area, in % per year", {
  core <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = "EPSG:3978")   # 1000 m cells = 100 ha
  core <- terra::setValues(core, ifelse(seq_len(100) <= 40, 1, NA))                                                 # 40 burnable cells
  b <- burnRate(core, data.frame(year = 2000:2003, area50 = c(0, 4, 0, 4)))
  expect_equal(b$burnableHa, 40 * 1000^2 / 1e4)                           # cells are 1000 m: 100 ha each
  expect_equal(b$burnRate, 100 * 2 / b$burnableHa)
})

test_that("observedRecord counts the fires that start in the core at their whole extent, buffer-start fires not at all, and records its escape definition", {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = "EPSG:3978")   # 100 ha cells
  col <- terra::colFromCell(r, seq_len(terra::ncell(r)))
  core <- terra::setValues(r, ifelse(col <= 5, 1, NA))                                                        # core: columns 1-5
  flammable <- terra::setValues(r, 1)                                                                          # burnable: core and buffer
  pts <- terra::vect(cbind(c(3000, 8500), 500), crs = "EPSG:3978"); pts$YEAR <- 2000L; pts$CAUSE <- "N"; pts$SIZE_HA <- 10
  poly <- terra::vect(c("POLYGON ((2000 0, 7000 0, 7000 1000, 2000 1000, 2000 0))",                          # 500 ha, starts in the core, 200 ha in the buffer
                        "POLYGON ((7000 0, 10000 0, 10000 1000, 7000 1000, 7000 0))"), crs = "EPSG:3978")     # 300 ha, starts in the buffer
  poly$YEAR <- 2000L; poly$FIRECAUS <- "Natural"
  testthat::local_mocked_bindings(load_nfdb_points = function(...) pts, load_nbac_polys = function(...) poly, .package = "fireregimetools")
  rec <- observedRecord(core, flammable, "nfdb", "nbac", years = 2000L)
  expect_equal(rec$counts$areaBurnedHa, 500, tolerance = 0.03)              # whole extent, core and buffer; the buffer-start fire is not counted
  expect_equal(rec$counts$ignitions, 1); expect_equal(rec$counts$escapes, 1)
  expect_equal(rec$sizes$ha[rec$sizes$cause == "natural"], 500, tolerance = 0.03)
  expect_equal(rec$escapeSizeHa, 50); expect_equal(rec$fireCauses, c("natural", "unknown"))
  big <- observedRecord(core, flammable, "nfdb", "nbac", years = 2000L, escapeSizeHa = 600)   # escapeSizeHa is the parameter, not a constant
  expect_equal(big$counts$escapes, 0); expect_equal(big$escapeSizeHa, 600)
})
