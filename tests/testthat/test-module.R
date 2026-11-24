## The module run on a toy fireSense hindcast: .inputObjects reads it (readFireSenseHindcast()), init scores it and writes
## wave1.rds, wave1_results.csv and validation.html to outputPath, and its figures through Plots(). Toy ELF: see helper-source.R.

toyDefinition <- new.env()                      # the escape definition hindcastObserved() was last called with

## The hindcast folders, inputs and observed-record cache of toyElf(): set 1 of every scored year, a 10 x 10 core.
toyHindcastFiles <- function(e, root) {
  hc <- file.path(root, "hindcast"); inp <- file.path(root, "inputs", e$elf); obs <- file.path(root, "outputs", "obs")
  dir.create(inp, recursive = TRUE); dir.create(obs, recursive = TRUE, showWarnings = FALSE)   # a run on the unfinished hindcast already made obs
  g <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = "EPSG:3978", vals = 0)
  terra::writeRaster(terra::setValues(g, 2), file.path(inp, paste0("rtml_", e$elf, ".tif")))  # 2: the ELF core
  rec <- e$obsYears$year
  for (y in unique(e$simYears$year)) {
    d <- file.path(hc, e$elf, y, 1); dir.create(file.path(d, "climate"), recursive = TRUE)
    s <- e$simYears[e$simYears$year == y, ]; z <- e$simSizes[e$simSizes$year == y, ]
    write.csv(data.frame(rep = s$run, ignitions = s$ignitions, nFires50 = s$escapes), file.path(d, "hindcast_reps.csv"), row.names = FALSE)
    write.csv(data.frame(rep = z$run, ha = z$ha), file.path(d, "hindcast_sizes.csv"), row.names = FALSE)
    o <- e$obsYears[e$obsYears$year == y, ]; os <- e$obsSizes[e$obsSizes$year == y, ]
    write.csv(data.frame(year = y, ignitions = o$ignitions, escapes = o$escapes), file.path(d, "hindcast_observed.csv"), row.names = FALSE)
    write.csv(data.frame(cause = rep(c("natural", "allCause"), each = nrow(os)), ha = rep(os$ha, 2)),
              file.path(d, "hindcast_observedSizes.csv"), row.names = FALSE)
    terra::writeRaster(g, file.path(d, "hindcast_burnMap.tif"))
    terra::writeRaster(terra::rast(lapply(seq_along(rec), function(i) terra::setValues(g, e$covariate$value[i]))) |>
                         stats::setNames(paste0("CMD_", rec)), file.path(d, "climate", paste0("CMD_historical_", e$elf, ".tif")))
  }
  sp <- data.frame(polygonID = e$elf); sp$params <- list(data.frame(maxAsymptote = e$fit$pc)); sp$objFunVal <- list(e$fit$objective)
  saveRDS(sp, file.path(inp, "fireSenseParams_1986-2024_toy_suffix.rds"))
  ## the NFDB and NBAC shapefiles are not read: hindcastObserved() answers with the toy record, one year at a time, and notes the definition it was given
  testthat::local_mocked_bindings(load_nfdb_points = function(...) NULL, load_nbac_polys = function(...) NULL, .package = "fireregimetools", .env = parent.frame())
  testthat::local_mocked_bindings(hindcastObserved = function(core, pts, pol, year, escapeSizeHa, fireCauses, ...) {
    toyDefinition$escapeSizeHa <- escapeSizeHa; toyDefinition$fireCauses <- fireCauses
    os <- e$obsSizes[e$obsSizes$year == year, ]
    list(counts = data.table::data.table(year = year, ignitions = e$obsYears$ignitions[e$obsYears$year == year], escapes = e$obsYears$escapes[e$obsYears$year == year]),
         sizes = data.table::data.table(cause = rep(c("natural", "allCause"), each = nrow(os)), ha = rep(os$ha, 2)), polys = NULL)
  }, .package = "fireSenseUtils", .env = parent.frame())
  hc
}

## Run the module on toyElf()'s hindcast; `...` are further module parameters.
runToyModule <- function(root, e, hc, ...) {
  mp <- file.path(root, "modules"); dir.create(mp, showWarnings = FALSE)
  link <- file.path(mp, "fireSimulator_validation")
  if (!file.exists(link)) file.symlink(normalizePath(testthat::test_path("..", "..")), link) # the checkout's folder name may differ; a second link would be made inside the checkout
  withr::local_options(spades.useRequire = FALSE, spades.moduleCodeChecks = FALSE, .local_envir = parent.frame())
  SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = 1), modules = "fireSimulator_validation",
    params = list(fireSimulator_validation = utils::modifyList(   # `...` may replace inputsPath
      list(hindcastPath = hc, hindcastYears = sort(unique(e$simYears$year)),
           inputsPath = file.path(root, "inputs"), ledgerSuffix = "_suffix",
           nBoot = 20L, nNull = 19L, nPerm = 99L, .plots = "png", pageTitle = "Toy page"), list(...))),
    paths = list(modulePath = mp, outputPath = file.path(root, "outputs"), inputPath = file.path(root, "inputs"),
                 cachePath = file.path(root, "cache"), scratchPath = file.path(root, "scratch")))
}

test_that("the module reads a fireSense hindcast and writes wave1.rds, wave1_results.csv and validation.html to outputPath", {
  skip_if_not_installed("SpaDES.core")
  root <- withr::local_tempdir()
  e <- toyElf(seed = 3, name = "9.9")
  hc <- toyHindcastFiles(e, root)
  dir.create(file.path(hc, "8.8"))                                      # an ELF with no finished year is left out
  out <- file.path(root, "outputs")
  sim <- runToyModule(root, e, hc)
  expect_equal(vapply(sim$validationELFs, `[[`, "", "elf"), "9.9")
  expect_equal(sim$validationELFs[[1]]$fit$pc, e$fit$pc)                 # the ledger was found by its suffix
  expect_equal(nrow(sim$validationELFs[[1]]$obsYears), nrow(e$obsYears)) # the observed record came from <outputPath>/obs
  expect_equal(readRDS(file.path(out, "wave1.rds")), sim$validationResults)
  expect_equal(nrow(utils::read.csv(file.path(out, "wave1_results.csv"))), nrow(sim$validationResults$results))
  h <- paste(readLines(file.path(out, "validation.html"), warn = FALSE), collapse = "\n")
  expect_match(h, "<title>Toy page</title>", fixed = TRUE); expect_match(h, "\"panel_9.9\"", fixed = TRUE)
  ## Plots() took ".9" of "9.9" as an extension: the panels of 11.3 and 11.4 were both written to panel_11.png
  expect_true(file.exists(file.path(out, "figures", "fireSimulator_validation", "panel_9_9.png")))
})

test_that("escapeSizeHa, fireCauses and minHa reach the reader: the record of that definition is read, and minHa sets the scored size", {
  skip_if_not_installed("SpaDES.core")
  root <- withr::local_tempdir()
  e <- toyElf(seed = 3, name = "9.9")
  hc <- toyHindcastFiles(e, root)
  sim <- runToyModule(root, e, hc, escapeSizeHa = 100, fireCauses = "natural")
  expect_equal(toyDefinition$escapeSizeHa, 100); expect_equal(toyDefinition$fireCauses, "natural")   # the module's parameters reached the record
  x <- sim$validationELFs[[1]]
  expect_equal(x$escapeSizeHa, 100); expect_equal(x$fireCauses, "natural")            # properties of the observed record
  escapesOf <- function(z, ha) vapply(seq_len(nrow(z$simYears)), function(i)
    sum(z$simSizes$year == z$simYears$year[i] & z$simSizes$run == z$simYears$run[i] & z$simSizes$ha >= ha), 0)
  expect_equal(as.numeric(x$simYears$escapes), escapesOf(x, 100))                    # minHa = NA scores the escape size
  sim2 <- runToyModule(root, e, hc, escapeSizeHa = 100, fireCauses = "natural", minHa = 150)
  y <- sim2$validationELFs[[1]]
  expect_equal(as.numeric(y$simYears$escapes), escapesOf(y, 150))                    # minHa is honoured, not only drawn
  expect_lt(sum(y$simYears$escapes), sum(x$simYears$escapes))
  expect_lt(sum(y$simYears$area50), sum(x$simYears$area50))
  p <- SpaDES.core::params(sim2)$fireSimulator_validation
  expect_equal(p$escapeSizeHa, 100); expect_equal(p$fireCauses, "natural")
})

test_that("a cached run reads the hindcast files present when it runs: an empty folder gives a clear error, and the same parameters score the hindcast once it is finished", {
  skip_if_not_installed("SpaDES.core")
  root <- withr::local_tempdir()
  e <- toyElf(seed = 3, name = "9.9")
  hc <- file.path(root, "hindcast"); dir.create(hc)
  expect_error(runToyModule(root, e, hc, .useCache = TRUE), "no finished hindcasts")   # was: a NULL-names error in runWave1()
  expect_error(runWave1(list(), minHa = toyMinHa), "no ELFs to score")
  toyHindcastFiles(e, root)                                              # the hindcast finishes; same parameters, same cachePath
  sim <- runToyModule(root, e, hc, .useCache = TRUE)                     # was: the .inputObjects of the empty folder, from the cache
  expect_equal(vapply(sim$validationELFs, `[[`, "", "elf"), "9.9")
  expect_true(file.exists(file.path(root, "outputs", "wave1.rds")))
})

test_that("the module writes wave1_bySet.csv only when bySet = TRUE", {
  skip_if_not_installed("SpaDES.core")
  e <- toyElf(seed = 3, name = "9.9")
  root <- withr::local_tempdir(); hc <- toyHindcastFiles(e, root); f <- file.path(root, "outputs", "wave1_bySet.csv")
  runToyModule(root, e, hc)
  expect_false(file.exists(f))
  runToyModule(root, e, hc, bySet = TRUE)
  b <- utils::read.csv(f)
  expect_equal(nrow(b), 1L); expect_equal(as.character(b$ELF), "9.9")
})

## move every climate folder of `hc` to the same place under `ref`
moveClimate <- function(hc, ref) for (d in list.dirs(hc)[basename(list.dirs(hc)) == "climate"]) {
  to <- file.path(ref, sub(paste0(hc, "/"), "", d, fixed = TRUE)); dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  file.rename(d, to)
}

test_that("a replay hindcast, which saves no climate, is scored with covariatePath's climate; without it the run stops naming the file", {
  skip_if_not_installed("SpaDES.core")
  root <- withr::local_tempdir()
  e <- toyElf(seed = 3, name = "9.9")
  hc <- toyHindcastFiles(e, root)
  ref <- file.path(root, "reference"); moveClimate(hc, ref)              # the reference hindcast keeps the climate
  expect_error(runToyModule(root, e, hc), "no CMD_historical_9.9.tif in")
  sim <- runToyModule(root, e, hc, covariatePath = ref)
  expect_equal(sim$validationELFs[[1]]$covariate$value, e$covariate$value, tolerance = 1e-6)   # a float32 raster
})

test_that("where the hindcast is read from and how many cores read it are not in the init event's cache key: a copy of the hindcast read with other cores is a cache hit", {
  skip_if_not_installed("SpaDES.core")
  ## a cache built on one machine is then a hit on another, where hindcastPath and inputsPath differ; the scored content (validationELFs) is in the key
  root <- withr::local_tempdir()
  e <- toyElf(seed = 3, name = "9.9")
  hc <- toyHindcastFiles(e, root)
  out <- file.path(root, "outputs", "wave1.rds")
  runToyModule(root, e, hc, .useCache = "init")
  expect_true(file.exists(out))
  file.remove(out)                                                       # writeOutputs writes it again, on a cache hit too
  root2 <- file.path(root, "copy"); dir.create(root2)
  file.copy(file.path(root, c("hindcast", "inputs")), root2, recursive = TRUE)
  sim <- runToyModule(root, e, file.path(root2, "hindcast"), .useCache = "init", inputsPath = file.path(root2, "inputs"),
                      covariatePath = file.path(root2, "hindcast"), nCores = 2L)
  expect_true(file.exists(out))
  expect_equal(vapply(sim$validationELFs, `[[`, "", "elf"), "9.9")
  p <- SpaDES.core::params(sim)$fireSimulator_validation
  expect_identical(p$nCores, 2L); expect_identical(p$hindcastPath, file.path(root2, "hindcast"))   # a hit keeps the current values
})

test_that("a run whose init is a cache hit still writes the page, figures and csv into an empty outputPath, and its results equal the uncached run's", {
  skip_if_not_installed("SpaDES.core")
  root <- withr::local_tempdir()
  e <- toyElf(seed = 3, name = "9.9")
  hc <- toyHindcastFiles(e, root)
  out <- file.path(root, "outputs")
  plain <- runToyModule(root, e, hc)                                     # no cache
  runToyModule(root, e, hc, .useCache = TRUE, bySet = TRUE)              # fills the cache
  unlink(c(file.path(out, c("wave1.rds", "wave1_results.csv", "wave1_bySet.csv", "validation.html")), file.path(out, "figures")), recursive = TRUE)
  hit <- runToyModule(root, e, hc, .useCache = TRUE, bySet = TRUE)       # init is restored from the cache
  for (f in c("wave1.rds", "wave1_results.csv", "wave1_bySet.csv", "validation.html")) expect_true(file.exists(file.path(out, f)), info = f)
  expect_true(file.exists(file.path(out, "figures", "fireSimulator_validation", "panel_9_9.png")))
  expect_equal(hit$validationResults, plain$validationResults)
  expect_equal(readRDS(file.path(out, "wave1.rds")), plain$validationResults)
})

test_that("the module writes the bias rows to wave1_results.csv and the bias table to wave1.rds and validation.html", {
  skip_if_not_installed("SpaDES.core")
  root <- withr::local_tempdir()
  e <- toyElf(seed = 3, name = "9.9")
  hc <- toyHindcastFiles(e, root)
  out <- file.path(root, "outputs")
  sim <- runToyModule(root, e, hc)
  b <- utils::read.csv(file.path(out, "wave1_results.csv")); b <- b[b$test == "B", ]
  expect_setequal(b$measure, c("ratio", "log2ratio", "meanRatio")); expect_setequal(b$null, c("area", "top10", "ignitions", "escapes"))
  expect_equal(readRDS(file.path(out, "wave1.rds"))$bias, sim$validationResults$bias)
  expect_match(paste(readLines(file.path(out, "validation.html"), warn = FALSE), collapse = "\n"), "all: median of ELFs", fixed = TRUE)
})

## readHindcastELFs(): each ELF's read is cached; a failed read stops the job naming the ELF
elfReadArgs <- function(root, elfs, reader, ...) {
  utils::modifyList(list(elfs, file.path(root, "hindcast"), years = 1986L, sets = 1L, inputsDir = file.path(root, "inputs"), ledgerSuffix = "_suffix",
         obsCacheDir = file.path(root, "obs"), elfCacheDir = file.path(root, "elfs"), fireCauses = "natural", escapeSizeHa = 50, minHa = 50,
         reader = reader), list(...))
}
## a reader that counts its calls in a file (it runs in forks), fails for the ELFs in `failFor`
countingReader <- function(counter, failFor = character(), nullFor = character()) {
  force(counter)
  function(elf, hindcastDir, ledgerFile = NULL, nfdbPath, nbacPath, obsCacheDir, years, sets, inputsDir, cores = 1L, fireCauses, escapeSizeHa, minHa, covariateDir) {
    cat(elf, "\n", file = counter, append = TRUE, sep = "")
    if (elf %in% failFor) stop("boom in ", elf)
    if (elf %in% nullFor) return(NULL)
    list(elf = elf, value = readLines(file.path(hindcastDir, elf, "1986", "1", "hindcast_reps.csv"))[2])
  }
}
mkFakeSet <- function(root, elf, y = 1986, s = 1) {
  d <- file.path(root, "hindcast", elf, y, s); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  for (f in c("hindcast_reps.csv", "hindcast_sizes.csv", "hindcast_observed.csv", "hindcast_observedSizes.csv")) writeLines(c("a", "1"), file.path(d, f))
}

test_that("readHindcastELFs reads an ELF once; a changed hindcast file re-reads only that ELF", {
  skip_if_not_installed("reproducible")
  root <- withr::local_tempdir(); counter <- file.path(root, "n.txt")
  for (e in c("1.1", "2.2")) mkFakeSet(root, e)
  rd <- countingReader(counter)
  a <- do.call(readHindcastELFs, elfReadArgs(root, c("1.1", "2.2"), rd))
  expect_equal(sort(readLines(counter)), c("1.1", "2.2"))
  b <- do.call(readHindcastELFs, elfReadArgs(root, c("1.1", "2.2"), rd))
  expect_equal(length(readLines(counter)), 2L)                     # the second call read nothing
  expect_identical(a, b)
  writeLines(c("a", "99"), file.path(root, "hindcast", "2.2", "1986", "1", "hindcast_reps.csv"))
  d <- do.call(readHindcastELFs, elfReadArgs(root, c("1.1", "2.2"), rd))
  expect_equal(sort(readLines(counter)), c("1.1", "2.2", "2.2"))   # only 2.2 again
  expect_equal(d[[2]]$value, "99"); expect_equal(d[[1]], a[[1]])
})

test_that("readHindcastELFs stops naming the ELF whose read failed or returned nothing, never drops it", {
  skip_if_not_installed("reproducible")
  root <- withr::local_tempdir(); counter <- file.path(root, "n.txt")
  for (e in c("1.1", "2.2", "3.3")) mkFakeSet(root, e)
  expect_error(suppressWarnings(do.call(readHindcastELFs, elfReadArgs(root, c("1.1", "2.2", "3.3"), countingReader(counter, failFor = "2.2"), nCores = 2L))),   # mclapply warns of the failed fork
               "2\\.2: .*boom in 2\\.2")
  expect_error(do.call(readHindcastELFs, elfReadArgs(root, c("1.1", "3.3"), countingReader(counter, nullFor = "3.3"), elfCacheDir = file.path(root, "elfs2"))),
               "3\\.3: returned NULL")
  ok <- do.call(readHindcastELFs, elfReadArgs(root, c("1.1", "2.2", "3.3"), countingReader(counter)))   # a failed read was not cached
  expect_length(ok, 3L)
})

test_that("readHindcastELFs: the same hindcast under another root is a cache hit; a hit equals a miss, data.table elements included; an unscored year changes nothing", {
  skip_if_not_installed("reproducible")
  root <- withr::local_tempdir(); counter <- file.path(root, "n.txt")
  mkFakeSet(root, "1.1")
  rd <- countingReader(counter)
  dtReader <- function(...) { x <- rd(...); x$dt <- data.table::data.table(a = 1:3, b = letters[1:3]); x }
  a <- do.call(readHindcastELFs, elfReadArgs(root, "1.1", dtReader))
  root2 <- withr::local_tempdir(); file.copy(file.path(root, "hindcast"), root2, recursive = TRUE, copy.date = TRUE)   # same files, same times, another place
  args2 <- elfReadArgs(root2, "1.1", dtReader, elfCacheDir = file.path(root, "elfs"))
  b <- do.call(readHindcastELFs, args2)
  expect_equal(length(readLines(counter)), 1L)                     # root2 is a hit
  expect_true(isTRUE(all.equal(a, b)))
  expect_s3_class(b[[1]]$dt, "data.table"); expect_equal(b[[1]]$dt$b, letters[1:3])
  mkFakeSet(root, "1.1", y = 2001)                                 # a year that is not scored
  do.call(readHindcastELFs, elfReadArgs(root, "1.1", dtReader))
  expect_equal(length(readLines(counter)), 1L)
})
