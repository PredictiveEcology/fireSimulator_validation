## The module run on a toy fireSense hindcast: .inputObjects reads it (readFireSenseHindcast()), init scores it and writes
## wave1.rds, wave1_results.csv and validation.html to outputPath, and its figures through Plots(). Toy ELF: see helper-source.R.

toyDefinition <- new.env()                      # the escape definition hindcastObserved() was last called with

## The hindcast folders, inputs and observed-record cache of toyElf(): set 1 of every scored year, a 10 x 10 core.
toyHindcastFiles <- function(e, root) {
  hc <- file.path(root, "hindcast"); inp <- file.path(root, "inputs", e$elf); obs <- file.path(root, "outputs", "obs")
  dir.create(inp, recursive = TRUE); dir.create(obs, recursive = TRUE)
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
    params = list(fireSimulator_validation = list(hindcastPath = hc, hindcastYears = sort(unique(e$simYears$year)),
                                                  inputsPath = file.path(root, "inputs"), ledgerSuffix = "_suffix",
                                                  nBoot = 20L, nNull = 19L, nPerm = 99L, .plots = "png", pageTitle = "Toy page", ...)),
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
