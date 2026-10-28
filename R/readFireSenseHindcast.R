## fireSense-specific reader. Everything else in R/ takes plain tables; this file turns FireSenseTesting hindcast outputs
## (/mnt/fast/outputs/hindcast<ledgerSuffix>/<ELF>/<year>/<set>/), the fit ledger and NFDB/NBAC observations into the ELF input
## list that runWave1() takes. Observed fire uses fireSenseUtils::hindcastObserved() (the fit's own definitions).

fread0 <- function(f) if (file.exists(f)) data.table::fread(f)

#' Which scored years of an ELF have every set finished?
#'
#' A set is finished when `hindcast_reps.csv`, `hindcast_sizes.csv`, `hindcast_observed.csv` and `hindcast_observedSizes.csv` exist.
#' @param hindcastDir,elf,years,sets Where to look. @return Integer years with all `sets` finished.
hindcastFinishedYears <- function(hindcastDir, elf, years = seq(1986L, 2021L, 5L), sets = 1:5) {
  need <- c("hindcast_reps.csv", "hindcast_sizes.csv", "hindcast_observed.csv", "hindcast_observedSizes.csv")
  years[vapply(years, function(y) all(vapply(sets, function(s) all(file.exists(file.path(hindcastDir, elf, y, s, need))), TRUE)), TRUE)]
}

#' Simulated side: one row per run-year, and the simulated fire sizes
#'
#' Reads `hindcast_reps.csv` and `hindcast_sizes.csv` of every set. A run is a rep of a set (`run = set * 1000 + rep`).
#' `area50` is the area of the run's fires of `minHa` or more. An escape is a fire of `minHa` or more, so `escapes` is
#' the hindcast's `nFires50` (not its `escapes` column, which is the fit's ignition-point definition and differs from the
#' count of fires by size).
#' @param hindcastDir,elf,years,sets As [hindcastFinishedYears()]. @param minHa Escape size.
#' @return `list(simYears, simSizes)`.
readHindcastSim <- function(hindcastDir, elf, years, sets = 1:5, minHa = 50) {
  both <- lapply(years, function(y) lapply(sets, function(s) {
    d <- file.path(hindcastDir, elf, y, s)
    r <- fread0(file.path(d, "hindcast_reps.csv")); z <- fread0(file.path(d, "hindcast_sizes.csv"))
    run <- function(x) s * 1000L + x$rep
    list(reps = data.frame(year = y, run = run(r), ignitions = r$ignitions, escapes = r$nFires50),
         sizes = data.frame(year = rep(y, NROW(z)), run = if (NROW(z)) run(z) else integer(0), ha = if (NROW(z)) z$ha else numeric(0)))
  }))
  flat <- unlist(both, recursive = FALSE)
  simYears <- do.call(rbind, lapply(flat, `[[`, "reps")); simSizes <- do.call(rbind, lapply(flat, `[[`, "sizes"))
  big <- simSizes[simSizes$ha >= minHa, ]
  key <- paste(simYears$year, simYears$run); a <- tapply(big$ha, paste(big$year, big$run), sum)
  simYears$area50 <- ifelse(is.na(a[key]), 0, a[key])
  list(simYears = simYears, simSizes = simSizes)
}

#' The ELF core on the hindcast grid (1 where the simulation can burn, NA elsewhere)
#' @param hindcastDir,elf As above. @param inputsDir Where the ELF inputs are (`<inputsDir>/<ELF>/NRV/370/1/rtml_<ELF>.tif`).
#' @return `SpatRaster`.
hindcastCore <- function(hindcastDir, elf, inputsDir = "/mnt/fast/inputs") {
  f <- Sys.glob(file.path(hindcastDir, elf, "*", "*", "hindcast_burnMap.tif"))[1]
  g <- terra::rast(f)
  core <- terra::crop(terra::rast(file.path(inputsDir, elf, "NRV/370/1", paste0("rtml_", elf, ".tif"))), g)
  terra::ifel(core == 2 & !is.na(g), 1, NA)
}

#' Observed fire of every year of the record, with the fit's definitions
#'
#' [fireSenseUtils::hindcastObserved()] for each year on NFDB points and NBAC polygons clipped to the core. Cached in `cacheFile`.
#' @param core [hindcastCore()]. @param nfdbPath,nbacPath The NFDB point shapefile and the NBAC polygon shapefile.
#' @param years The record. @param cacheFile An `.rds` to read if present and write otherwise; `_burnable` is added to its name
#'   (observed area is counted on the burnable cells of `core`, `flammableRTM = core`), so caches of the earlier whole-polygon count are not read.
#' @param cores Forked workers.
#' @return `list(counts, sizes, polys)`: counts per year, sizes with `cause` and `year`, polygons (wrapped).
observedRecord <- function(core, nfdbPath, nbacPath, years = 1986:2024, cacheFile = NULL, cores = 1L) {
  if (!is.null(cacheFile)) cacheFile <- sub("(\\.rds)?$", "_burnable\\1", cacheFile)
  if (!is.null(cacheFile) && file.exists(cacheFile)) return(readRDS(cacheFile))
  sa <- terra::as.polygons(terra::ext(core), crs = terra::crs(core))
  pts <- fireregimetools::load_nfdb_points(nfdbPath, sa, years, min_size_ha = 0)
  pol <- fireregimetools::load_nbac_polys(nbacPath, sa, years)
  o <- parallel::mclapply(years, function(y) { x <- fireSenseUtils::hindcastObserved(core, pts, pol, y, flammableRTM = core)
    x$polys <- if (!is.null(x$polys) && NROW(x$polys)) terra::wrap(x$polys) else NULL   # SpatVectors cannot leave a fork unwrapped
    x }, mc.cores = cores)
  r <- list(counts = data.table::rbindlist(lapply(o, `[[`, "counts")),
            sizes = data.table::rbindlist(Map(function(x, y) cbind(year = y, x$sizes), o, years)),
            polys = stats::setNames(lapply(o, `[[`, "polys"), years))
  if (!is.null(cacheFile)) saveRDS(r, cacheFile)
  r
}

#' Observed side in the ELF input format
#'
#' For the scored years the hindcast's own observed files are used (their core is year-specific); the other years come from
#' [observedRecord()].
#' @param rec [observedRecord()]. @param hindcastDir,elf,scored Where the hindcast's own observed files are, and for which years. @param minHa Escape size.
#' @return `list(obsYears, obsSizes, nbacAreaHa, fitAreaHa)`; the last two are the area (ha) of all-cause observed fires and
#'   of those of `minHa` or more (the fires the fit used), over the record.
readHindcastObserved <- function(rec, hindcastDir, elf, scored, minHa = 50) {
  own <- function(y) list(counts = fread0(file.path(hindcastDir, elf, y, 1, "hindcast_observed.csv")),
                          sizes = fread0(file.path(hindcastDir, elf, y, 1, "hindcast_observedSizes.csv")))
  oc <- lapply(scored, own)
  counts <- rbind(as.data.frame(rec$counts)[!rec$counts$year %in% scored, c("year", "ignitions", "nFires50")],
                  do.call(rbind, lapply(oc, function(x) as.data.frame(x$counts)[, c("year", "ignitions", "nFires50")])))
  names(counts)[names(counts) == "nFires50"] <- "escapes"             # an escape is a fire of minHa or more
  ## escapes come from burned-area polygons, ignitions from ignition points; a year can have more of the
  ## former (6.3.2: 5 years), which no escape-rate model can fit, so cap escapes at ignitions
  over <- which(counts$escapes > counts$ignitions)
  if (length(over)) {
    message(elf, ": ", length(over), " record year(s) with more escapes (fires >= ", minHa,
            " ha) than ignitions; escapes capped at ignitions for the escape rate")
    counts$escapes[over] <- counts$ignitions[over]
  }
  sizes <- rbind(as.data.frame(rec$sizes)[!rec$sizes$year %in% scored, c("year", "cause", "ha")],
                 do.call(rbind, Map(function(x, y) data.frame(year = rep(y, NROW(x$sizes)), cause = x$sizes$cause, ha = x$sizes$ha), oc, scored)))
  sizes <- sizes[!is.na(sizes$ha), ]                      # a year without fire is one NA row in the record
  nat <- sizes[sizes$cause == "natural" & sizes$ha >= minHa, ]
  counts <- counts[order(counts$year), ]
  counts$area50 <- as.numeric(tapply(nat$ha, factor(nat$year, counts$year), sum)); counts$area50[is.na(counts$area50)] <- 0
  all <- sizes[sizes$cause == "allCause", ]
  list(obsYears = counts, obsSizes = nat[, c("year", "ha")], nbacAreaHa = sum(all$ha), fitAreaHa = sum(all$ha[all$ha >= minHa]))
}

#' ELF median of the fit's climate covariate (CMD) for every year of the record
#'
#' Reads `climate/CMD_historical_<ELF>.tif` (one layer per year) from the first set that has it and averages it over the core.
#' @param hindcastDir,elf As above. @param core [hindcastCore()]. @param variable Climate variable of the fit.
#' @return `data.frame(year, value)`.
readHindcastCovariate <- function(hindcastDir, elf, core, variable = "CMD") {
  f <- Sys.glob(file.path(hindcastDir, elf, "*", "*", "climate", paste0(variable, "_historical_", elf, ".tif")))[1]
  r <- terra::rast(f)
  cr <- terra::crop(r, terra::ext(core))
  m <- if (terra::compareGeom(cr, core, stopOnError = FALSE)) core else terra::resample(core, cr, "near")
  v <- apply(terra::values(terra::mask(cr, m)), 2, stats::median, na.rm = TRUE)
  data.frame(year = as.integer(sub(".*_(\\d{4})$", "\\1", names(r))), value = v)
}

#' The ledger row of an ELF: what the gates need
#'
#' @param ledgerFile A spread-fit ledger `.rds` (`fireSenseParams_*.rds`; the caller downloads it). @param elf ELF (`polygonID`).
#' @return List: `objective` (`objFunVal` of the kept sets), `pc` (`maxAsymptote` of each kept set) and `p50`
#'   (`NULL` when the ledger has no power-law columns).
readLedgerRow <- function(ledgerFile, elf) {
  sp <- readRDS(ledgerFile); i <- which(sp$polygonID == elf)[1]
  p <- sp$params[[i]]
  list(objective = as.numeric(unlist(sp$objFunVal[[i]])), pc = as.numeric(p[["maxAsymptote"]]),
       p50 = if ("powerLawPc50_spread" %in% names(sp)) sp$powerLawPc50_spread[[i]] else NULL)
}

#' One ELF's hindcast as the input of [runWave1()]
#'
#' @param elf ELF name. @param hindcastDir `/mnt/fast/outputs/hindcast<ledgerSuffix>`. @param ledgerFile See [readLedgerRow()].
#' @param nfdbPath,nbacPath Observation shapefiles. @param obsCacheDir Folder for the [observedRecord()] caches.
#' @param years Scored years (only those with every set finished are used). @param sets Sets. @param inputsDir See [hindcastCore()].
#' @param cores Forked workers for the observed record. @return The list [runWave1()] documents; `NULL` if no year is finished.
readFireSenseHindcast <- function(elf, hindcastDir, ledgerFile = NULL, nfdbPath, nbacPath, obsCacheDir,
                                  years = seq(1986L, 2021L, 5L), sets = 1:5, inputsDir = "/mnt/fast/inputs", cores = 1L) {
  done <- hindcastFinishedYears(hindcastDir, elf, years, sets)
  if (!length(done)) return(NULL)
  sim <- readHindcastSim(hindcastDir, elf, done, sets)
  core <- hindcastCore(hindcastDir, elf, inputsDir)
  rec <- observedRecord(core, nfdbPath, nbacPath, cacheFile = file.path(obsCacheDir, paste0(elf, ".rds")), cores = cores)
  obs <- readHindcastObserved(rec, hindcastDir, elf, done)
  fit <- if (is.null(ledgerFile)) list() else readLedgerRow(ledgerFile, elf)
  fit$fitAreaHa <- obs$fitAreaHa; fit$totalAreaHa <- obs$nbacAreaHa
  c(list(elf = elf), sim, obs[c("obsYears", "obsSizes")],
    list(covariate = readHindcastCovariate(hindcastDir, elf, core), fit = fit))
}
