## fireSense-specific reader. Everything else in R/ takes plain tables; this file turns FireSenseTesting hindcast outputs
## (/mnt/fast/outputs/hindcast<ledgerSuffix>/<ELF>/<year>/<set>/), the fit ledger and NFDB/NBAC observations into the ELF input
## list that runWave1() takes. Observed fire uses fireSenseUtils::hindcastObserved(coreStart = TRUE), the fit's own escape
## definition (escapeSizeHa, fireCauses) and start rule: a fire belongs to the ELF when its start is in the core, at its whole extent.

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
#' `area50` is the area of the run's fires of `minHa` or more. An escape is a fire of `minHa` or more, so `escapes` counts the
#' fires of `hindcast_sizes.csv` of that size (not the hindcast's `escapes` column, which is the fit's ignition-point definition
#' and differs from the count of fires by size). A run-year with no such fire has 0 escapes.
#' @param hindcastDir,elf,years,sets As [hindcastFinishedYears()]. @param minHa Escape size.
#' @return `list(simYears, simSizes)`.
readHindcastSim <- function(hindcastDir, elf, years, sets = 1:5, minHa) {
  both <- lapply(years, function(y) lapply(sets, function(s) {
    d <- file.path(hindcastDir, elf, y, s)
    r <- fread0(file.path(d, "hindcast_reps.csv")); z <- fread0(file.path(d, "hindcast_sizes.csv"))
    run <- function(x) s * 1000L + x$rep
    list(reps = data.frame(year = y, run = run(r), ignitions = r$ignitions),
         sizes = data.frame(year = rep(y, NROW(z)), run = if (NROW(z)) run(z) else integer(0), ha = if (NROW(z)) z$ha else numeric(0)))
  }))
  flat <- unlist(both, recursive = FALSE)
  simYears <- do.call(rbind, lapply(flat, `[[`, "reps")); simSizes <- do.call(rbind, lapply(flat, `[[`, "sizes"))
  big <- simSizes[simSizes$ha >= minHa, ]
  key <- paste(simYears$year, simYears$run); a <- tapply(big$ha, paste(big$year, big$run), sum)
  simYears$area50 <- ifelse(is.na(a[key]), 0, a[key])
  n <- table(paste(big$year, big$run)); simYears$escapes <- ifelse(is.na(n[key]), 0L, as.integer(n[key]))
  list(simYears = simYears, simSizes = simSizes)
}

#' The hindcast burn map: NA where the simulation cannot burn, over the whole simulated landscape
#' @param hindcastDir,elf As above. @return `SpatRaster`, read once for [hindcastCore()] and [hindcastFlammable()].
hindcastBurnMap <- function(hindcastDir, elf) terra::rast(Sys.glob(file.path(hindcastDir, elf, "*", "*", "hindcast_burnMap.tif"))[1])

#' The ELF core on the hindcast grid (1 where the simulation can burn, NA elsewhere)
#'
#' The core is the cells of the ELF's `rtml` raster that [fireSenseUtils::inELFCore()] puts in the core (the rule of the fit and the burn).
#' @param hindcastDir,elf As above. @param inputsDir Where the ELF inputs are (`<inputsDir>/<ELF>/rtml_<ELF>.tif` on the GPSC, `<inputsDir>/<ELF>/NRV/370/1/rtml_<ELF>.tif` where an NRV run has been made).
#' @param burnMap [hindcastBurnMap()].
#' @return `SpatRaster`.
hindcastCore <- function(hindcastDir, elf, inputsDir = "/mnt/fast/inputs", burnMap = hindcastBurnMap(hindcastDir, elf)) {
  rtml <- file.path(inputsDir, elf, c("", "NRV/370/1"), paste0("rtml_", elf, ".tif"))
  zone <- terra::crop(terra::rast(rtml[file.exists(rtml)][1]), burnMap)
  terra::setValues(burnMap, ifelse(fireSenseUtils::inELFCore(seq_len(terra::ncell(zone)), zone) & !is.na(terra::values(burnMap)[, 1]), 1, NA))
}

#' The burnable cells of the hindcast grid, core and buffer (1 where the simulation can burn, NA elsewhere)
#'
#' @param burnMap [hindcastBurnMap()]. @return `SpatRaster`, the `flammableRTM` [observedRecord()] needs.
hindcastFlammable <- function(burnMap) terra::ifel(is.na(burnMap), NA, 1)

#' Observed fire of every year of the record, with the fit's definitions
#'
#' [fireSenseUtils::hindcastObserved()] (`coreStart = TRUE`) for each year on NFDB points and NBAC polygons: the fires that
#' start in the core, each at its whole extent (core and buffer); escapes are the polygons of `fireCauses` of `escapeSizeHa` or more
#' (whole polygon). A [reproducible::Cache()] in `cachePath`: a record made by another set of arguments, or by older escape code of
#' fireSenseUtils ([fireSenseUtils::escapeFireSetDeps()]), is not read.
#' @param core [hindcastCore()]. @param flammableRTM The burnable cells of the core AND the buffer ([hindcastFlammable()]); fires are
#'   started and measured on them.
#' @param nfdbPath,nbacPath The NFDB point shapefile and the NBAC polygon shapefile.
#' @param years The record. @param cachePath The cache folder (`NULL`: reproducible's default).
#' @param cores Forked workers (not part of the cache key).
#' @param fireCauses Cause groups (`fireSenseUtils::fireCauseCodes()`) whose fires are counted (labelled `"natural"` in `sizes`);
#'   default `fireSenseUtils::fireSenseFireCauses` (`c("natural", "unknown")`), the default of `fireSense_dataPrepFit` and `fireSense_summary`. Use the `fireCauses` the fit
#'   and the hindcast ran with.
#' @param escapeSizeHa Smallest whole polygon (ha) that is an escape; default `fireSenseUtils::fireSenseEscapeSizeHa`. Use the one the fit and
#'   the hindcast ran with.
#' @return `list(counts, sizes, polys, escapeSizeHa, fireCauses)`: counts per year (`escapes` is the number of escape polygons), sizes with
#'   `cause` and `year`, polygons (wrapped), and the escape definition the record was made with.
observedRecord <- function(core, flammableRTM, nfdbPath, nbacPath, years = 1986:2024, cachePath = NULL, cores = 1L,
                           fireCauses = fireSenseUtils::fireSenseFireCauses, escapeSizeHa = fireSenseUtils::fireSenseEscapeSizeHa) {
  reproducible::Cache(buildObservedRecord, core, flammableRTM, nfdbPath, nbacPath, years, cores, fireCauses, escapeSizeHa,
                      omitArgs = "cores", .cacheExtra = fireSenseUtils::escapeFireSetDeps(), cachePath = cachePath)
}

buildObservedRecord <- function(core, flammableRTM, nfdbPath, nbacPath, years, cores, fireCauses, escapeSizeHa) {
  sa <- terra::as.polygons(terra::ext(core), crs = terra::crs(core))
  pts <- fireregimetools::load_nfdb_points(nfdbPath, sa, years, min_size_ha = 0)
  pol <- fireregimetools::load_nbac_polys(nbacPath, sa, years)
  o <- parallel::mclapply(years, function(y) { x <- fireSenseUtils::hindcastObserved(core, pts, pol, y, escapeSizeHa = escapeSizeHa, flammableRTM = flammableRTM,
                                                  fireCauses = fireCauses, coreStart = TRUE)
    x$polys <- if (!is.null(x$polys) && NROW(x$polys)) terra::wrap(x$polys) else NULL   # SpatVectors cannot leave a fork unwrapped
    x }, mc.cores = cores)
  list(counts = data.table::rbindlist(lapply(o, `[[`, "counts")),
       sizes = data.table::rbindlist(Map(function(x, y) cbind(year = y, x$sizes), o, years)),
       polys = stats::setNames(lapply(o, `[[`, "polys"), years), escapeSizeHa = escapeSizeHa, fireCauses = fireCauses)
}

#' Observed side in the ELF input format
#'
#' For the scored years the hindcast's own observed files are used (their core is year-specific); the other years come from
#' [observedRecord()].
#' @param rec [observedRecord()]. @param hindcastDir,elf,scored Where the hindcast's own observed files are, and for which years.
#' @param minHa Escape size scored. When it is the record's `escapeSizeHa`, `escapes` are the record's escape polygons (the escape set,
#'   one per core-start escape, no cap by ignitions); otherwise the natural fires of `minHa` or more.
#' @return `list(obsYears, obsSizes, nbacAreaHa, fitAreaHa, escapeSizeHa, fireCauses)`; `nbacAreaHa` and `fitAreaHa` are the area (ha) of all-cause
#'   observed fires and of those of `minHa` or more (the fires the fit used), over the record; the last two are the record's escape definition.
readHindcastObserved <- function(rec, hindcastDir, elf, scored, minHa) {
  own <- function(y) list(counts = fread0(file.path(hindcastDir, elf, y, 1, "hindcast_observed.csv")),
                          sizes = fread0(file.path(hindcastDir, elf, y, 1, "hindcast_observedSizes.csv")))
  oc <- lapply(scored, own)
  cols <- c("year", "ignitions", "escapes")
  counts <- rbind(as.data.frame(rec$counts)[!rec$counts$year %in% scored, cols],
                  do.call(rbind, lapply(oc, function(x) as.data.frame(x$counts)[, cols])))
  sizes <- rbind(as.data.frame(rec$sizes)[!rec$sizes$year %in% scored, c("year", "cause", "ha")],
                 do.call(rbind, Map(function(x, y) data.frame(year = rep(y, NROW(x$sizes)), cause = x$sizes$cause, ha = x$sizes$ha), oc, scored)))
  sizes <- sizes[!is.na(sizes$ha), ]                      # a year without fire is one NA row in the record
  nat <- sizes[sizes$cause == "natural" & sizes$ha >= minHa, ]
  counts <- counts[order(counts$year), ]
  if (!isTRUE(minHa == rec$escapeSizeHa))                           # another size than the escape set's: count the fires
    counts$escapes <- as.integer(table(factor(nat$year, counts$year)))
  counts$area50 <- as.numeric(tapply(nat$ha, factor(nat$year, counts$year), sum)); counts$area50[is.na(counts$area50)] <- 0
  all <- sizes[sizes$cause == "allCause", ]
  list(obsYears = counts, obsSizes = nat[, c("year", "ha")], nbacAreaHa = sum(all$ha), fitAreaHa = sum(all$ha[all$ha >= minHa]),
       escapeSizeHa = rec$escapeSizeHa, fireCauses = rec$fireCauses)
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

#' Burnable area and burn rate of an ELF
#'
#' The burn rate is the mean annual observed area burned over the record divided by the burnable area, in % per year (the
#' conventional fire-cycle rate; its inverse is the fire cycle in years). The area is `area50` of `obsYears`, the area of
#' the natural fires of `minHa` or more, the quantity the area plots and the area gate use; smaller fires are not in the record's annual area.
#' @param core Burnable-cell raster ([hindcastCore()]; NA = not burnable). @param obsYears `data.frame(year, area50)`, every observed year.
#' @return `list(burnableHa, burnRate)`: area (ha) of the non-NA cells and the rate (% per year).
burnRate <- function(core, obsYears) {
  ha <- as.numeric(terra::global(!is.na(core), "sum")) * prod(terra::res(core)) / 1e4          # non-NA cells x cell area (the cores are projected, in m)
  list(burnableHa = ha, burnRate = 100 * mean(obsYears$area50) / ha)
}

#' One ELF's hindcast as the input of [runWave1()]
#'
#' @param elf ELF name. @param hindcastDir `/mnt/fast/outputs/hindcast<ledgerSuffix>`. @param ledgerFile See [readLedgerRow()].
#' @param nfdbPath,nbacPath Observation shapefiles. @param obsCacheDir Cache folder of [observedRecord()].
#' @param years Scored years (only those with every set finished are used). @param sets Sets. @param inputsDir See [hindcastCore()].
#' @param cores Forked workers for the observed record. @param fireCauses,escapeSizeHa See [observedRecord()]; use the ones the hindcast ran with.
#' @param minHa Escape size scored (the simulated and observed `escapes`, `area50` and sizes); default `escapeSizeHa`.
#' @return The list [runWave1()] documents, with `burnableHa` and `burnRate` ([burnRate()]) and the escape definition of the observed record
#'   (`escapeSizeHa`, `fireCauses`); `NULL` if no year is finished.
readFireSenseHindcast <- function(elf, hindcastDir, ledgerFile = NULL, nfdbPath, nbacPath, obsCacheDir,
                                  years = seq(1986L, 2021L, 5L), sets = 1:5, inputsDir = "/mnt/fast/inputs", cores = 1L,
                                  fireCauses = fireSenseUtils::fireSenseFireCauses, escapeSizeHa = fireSenseUtils::fireSenseEscapeSizeHa,
                                  minHa = escapeSizeHa) {
  done <- hindcastFinishedYears(hindcastDir, elf, years, sets)
  if (!length(done)) return(NULL)
  sim <- readHindcastSim(hindcastDir, elf, done, sets, minHa)
  burnMap <- hindcastBurnMap(hindcastDir, elf)
  core <- hindcastCore(hindcastDir, elf, inputsDir, burnMap)
  rec <- observedRecord(core, hindcastFlammable(burnMap), nfdbPath, nbacPath, cachePath = obsCacheDir,
                        cores = cores, fireCauses = fireCauses, escapeSizeHa = escapeSizeHa)
  obs <- readHindcastObserved(rec, hindcastDir, elf, done, minHa)
  fit <- if (is.null(ledgerFile)) list() else readLedgerRow(ledgerFile, elf)
  fit$fitAreaHa <- obs$fitAreaHa; fit$totalAreaHa <- obs$nbacAreaHa
  c(list(elf = elf), sim, obs[c("obsYears", "obsSizes", "escapeSizeHa", "fireCauses")], burnRate(core, obs$obsYears),
    list(covariate = readHindcastCovariate(hindcastDir, elf, core), fit = fit))
}
