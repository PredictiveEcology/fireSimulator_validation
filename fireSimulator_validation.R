## Everything in this file and any files in the R directory are sourced during `simInit()`.
## The scoring functions in R/ take plain tables only; the module file hands its input to the runner and saves the results.
## With hindcastPath set and no validationELFs supplied, .inputObjects reads fireSense hindcasts (R/readFireSenseHindcast.R).
defineModule(sim, list(
  name = "fireSimulator_validation",
  description = paste("Scores a fire simulator against observed fire: gates and skill tests on annual ignitions, escape rate,",
                      "area burned, fire sizes and ranking of years, against nulls fitted to the observations (wave 1 of the hindcast test plan)."),
  keywords = c("fire", "validation", "hindcast", "CRPS"),
  authors = c(person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(fireSimulator_validation = "0.0.0.9013"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  documentation = list("README.md"),
  reqdPkgs = list("data.table", "MASS", "withr", "ggplot2", "patchwork", "scales", "jsonlite", "terra",
                  "FOR-CAST/fireregimetools@main (>= 0.1.0.9001)", "PredictiveEcology/fireSenseUtils@development (>= 0.2.3.9139)"),
  parameters = bindrows(
    defineParameter("nBoot", "integer", 999L, 1L, NA, "Bootstrap resamples of ELFs and years for the pooled verdicts."),
    defineParameter("nNull", "integer", 199L, 1L, NA, "Draws per year of the Anderson-Darling null."),
    defineParameter("nPerm", "integer", 9999L, 1L, NA, "Permutations of the ranking floor."),
    defineParameter("seed", "integer", 1L, NA, NA, "Seed of the bootstrap and permutations."),
    defineParameter("simLabel", "character", "simulated", NA, NA, "Name of the simulated series in the figures."),
    defineParameter("escapeSizeHa", "numeric", fireSenseUtils::fireSenseEscapeSizeHa, 0, NA,
                    paste("Smallest whole fire polygon (ha) that is an escape: the escape definition of the fit and the hindcast, recorded as a property of",
                          "the observed record. Checked against `escapeSizeHa` of the other modules.")),
    defineParameter("fireCauses", "character", fireSenseUtils::fireSenseFireCauses, NA, NA,
                    paste("Cause groups (`fireSenseUtils::fireCauseCodes()`) whose fires are counted as observed fire: the fit's and the hindcast's `fireCauses`,",
                          "recorded as a property of the observed record. Checked against `fireCauses` of the other modules.")),
    defineParameter("minHa", "numeric", NA_real_, 0, NA,
                    "Smallest fire (ha) scored and drawn: escapes, area and sizes of the simulated and observed series. NA: `escapeSizeHa`."),
    defineParameter(".plots", "character", "screen", NA, NA, "Where the figures go, as for SpaDES.core::Plots(): 'screen', 'png', both, or NA for none."),
    defineParameter(".plotInitialTime", "numeric", 1, NA, NA, "When the figures are drawn on screen; NA for never."),
    ## the fireSense reader: used only when validationELFs is not supplied
    defineParameter("hindcastPath", "character", NA_character_, NA, NA,
                    "FireSense hindcast outputs, <hindcastPath>/<ELF>/<year>/<set>/. When set and validationELFs is not supplied, every ELF there that has finished all hindcastYears is read with readFireSenseHindcast()."),
    defineParameter("hindcastYears", "integer", 1986:2024, NA, NA, "Scored years; an ELF is scored only when every one is finished."),
    defineParameter("hindcastSets", "integer", 1L, NA, NA, "Hindcast parameter sets (<set> folders) read."),
    defineParameter("skipELFs", "character", character(0), NA, NA, "ELFs of hindcastPath left out."),
    defineParameter("inputsPath", "character", NA_character_, NA, NA,
                    "Folder of the ELF input folders: <inputsPath>/<ELF>/ has rtml_<ELF>.tif, the spread-fit ledger fireSenseParams_*<ledgerSuffix>.rds and the NFDB_point_*.shp and NBAC_merged_*.shp observations (the latest of each)."),
    defineParameter("ledgerSuffix", "character", "", NA, NA, "Suffix of the spread-fit ledger, e.g. '_crpsCauseAnnD03WholeLong_20261108'."),
    defineParameter("nCores", "integer", 1L, 1L, NA, "ELFs read at once (forked)."),
    ## the page: <outputPath>/validation.html, beside wave1.rds and wave1_results.csv
    defineParameter("pageTitle", "character", "Hindcast validation", NA, NA, "Title of validation.html."),
    defineParameter("pageSubtitle", "character", "", NA, NA, "Subtitle of validation.html."),
    defineParameter("pageDpi", "numeric", 50, 1, NA, "Resolution of the page's figures.")
  ),
  inputObjects = bindrows(
    expectsInput("validationELFs", "list",
                 "One plain list per ELF: elf, simYears, simSizes, obsYears, obsSizes, covariate, fit (see runWave1()). A simulator-specific reader builds them; for fireSense, readFireSenseHindcast().")
  ),
  outputObjects = bindrows(
    createsOutput("validationResults", "list", "runWave1() result: results, gates, perELF, flagged. Also saved to <outputPath>/wave1.rds, its results table to wave1_results.csv, and the page to validation.html.")
  )
))

## the size scored: minHa, or the escape size when it is NA
scoredMinHa <- function(p) if (is.na(p$minHa)) p$escapeSizeHa else p$minHa

.inputObjects <- function(sim) {
  fireSenseUtils::escapeDefinitionOf(sim)     # stops if another module has another escapeSizeHa or fireCauses
  ## fireSense hindcasts: every ELF that has finished all its years; the observed records are cached in <outputPath>/obs
  if (!suppliedElsewhere("validationELFs", sim) && !is.na(P(sim)$hindcastPath)) {
    hc <- P(sim)$hindcastPath; yrs <- P(sim)$hindcastYears; sets <- P(sim)$hindcastSets; inp <- P(sim)$inputsPath
    obsDir <- file.path(outputPath(sim), "obs"); dir.create(obsDir, recursive = TRUE, showWarnings = FALSE)
    done <- Filter(function(e) length(hindcastFinishedYears(hc, e, yrs, sets)) == length(yrs),
                   setdiff(list.files(hc), P(sim)$skipELFs))
    elfs <- parallel::mclapply(done, mc.cores = P(sim)$nCores, mc.preschedule = FALSE, function(e)
      readFireSenseHindcast(e, hc, years = yrs, sets = sets, inputsDir = inp,
                            ledgerFile = Sys.glob(file.path(inp, e, paste0("fireSenseParams_*", P(sim)$ledgerSuffix, ".rds")))[1],
                            nfdbPath = tail(Sys.glob(file.path(inp, e, "NFDB_point_*.shp")), 1),
                            nbacPath = tail(Sys.glob(file.path(inp, e, "NBAC_merged_*.shp")), 1),
                            obsCacheDir = obsDir, fireCauses = P(sim)$fireCauses,
                            escapeSizeHa = P(sim)$escapeSizeHa, minHa = scoredMinHa(P(sim))))
    sim$validationELFs <- Filter(Negate(is.null), elfs)
    message("scoring ", length(sim$validationELFs), " ELFs")
  }
  return(invisible(sim))
}

doEvent.fireSimulator_validation <- function(sim, eventTime, eventType) {
  switch(eventType,
    init = {
      sim$validationResults <- runWave1(sim$validationELFs, minHa = scoredMinHa(P(sim)), nBoot = P(sim)$nBoot, nNull = P(sim)$nNull,
                                        nPerm = P(sim)$nPerm, seed = P(sim)$seed)
      summ <- summariseGates(sim$validationResults$gates)
      pl <- stats::setNames(lapply(sim$validationELFs, elfPlots, label = P(sim)$simLabel, minHa = scoredMinHa(P(sim))),
                            vapply(sim$validationELFs, `[[`, "", "elf"))
      for (e in sim$validationELFs)
        Plots(data = e, fn = plotELFPanel, summary = summ, label = P(sim)$simLabel, minHa = scoredMinHa(P(sim)), plots = pl[[e$elf]],
              filename = paste0("panel_", gsub(".", "_", e$elf, fixed = TRUE)), ggsaveArgs = panelSize) # Plots() takes ".3" of "11.3" as an extension
      Plots(data = summ, fn = plotOverview, filename = "overview", ggsaveArgs = overviewSize)
      for (m in names(metrics))
        Plots(data = pl, fn = plotMetricPanel, summary = summ, metric = m,
              filename = paste0("metric_", m), ggsaveArgs = metricSize(length(pl)))
      out <- outputPath(sim); dir.create(out, recursive = TRUE, showWarnings = FALSE)
      saveRDS(sim$validationResults, file.path(out, "wave1.rds"))
      utils::write.csv(sim$validationResults$results, file.path(out, "wave1_results.csv"), row.names = FALSE)
      buildValidationPage(sim$validationELFs, sim$validationResults$gates, file.path(out, "validation.html"),
                          template = file.path(modulePath(sim), currentModule(sim), "pageTemplate.html"),
                          title = P(sim)$pageTitle, subtitle = P(sim)$pageSubtitle, label = P(sim)$simLabel,
                          minHa = scoredMinHa(P(sim)), dpi = P(sim)$pageDpi)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}
