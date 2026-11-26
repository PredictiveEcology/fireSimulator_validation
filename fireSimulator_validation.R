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
  version = list(fireSimulator_validation = "0.0.0.9025"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  documentation = list("README.md"),
  reqdPkgs = list("data.table", "reproducible", "MASS", "withr", "ggplot2", "patchwork", "scales", "jsonlite", "terra",
                  "PredictiveEcology/SpaDES.core@development (>= 3.2.1.9047)", # 9037: reads `.neverCache`; 9047: defineParameter(cacheKey = FALSE)
                  "FOR-CAST/fireregimetools@main (>= 0.1.0.9001)", "PredictiveEcology/fireSenseUtils@development (>= 0.2.3.9158)"),
  parameters = bindrows(
    defineParameter("nBoot", "integer", 999L, 1L, NA, "Bootstrap resamples of ELFs and years for the pooled verdicts."),
    defineParameter("nNull", "integer", 199L, 1L, NA, "Draws per year of the Anderson-Darling null."),
    defineParameter("nPerm", "integer", 9999L, 1L, NA, "Permutations of the ranking floor."),
    defineParameter("seed", "integer", 1L, NA, NA, "Seed of the bootstrap and permutations."),
    defineParameter("simLabel", "character", "simulated", NA, NA, "Name of the simulated series in the figures.",
                    cacheKey = FALSE), # used only by writeOutputs
    defineParameter("escapeSizeHa", "numeric", fireSenseUtils::fireSenseEscapeSizeHa, 0, NA,
                    paste("Smallest whole fire polygon (ha) that is an escape: the escape definition of the fit and the hindcast, recorded as a property of",
                          "the observed record. Checked against `escapeSizeHa` of the other modules.")),
    defineParameter("fireCauses", "character", fireSenseUtils::fireSenseFireCauses, NA, NA,
                    paste("Cause groups (`fireSenseUtils::fireCauseCodes()`) whose fires are counted as observed fire: the fit's and the hindcast's `fireCauses`,",
                          "recorded as a property of the observed record. Checked against `fireCauses` of the other modules.")),
    defineParameter("minHa", "numeric", NA_real_, 0, NA,
                    "Smallest fire (ha) scored and drawn: escapes, area and sizes of the simulated and observed series. NA: `escapeSizeHa`."),
    defineParameter(".neverCache", "character", c(".inputObjects", "writeOutputs"), NA, NA,
                    paste("Events that are never cached, whatever `.useCache` says (read by SpaDES.core, 3.2.1.9037 or later):",
                          "`.inputObjects` reads the hindcast files on disk, which the cache key does not hold, so a cached one returns the ELFs of an earlier, unfinished hindcast.",
                          "`writeOutputs` writes the figures, wave1 files and validation.html from `validationResults`, which a cached `init` restores without them.")),
    defineParameter(".plots", "character", "screen", NA, NA, "Where the figures go, as for SpaDES.core::Plots(): 'screen', 'png', both, or NA for none.",
                    cacheKey = FALSE), # used only by writeOutputs
    defineParameter(".plotInitialTime", "numeric", 1, NA, NA, "When the figures are drawn on screen; NA for never.",
                    cacheKey = FALSE), # used only by writeOutputs
    ## the fireSense reader: used only when validationELFs is not supplied
    defineParameter("hindcastPath", "character", NA_character_, NA, NA,
                    "FireSense hindcast outputs, <hindcastPath>/<ELF>/<year>/<set>/. When set and validationELFs is not supplied, every ELF there that has finished all hindcastYears is read with readFireSenseHindcast().",
                    cacheKey = FALSE), # only .inputObjects (never cached) reads it; the ELFs it gives are in init's key
    defineParameter("covariatePath", "character", NA_character_, NA, NA,
                    "Hindcast outputs the climate covariate (CMD) is read from, laid out as hindcastPath; NA: hindcastPath. A replay hindcast (fireSense_burn ignitionSource = 'replay') saves no climate: use its reference hindcast's folder, whose climate it replayed.",
                    cacheKey = FALSE), # as hindcastPath
    defineParameter("hindcastYears", "integer", 1986:2024, NA, NA, "Scored years; an ELF is scored only when every one is finished."),
    defineParameter("hindcastSets", "integer", 1L, NA, NA, "Hindcast parameter sets (<set> folders) read."),
    defineParameter("bySet", "logical", FALSE, NA, NA,
                    paste("Also score each parameter set of each ELF on its own (the `paramSet` column of hindcast_reps.csv and hindcast_sizes.csv, else the set folder):",
                          "gate ratios with single-run 2.5-97.5% ranges, annual-area skill and gate results, one row per set, written to wave1_bySet.csv.")),
    defineParameter("mode", "character", "score", NA, NA,
                    paste("'score': score the hindcast (the default). 'weights': tune objective weights from `weightsTerms`, `weightsBySet` and `weightsGrid` (`tuneWeights()`):",
                          "no hindcast is read; the tables `weights_report.csv` (the held-out ELFs' rows: the result), `weights_tuning.csv` (tune ELFs, labelled) and `weights_losses.csv` go to outputPath.",
                          "'anova': partition the variation of the per-ELF skill scores of several arms among ELF, the arms' factors and their interactions, from `anovaScores`, `anovaFactors` and `anovaMetrics` (`computeAnova()`):",
                          "no hindcast is read; `anova_terms.csv` and `anova_cells.csv` go to outputPath, with one figure.")),
    defineParameter("anovaScores", "character", NA_character_, NA, NA,
                    paste("Mode 'anova': score folders (outputPath of mode 'score': wave1_results.csv), one arm each. The arm is the name in this vector, else the folder's basename;",
                          "folders with one arm name are combined (the ELFs of one arm scored in several folders)."),
                    cacheKey = FALSE), # only .inputObjects (never cached) reads the files; what it read is in init's key
    defineParameter("anovaFactors", "ANY", NA_character_, NA, NA,
                    "Mode 'anova': csv path or data.frame, column `arm` (as named by `anovaScores`) and one column per factor.",
                    cacheKey = FALSE), # only .inputObjects (never cached) reads the files; what it read is in init's key
    defineParameter("anovaMetrics", "character", NA_character_, NA, NA,
                    "Mode 'anova': per-ELF skill metrics analysed, `<test>_<null>_skill` (e.g. S1_N1_skill, S3_N1_skill, S4largest_N1_skill). NA: every `*_N1_skill` in the scores."),
    defineParameter("weightsTerms", "character", NA_character_, NA, NA,
                    paste("Mode 'weights': csv files of raw objective terms, one per ELF as fireSense_spreadFit writes them (`spreadFitTerms_<runName>.csv`: `member`, `rep`, parameters, terms; the ELF is read from the file name (the one token like `4.3` or `6.6.2`; an ambiguous name must be given by naming the file with its ELF in this vector)),",
                          "or one file with an `ELF` column. `member` must be the `paramSet` of the bySet scores; replicates (`rep`) are reduced to their median."),
                    cacheKey = FALSE), # only .inputObjects (never cached) reads the files; what it read is in init's key
    defineParameter("weightsBySet", "character", NA_character_, NA, NA,
                    "Mode 'weights': the bySet scores of the same members (`wave1_bySet.csv`); the report rows come from it. NA: wave1_bySet.csv in outputPath.",
                    cacheKey = FALSE), # only .inputObjects (never cached) reads the files; what it read is in init's key
    defineParameter("weightsGrid", "character", NA_character_, NA, NA,
                    "Mode 'weights': csv of weight vectors fixed in advance, at most 50 rows, one column per objective term weighted (term names, not parameters; optional `weightId`).",
                    cacheKey = FALSE), # only .inputObjects (never cached) reads the files; what it read is in init's key
    defineParameter("weightsTuneELFs", "character", character(0), NA, NA, "Mode 'weights': ELFs the weights are selected on. Empty: every ELF of the terms that is not reported."),
    defineParameter("weightsReportELFs", "character", character(0), NA, NA,
                    "Mode 'weights': ELFs the result is reported on, none of them tuned on. Empty is allowed only with `weightsLeaveOneOut` (then the tune ELFs)."),
    defineParameter("weightsLeaveOneOut", "logical", FALSE, NA, NA, "Mode 'weights': hold out each report ELF in turn; its result uses weights selected on the other tune ELFs."),
    defineParameter("weightsLoss", "character", "area_ratio", NA, NA,
                    "Mode 'weights': bySet column the weights are selected on, as the median over the tune ELFs: distance from 1 (log scale) for a `_ratio` column, else maximised."),
    defineParameter("weightsTuneBySet", "character", NA_character_, NA, NA,
                    "Mode 'weights': bySet scores of other years to select on (to tune and report on the same ELFs). NA: `weightsBySet`.",
                    cacheKey = FALSE), # only .inputObjects (never cached) reads the files; what it read is in init's key
    defineParameter("skipELFs", "character", character(0), NA, NA, "ELFs of hindcastPath left out."),
    defineParameter("inputsPath", "character", NA_character_, NA, NA,
                    "Folder of the ELF input folders: <inputsPath>/<ELF>/ has rtml_<ELF>.tif, the spread-fit ledger fireSenseParams_*<ledgerSuffix>.rds and the NFDB_point_*.shp and NBAC_merged_*.shp observations (the latest of each).",
                    cacheKey = FALSE), # as hindcastPath
    defineParameter("ledgerSuffix", "character", "", NA, NA, "Suffix of the spread-fit ledger, e.g. '_crpsCauseAnnD03WholeLong_20261108'."),
    defineParameter("nCores", "integer", 1L, 1L, NA, "ELFs read at once (forked).", cacheKey = FALSE), # how many read at once, not what is read
    ## the page: <outputPath>/validation.html, beside wave1.rds and wave1_results.csv
    defineParameter("pageTitle", "character", "Hindcast validation", NA, NA, "Title of validation.html.",
                    cacheKey = FALSE), # used only by writeOutputs
    defineParameter("pageSubtitle", "character", "", NA, NA, "Subtitle of validation.html.",
                    cacheKey = FALSE), # used only by writeOutputs
    defineParameter("pageDpi", "numeric", 50, 1, NA, "Resolution of the page's figures.",
                    cacheKey = FALSE) # used only by writeOutputs
  ),
  inputObjects = bindrows(
    expectsInput("validationELFs", "list",
                 "One plain list per ELF: elf, simYears, simSizes, obsYears, obsSizes, covariate, fit (see runWave1()). A simulator-specific reader builds them; for fireSense, readFireSenseHindcast()."),
    expectsInput("anovaLong", "data.frame",
                 "Mode 'anova': the scores of the arms with their factors, readAnovaScores() of `anovaScores`, `anovaFactors` and `anovaMetrics`; read by .inputObjects, so that the file content is in init's cache key."),
    expectsInput("weightsTables", "list",
                 "Mode 'weights': `terms`, `bySet`, `grid` and `tuneBySet` (NULL when not given), read by .inputObjects from `weightsTerms`, `weightsBySet`, `weightsGrid` and `weightsTuneBySet`, so that the file content is in init's cache key.")
  ),
  outputObjects = bindrows(
    createsOutput("anovaResults", "list", "Mode 'anova': `terms` (per metric and term: df, pctSS, dropOnePctSS, F, p, note), `cells` and `long`, from computeAnova(); the first two are also written to <outputPath>/anova_terms.csv and anova_cells.csv."),
    createsOutput("weightsResults", "list", "Mode 'weights': `report`, `tuning` and `losses`, from tuneWeights(); also written to <outputPath>/weights_report.csv, weights_tuning.csv and weights_losses.csv."),
    createsOutput("bySetResults", "data.frame", "Mode 'score' with bySet: scoreBySet() of every ELF, written to <outputPath>/wave1_bySet.csv. NULL otherwise."),
    createsOutput("validationResults", "list", "runWave1() result: results, gates, bias, climate, perELF, flagged. Also saved to <outputPath>/wave1.rds, its results table to wave1_results.csv, and the page to validation.html.")
  )
))

## the size scored: minHa, or the escape size when it is NA
scoredMinHa <- function(p) if (is.na(p$minHa)) p$escapeSizeHa else p$minHa

.inputObjects <- function(sim) {
  fireSenseUtils::escapeDefinitionOf(sim)     # stops if another module has another escapeSizeHa or fireCauses
  ## fireSense hindcasts: every ELF that has finished all its years; each ELF's read is cached in <outputPath>/elfs, the observed records in <outputPath>/obs
  if (P(sim)$mode == "score" && !suppliedElsewhere("validationELFs", sim) && !is.na(P(sim)$hindcastPath)) {
    hc <- P(sim)$hindcastPath; yrs <- P(sim)$hindcastYears; sets <- P(sim)$hindcastSets; inp <- P(sim)$inputsPath
    obsDir <- file.path(outputPath(sim), "obs"); dir.create(obsDir, recursive = TRUE, showWarnings = FALSE)
    elfDir <- file.path(outputPath(sim), "elfs"); dir.create(elfDir, recursive = TRUE, showWarnings = FALSE)
    done <- Filter(function(e) length(hindcastFinishedYears(hc, e, yrs, sets)) == length(yrs),
                   setdiff(list.files(hc), P(sim)$skipELFs))
    sim$validationELFs <- readHindcastELFs(done, hc, years = yrs, sets = sets, inputsDir = inp, ledgerSuffix = P(sim)$ledgerSuffix,
                                           covariateDir = if (is.na(P(sim)$covariatePath)) hc else P(sim)$covariatePath,
                                           obsCacheDir = obsDir, elfCacheDir = elfDir, fireCauses = P(sim)$fireCauses,
                                           escapeSizeHa = P(sim)$escapeSizeHa, minHa = scoredMinHa(P(sim)), nCores = P(sim)$nCores)
    if (!length(sim$validationELFs))
      stop("no finished hindcasts in ", hc, ": no ELF has every one of hindcastYears ", min(yrs), "-", max(yrs), " in hindcastSets ", paste(sets, collapse = ", "),
           if (length(P(sim)$skipELFs)) paste0(" (skipELFs: ", paste(P(sim)$skipELFs, collapse = ", "), ")"))
    message("scoring ", length(sim$validationELFs), " ELFs")
  }
  if (P(sim)$mode == "anova" && !suppliedElsewhere("anovaLong", sim)) {
    p <- P(sim)
    sim$anovaLong <- readAnovaScores(p$anovaScores, p$anovaFactors, if (!anyNA(p$anovaMetrics)) p$anovaMetrics)
  }
  if (P(sim)$mode == "weights" && !suppliedElsewhere("weightsTables", sim)) {
    p <- P(sim)
    sim$weightsTables <- readWeightsTables(p$weightsTerms, if (is.na(p$weightsBySet)) file.path(outputPath(sim), "wave1_bySet.csv") else p$weightsBySet,
                                           p$weightsGrid, p$weightsTuneBySet)
  }
  return(invisible(sim))
}

doEvent.fireSimulator_validation <- function(sim, eventTime, eventType) {
  switch(eventType,
    init = {
      p <- P(sim)
      ## init only computes; the files are written by `writeOutputs`, because a cached init restores only the declared outputs
      if (p$mode == "anova") {
        sim$anovaResults <- computeAnova(sim$anovaLong)
      } else if (p$mode == "weights") {
        orNULL <- function(x) if (length(x)) x
        w <- sim$weightsTables
        sim$weightsResults <- tuneWeights(w$terms, w$bySet, w$grid, tuneBySet = w$tuneBySet, tuneELFs = orNULL(p$weightsTuneELFs),
                                          reportELFs = orNULL(p$weightsReportELFs), leaveOneOut = p$weightsLeaveOneOut, loss = p$weightsLoss)
      } else {
        sim$validationResults <- runWave1(sim$validationELFs, minHa = scoredMinHa(p), nBoot = p$nBoot, nNull = p$nNull, nPerm = p$nPerm, seed = p$seed)
        sim$bySetResults <- if (p$bySet) do.call(rbind, lapply(sim$validationELFs, scoreBySet, nBoot = p$nBoot, seed = p$seed))
      }
      params(sim)$fireSimulator_validation$.neverCache <- union(p$.neverCache, "writeOutputs")   # a user's `.neverCache` replaces the default; writeOutputs is never cached
      sim <- scheduleEvent(sim, start(sim), "fireSimulator_validation", "writeOutputs", eventPriority = 6)
    },
    writeOutputs = if (P(sim)$mode == "anova") {
      writeAnova(sim$anovaResults, outputPath(sim))
      Plots(data = sim$anovaResults$terms, fn = plotAnova, filename = "anova", ggsaveArgs = anovaSize)
    } else if (P(sim)$mode == "weights") {
      writeWeights(sim$weightsResults, outputPath(sim))
    } else {
      if (!all(is.na(P(sim)$.plots))) {                                  # the figures are built only when they are drawn
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
      }
      out <- outputPath(sim); dir.create(out, recursive = TRUE, showWarnings = FALSE)
      saveRDS(sim$validationResults, file.path(out, "wave1.rds"))
      utils::write.csv(sim$validationResults$results, file.path(out, "wave1_results.csv"), row.names = FALSE)
      if (P(sim)$bySet) utils::write.csv(sim$bySetResults, file.path(out, "wave1_bySet.csv"), row.names = FALSE)
      buildValidationPage(sim$validationELFs, sim$validationResults$gates, file.path(out, "validation.html"),
                          template = file.path(modulePath(sim), currentModule(sim), "pageTemplate.html"),
                          title = P(sim)$pageTitle, subtitle = P(sim)$pageSubtitle, label = P(sim)$simLabel,
                          minHa = scoredMinHa(P(sim)), dpi = P(sim)$pageDpi, bias = sim$validationResults$bias, climate = sim$validationResults$climate)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}
