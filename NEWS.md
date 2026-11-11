# fireSimulator_validation (development version)

* `observedRecord()` is a `reproducible::Cache()` in `cachePath` (replacing the hand-made `cacheFile` named by its arguments), keyed on its arguments and on `fireSenseUtils::escapeFireSetDeps()`, so a record made by older escape code is no longer read back. `hindcastBurnMap()` reads the hindcast burn map once for `hindcastCore()` and `hindcastFlammable()` (which now takes it); the core is `fireSenseUtils::inELFCore()`, and `.inputObjects` checks the escape definition with `fireSenseUtils::escapeDefinitionOf()`. Needs fireSenseUtils #222.

* Scoring uses the shared escape definition and the core-start rule of fireSenseUtils (needs the `escape-core-perfire` changes: `hindcastObserved(coreStart = TRUE)`). New parameters `escapeSizeHa` (default `fireSenseUtils::fireSenseEscapeSizeHa`, 50) and `fireCauses` (default `fireSenseUtils::fireSenseFireCauses`), checked against the other modules' with `paramCheckOtherMods()` (error on a difference), passed to `readFireSenseHindcast()` (`.inputObjects` did not pass `fireCauses`) and recorded as `escapeSizeHa` and `fireCauses` on the observed record and the ELF input. The observed record counts fires that start in the core at their whole extent, with a `flammableRTM` of core and buffer (`hindcastFlammable()`); `observedRecord()` takes it as its second argument, and its cache file name has `_coreStart_<causes>_<escapeSizeHa>ha` (old `_burnable` caches are not read). Observed `escapes` are the escape set's, one per core-start escape polygon, no longer capped at the ignitions.
* `minHa` is honoured. It was passed only to the figures, while the readers and the size tests used 50: now `readHindcastSim()` counts escapes and area, `readHindcastObserved()`, `runTablesToELF()`, `scoreS4()` (through `runWave1(minHa =)`), `scoreSizesVsNull()`, `shuffledSizesAD()`, the plots and `buildValidationPage()` all take it, with no default. The parameter `minHa` now defaults to `NA`: the escape size, `escapeSizeHa`.

* The module scores fireSense hindcasts on its own: with `hindcastPath` set and no `validationELFs` supplied, `.inputObjects` reads every ELF of `hindcastPath` that has finished all `hindcastYears` with `readFireSenseHindcast()` (parameters `hindcastSets`, `skipELFs`, `inputsPath`, `ledgerSuffix`, `nCores`; observed-record caches in `<outputPath>/obs`). `init` writes `wave1.rds`, `wave1_results.csv` and `validation.html` (`pageTitle`, `pageSubtitle`, `pageDpi`) to `outputPath`. This replaces the scoring code in FireSenseTesting's `global.R`. Adds `jsonlite`, `terra` and `fireregimetools` to `reqdPkgs` (already used).

* The per-ELF panel figures are written as `panel_<ELF>` with the ELF's dots as underscores (`panel_11_3.png`): `Plots()` took `.3` of `11.3` as an extension, so 11.3 and 11.4 both wrote `panel_11.png`.

* The hindcast validation page orders ELFs by observed burn rate, most burned at the top: the scores table has a "Burned %/yr" column (the default sort, descending), the ELF picker and the first ELF shown follow it, and `elfOrder()` (so `plotOverview()` and `plotMetricPanel()`) orders by `burnRate` instead of the area ratio. `burnRate()` is the mean annual observed area (`area50`, natural fires of 50 ha or more, as the area plots) over the burnable area of the core, in % per year; `readFireSenseHindcast()` adds `burnableHa` and `burnRate` to the ELF input. An ELF input without `burnRate` sorts last.

* `observedRanks()` (the page's "obs N%" and observed-rank chart) compares the runs with the observed record of the scored years only. It summed the whole record (1986-2024) against the runs' 8 hindcast years, so the observed record ranked above nearly every run.

* The per-metric figures of `buildValidationPage()` (about 4 inches per row of three ELFs) are saved with `limitsize = FALSE`; with more than 36 ELFs ggsave refused them (over 50 inches) and the page was not written.

* Needs fireSenseUtils >= 0.2.3.9139: the observed-fire overlay of a very fragmented ELF core (5.1.2, 5.1.3) took hours per ELF before.

* `buildValidationPage()` takes `title` (default `"Hindcast validation"`), used for the browser tab and the heading; the template no longer hard-codes "burnable fits".

* `hindcastCore()` finds `rtml_<ELF>.tif` directly in `<inputsDir>/<ELF>/` (the GPSC layout) as well as in `<ELF>/NRV/370/1/`; scoring on the GPSC stopped with a missing file before.

* `runTablesToELF()` builds an ELF input from plain per-run tables (fire sizes and ignitions per run-year), for simulators that do not write FireSense hindcast folders, e.g. a fast hindcast. `buildValidationPage()` writes the scores table, the observed-rank chart, the overview and the per-metric and per-ELF figures as one self-contained HTML file from `pageTemplate.html`.

* `fireCauses` of `observedRecord()` and `readFireSenseHindcast()` now defaults to `fireSenseUtils::fireSenseFireCauses` (`c("natural", "unknown")`), the default of `fireSense_dataPrepFit` and `fireSense_summary`, instead of `NULL` (NFDB L/N, NBAC Natural), so the observed record has the causes the model was fitted to. The `"natural"` label of the observed sizes names the fires of these causes. `NULL` still gives the old record. Observed-record caches of the default now carry `_natural-unknown` in the file name, so old default caches are not read.

* Gates G2, G4 and G6 report where the observed record ranks among the single runs: `bootRuns()` also returns `single` (the statistic of each run), `obsRank` (mid-rank percentile of the observed, 1, among the runs: 0 = every run above it) and `nRuns`; `gates` gains the columns `obsRank` and `nRuns` (G2, G4; NA elsewhere), and every `detail` ratio ends "; observed above 35% of 20 runs" (G6: one per link), after the unchanged interval text. Values and pass/fail unchanged.

* `observedRecord()` and `readFireSenseHindcast()` take `fireCauses` (default `NULL`: unchanged, NFDB L/N and NBAC Natural), passed to `fireSenseUtils::hindcastObserved()`; e.g. `c("natural", "unknown")` scores the causes the model was fitted to. The groups are in the observed-record cache file's name, so records of other causes are not read. Needs fireSenseUtils with `fireCauseCodes()` (Depends-on: PredictiveEcology/fireSenseUtils#193).

* `observedRecord()` passes `flammableRTM = core` to `fireSenseUtils::hindcastObserved()`, so observed area burned and fire sizes of the unscored years are counted on the burnable cells of the core, as the scored years of the hindcast are. Its cache file name gets `_burnable` (`<ELF>_burnable.rds`), so caches of the whole-polygon count are not reused and are rebuilt. Needs fireSenseUtils#190.

* Hindcast figures come from the module: `plotAreaByYear()`, `plotAnnualAreaDensity()`, `plotFireSizes()`, `plotIgnitions()`, `plotEscapes()` (density curves; one simulated series against the observed record, medians, single-run envelopes), `plotELFPanel()` (one patchwork figure per ELF), `plotMetricPanel()` (one metric, a pane per ELF), `elfPlots()` and `plotOverview()` (the five gate ratios of all ELFs), drawn after scoring through `SpaDES.core::Plots()`, so the new `.plots` parameter decides screen or png. New `summariseGates()` gives the gate ratios as one row per ELF. New parameters `.plots`, `.plotInitialTime`, `simLabel`, `minHa`; `ggplot2`, `patchwork` and `scales` are in `reqdPkgs`.

* Gates G2, G4 and G6 report their run-to-run uncertainty: new columns `lo` and `hi` (2.5% and 97.5% of a bootstrap that resamples whole runs, with replacement, within each year, `nBoot = 1000`, `seed = 1`) beside `value` in `gates`, and "0.72 (0.63-1.03)" in `detail`. The point value and the resamples go through the same statistic (`bootRuns()`). Pass/fail logic unchanged.

* Gates G2, G4 and G6 also report the spread of single runs: new columns `runLo` and `runHi` (2.5% and 97.5% of the gate's ratio computed for each run on its own, across runs), and "0.72 (median-run 95%: 0.63-1.03; single runs 95%: 0.31-1.85)" in `detail`, so an observed record can be judged against what one run could produce. Computed in `bootRuns()` through the same statistic as `value`. `lo`, `hi` and pass/fail unchanged.

# fireSimulator_validation 0.0.0.9000

* `readHindcastSim()` no longer fails when a simulated set-year has no fires (a header-only `hindcast_sizes.csv`); that set-year reads as zero area.
* G4 compares the median of the 10 largest fires (each run's, median over runs) with the observed median of the 10 largest, not the single largest fire: one fire in 39 years was too noisy a target. The single largest fire stays in `detail` for information. New internal `topFires()`. Pass rule unchanged (within 10x).

* Wave 1 of the hindcast test plan: gates G1-G6, skill tests S1-S5, nulls N0-N2, pooled bootstrap, `runWave1()`, and a fireSense hindcast reader.
