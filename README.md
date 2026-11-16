# fireSimulator_validation

A SpaDES module that scores a fire simulator against observed fire. The scoring functions in `R/` are
model-agnostic: they take plain tables (and rasters) only, so any simulator can use them. The module file
calls `runWave1()` on `validationELFs` and writes `wave1.rds`, `wave1_results.csv` and `validation.html` to its
`outputPath`. A fireSense hindcast needs no `validationELFs`: set `hindcastPath`, `hindcastYears`, `inputsPath` and
`ledgerSuffix`, and `.inputObjects` reads every ELF that has finished all its years.

Wave 1 of the hindcast test plan:

- Gates G1-G6 (`R/gates.R`): pass/fail, no null.
- Skill tests S1-S5 (`R/skillTests.R`): CRPS skill on log(1 + x) against N1 (the ELF's own history) and N2 (a
  regression on the ELF-mean climate covariate), N0 (shuffled years) as a floor; pooled over ELFs with a bootstrap
  over ELFs and years (`R/bootstrap.R`).
- `runWave1()` (`R/runWave1.R`) takes a list of ELF inputs and returns a tidy results table, the gates and per-ELF tables.
- `runTablesToELF()` (`R/readRunTables.R`) builds the same input from plain per-run tables (fire sizes, ignitions per run-year), and `buildValidationPage()` (`R/buildPage.R`, layout `pageTemplate.html`) writes the figures and the scores table as one HTML file.
- `R/readFireSenseHindcast.R` is the only fireSense-specific file: it turns FireSenseTesting hindcast outputs, the
  spread-fit ledger and NFDB/NBAC observations into the ELF input list.
- Figures (`R/plots.R`): `plotAreaByYear()`, `plotAnnualAreaDensity()`, `plotFireSizes()`, `plotIgnitions()`, `plotEscapes()`, one ELF each; `plotELFPanel()` puts them in one figure per ELF and `plotMetricPanel()` lays one metric out for all ELFs; `plotOverview()` shows the five gate ratios of all ELFs. The module draws the panels and the overview after scoring through `SpaDES.core::Plots()`; parameter `.plots` (`"screen"`, `"png"`) decides where.

- Weight tuning (`R/weights.R`, parameter `mode = "weights"`): `tuneWeights()` picks, for each weight vector of a fixed grid and each ELF, the spread-fit member with the lowest weighted objective terms, reports that member's bySet indicators, and selects the vector on tune ELFs; the result table has the held-out ELFs only. Terms table: `ELF`, `member` (= `paramSet`), one column per raw term.
- Analysis of variance (`R/anova.R`, parameter `mode = "anova"`): `runAnova()` reads the `wave1_results.csv` of several score folders (one arm each) and a table of the arms' factors, and `anovaTerms()` partitions the variation of each per-ELF skill score among ELF, the factors and their interactions (`anova_terms.csv`, `anova_cells.csv` and a figure); aliased terms are reported, not dropped.

Tests: `testthat::test_dir("tests/testthat")`.
