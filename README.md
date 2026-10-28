# fireSimulator_validation

A SpaDES module that scores a fire simulator against observed fire. The scoring functions in `R/` are
model-agnostic: they take plain tables (and rasters) only, so any simulator can use them. The module file just
calls `runWave1()`.

Wave 1 of the hindcast test plan:

- Gates G1-G6 (`R/gates.R`): pass/fail, no null.
- Skill tests S1-S5 (`R/skillTests.R`): CRPS skill on log(1 + x) against N1 (the ELF's own history) and N2 (a
  regression on the ELF-mean climate covariate), N0 (shuffled years) as a floor; pooled over ELFs with a bootstrap
  over ELFs and years (`R/bootstrap.R`).
- `runWave1()` (`R/runWave1.R`) takes a list of ELF inputs and returns a tidy results table, the gates and per-ELF tables.
- `R/readFireSenseHindcast.R` is the only fireSense-specific file: it turns FireSenseTesting hindcast outputs, the
  spread-fit ledger and NFDB/NBAC observations into the ELF input list.
- Figures (`R/plots.R`): `plotAreaByYear()`, `plotAnnualAreaDensity()`, `plotFireSizes()`, `plotIgnitions()`, `plotEscapes()`, one ELF each; `plotELFPanel()` puts them in one figure per ELF; `plotOverview()` shows the five gate ratios of all ELFs. The module draws the panels and the overview after scoring through `SpaDES.core::Plots()`; parameter `.plots` (`"screen"`, `"png"`) decides where.

Tests: `testthat::test_dir("tests/testthat")`.
