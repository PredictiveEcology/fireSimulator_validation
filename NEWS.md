# fireSimulator_validation (development version)

* Hindcast figures come from the module: `plotAreaByYear()`, `plotAnnualAreaDensity()`, `plotFireSizes()`, `plotIgnitions()`, `plotEscapes()` (density curves; one simulated series against the observed record, medians, single-run envelopes), `plotELFPanel()` (one patchwork figure per ELF), `plotMetricPanel()` (one metric, a pane per ELF), `elfPlots()` and `plotOverview()` (the five gate ratios of all ELFs), drawn after scoring through `SpaDES.core::Plots()`, so the new `.plots` parameter decides screen or png. New `summariseGates()` gives the gate ratios as one row per ELF. New parameters `.plots`, `.plotInitialTime`, `simLabel`, `minHa`; `ggplot2`, `patchwork` and `scales` are in `reqdPkgs`.

* Gates G2, G4 and G6 report their run-to-run uncertainty: new columns `lo` and `hi` (2.5% and 97.5% of a bootstrap that resamples whole runs, with replacement, within each year, `nBoot = 1000`, `seed = 1`) beside `value` in `gates`, and "0.72 (0.63-1.03)" in `detail`. The point value and the resamples go through the same statistic (`bootRuns()`). Pass/fail logic unchanged.

* Gates G2, G4 and G6 also report the spread of single runs: new columns `runLo` and `runHi` (2.5% and 97.5% of the gate's ratio computed for each run on its own, across runs), and "0.72 (median-run 95%: 0.63-1.03; single runs 95%: 0.31-1.85)" in `detail`, so an observed record can be judged against what one run could produce. Computed in `bootRuns()` through the same statistic as `value`. `lo`, `hi` and pass/fail unchanged.

# fireSimulator_validation 0.0.0.9000

* `readHindcastSim()` no longer fails when a simulated set-year has no fires (a header-only `hindcast_sizes.csv`); that set-year reads as zero area.
* G4 compares the median of the 10 largest fires (each run's, median over runs) with the observed median of the 10 largest, not the single largest fire: one fire in 39 years was too noisy a target. The single largest fire stays in `detail` for information. New internal `topFires()`. Pass rule unchanged (within 10x).

* Wave 1 of the hindcast test plan: gates G1-G6, skill tests S1-S5, nulls N0-N2, pooled bootstrap, `runWave1()`, and a fireSense hindcast reader.
