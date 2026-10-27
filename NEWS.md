# fireSimulator_validation 0.0.0.9000

* `readHindcastSim()` no longer fails when a simulated set-year has no fires (a header-only `hindcast_sizes.csv`); that set-year reads as zero area.
* Wave 1 of the hindcast test plan: gates G1-G6, skill tests S1-S5, nulls N0-N2, pooled bootstrap, `runWave1()`, and a fireSense hindcast reader.
