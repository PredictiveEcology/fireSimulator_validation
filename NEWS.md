# fireSimulator_validation (development version)

* `observedRecord()` passes `flammableRTM = core` to `fireSenseUtils::hindcastObserved()`, so observed area burned and fire sizes of the unscored years are counted on the burnable cells of the core, as the scored years of the hindcast are. Its cache file name gets `_burnable` (`<ELF>_burnable.rds`), so caches of the whole-polygon count are not reused and are rebuilt. Needs fireSenseUtils#190.

# fireSimulator_validation 0.0.0.9000

* Wave 1 of the hindcast test plan: gates G1-G6, skill tests S1-S5, nulls N0-N2, pooled bootstrap, `runWave1()`, and a fireSense hindcast reader.
