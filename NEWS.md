# fireSimulator_validation 0.0.0.9000

* G4 compares the median of the 10 largest fires (each run's, median over runs) with the observed median of the 10 largest, not the single largest fire: one fire in 39 years was too noisy a target. The single largest fire stays in `detail` for information. New internal `topFires()`. Pass rule unchanged (within 10x).

* Wave 1 of the hindcast test plan: gates G1-G6, skill tests S1-S5, nulls N0-N2, pooled bootstrap, `runWave1()`, and a fireSense hindcast reader.
