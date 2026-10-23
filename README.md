# fireSimulator_validation

A SpaDES module that scores a fire simulator against observed fire. Metrics are model-agnostic: the scoring
functions in `R/` take plain tables only. `R/readFireSenseHindcast.R` is the one fireSense-specific file.
