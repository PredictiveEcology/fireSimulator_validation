## Everything in this file and any files in the R directory are sourced during `simInit()`.
## The scoring functions in R/ take plain tables only; the module file just hands its input to the runner.
defineModule(sim, list(
  name = "fireSimulator_validation",
  description = paste("Scores a fire simulator against observed fire: gates and skill tests on annual ignitions, escape rate,",
                      "area burned, fire sizes and ranking of years, against nulls fitted to the observations (wave 1 of the hindcast test plan)."),
  keywords = c("fire", "validation", "hindcast", "CRPS"),
  authors = c(person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(fireSimulator_validation = "0.0.0.9002"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  documentation = list("README.md"),
  reqdPkgs = list("data.table", "MASS", "withr", "PredictiveEcology/fireSenseUtils@development (>= 0.2.3.9118)"),
  parameters = bindrows(
    defineParameter("nBoot", "integer", 999L, 1L, NA, "Bootstrap resamples of ELFs and years for the pooled verdicts."),
    defineParameter("nNull", "integer", 199L, 1L, NA, "Draws per year of the Anderson-Darling null."),
    defineParameter("nPerm", "integer", 9999L, 1L, NA, "Permutations of the ranking floor."),
    defineParameter("seed", "integer", 1L, NA, NA, "Seed of the bootstrap and permutations.")
  ),
  inputObjects = bindrows(
    expectsInput("validationELFs", "list",
                 "One plain list per ELF: elf, simYears, simSizes, obsYears, obsSizes, covariate, fit (see runWave1()). A simulator-specific reader builds them; for fireSense, readFireSenseHindcast().")
  ),
  outputObjects = bindrows(
    createsOutput("validationResults", "list", "runWave1() result: results, gates, perELF, flagged.")
  )
))

doEvent.fireSimulator_validation <- function(sim, eventTime, eventType) {
  switch(eventType,
    init = {
      sim$validationResults <- runWave1(sim$validationELFs, nBoot = P(sim)$nBoot, nNull = P(sim)$nNull,
                                        nPerm = P(sim)$nPerm, seed = P(sim)$seed)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}
