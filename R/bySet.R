## Scoring per parameter set: the indicators of the gates and of the annual-area test for each set on its own.

#' One ELF's input restricted to one parameter set
#'
#' @param elf An ELF input list ([runWave1()]) whose `simYears` and `simSizes` have a `paramSet` column ([readHindcastSim()]).
#' @param set The `paramSet` to keep. @return The same list with only that set's runs.
elfOfSet <- function(elf, set) {
  elf$simYears <- elf$simYears[elf$simYears$paramSet == set, , drop = FALSE]
  elf$simSizes <- elf$simSizes[elf$simSizes$paramSet == set, , drop = FALSE]
  elf
}

#' Score each parameter set of an ELF on its own
#'
#' For every `paramSet` of `elf$simYears`: the gates ([runGates()]) and the annual-area skill ([scoreS3()], [skillScore()]), on that set's runs only.
#' The ratios are those of [summariseGates()]: the median run over the observation (area, median of the 10 largest fires, ignitions, escape rate,
#' area per escape), each with the 2.5% and 97.5% of single runs (`_runLo`, `_runHi`).
#'
#' @param elf An ELF input list with a `paramSet` column in `simYears` and `simSizes`. @param nBoot,seed As [runGates()].
#' @return `data.frame`, one row per set: `ELF`, `paramSet`, `nRuns`, the ratio columns of [summariseGates()], `annualArea_skillN1`,
#'   `annualArea_skillN2` (skill of the set's annual area against the null N1 and N2), and `G2`, `G4`, `G6` (pass).
scoreBySet <- function(elf, nBoot = 1000, seed = 1) {
  do.call(rbind, lapply(sort(unique(elf$simYears$paramSet)), function(s) {
    e <- elfOfSet(elf, s)
    g <- cbind(ELF = e$elf, runGates(e, nBoot = nBoot, seed = seed))
    sk <- scoreS3(e)$scores
    skill <- function(nl) { x <- sk[sk$null == nl, ]; skillScore(x$lossModel, x$lossNull) }
    pass <- stats::setNames(lapply(c("G2", "G4", "G6"), function(id) g$pass[g$gate == id]), c("G2", "G4", "G6"))
    cbind(summariseGates(g), paramSet = s, nRuns = length(unique(e$simYears$run)),
          annualArea_skillN1 = skill("N1"), annualArea_skillN2 = skill("N2"), as.data.frame(pass))
  }))
}
