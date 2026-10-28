## The scoring functions are plain functions in R/; source them as simInit() would (reqdPkgs of the plots are attached by simInit()).
suppressMessages({library(ggplot2); library(patchwork)})
for (f in list.files(testthat::test_path("..", "..", "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

## A toy ELF in the runWave1() input format. `perfect = TRUE`: every simulated run is the observed year;
## otherwise runs are drawn from the whole record (the climatological forecast).
toyElf <- function(name = "toy", perfect = FALSE, seed = 1, scored = seq(1986L, 2021L, 5L), nRuns = 20L, years = 1986:2024) {
  set.seed(seed)
  cov <- rnorm(length(years))
  ign <- rnbinom(length(years), mu = exp(3 + 0.4 * cov), size = 5)
  esc <- rbinom(length(years), ign, plogis(-1.5 + 0.5 * cov))
  sizes <- lapply(seq_along(years), function(i) 50 * exp(rexp(esc[i], 1 / 1.2)))
  obsSizes <- data.frame(year = rep(years, esc), ha = unlist(sizes))
  area <- vapply(sizes, sum, 0)
  obsYears <- data.frame(year = years, ignitions = ign, escapes = esc, area50 = area)
  pick <- function(y) if (perfect) match(y, years) else sample.int(length(years), 1L)
  runs <- do.call(rbind, lapply(scored, function(y) do.call(rbind, lapply(seq_len(nRuns), function(r) {
    i <- pick(y)
    data.frame(year = y, run = r, ignitions = ign[i], escapes = esc[i], area50 = area[i], idx = i)
  }))))
  simSizes <- do.call(rbind, lapply(seq_len(nrow(runs)), function(k) {
    ha <- sizes[[runs$idx[k]]]
    if (length(ha)) data.frame(year = runs$year[k], run = runs$run[k], ha = ha)
  }))
  list(elf = name, simYears = runs[, c("year", "run", "ignitions", "escapes", "area50")], simSizes = simSizes,
       obsYears = obsYears, obsSizes = obsSizes, covariate = data.frame(year = years, value = cov),
       fit = list(objective = c(100, 101, 102), pipelineObjective = c(100.5, 101.5, 99.5, 101), pc = c(0.2, 0.21), p50 = 0.25,
                  fitAreaHa = sum(area), totalAreaHa = sum(area) * 1.2))
}
