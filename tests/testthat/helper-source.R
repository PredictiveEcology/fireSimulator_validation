## The scoring functions are plain functions in R/; source them as simInit() would (reqdPkgs of the plots are attached by simInit()).
suppressMessages({library(ggplot2); library(patchwork)})
for (f in list.files(testthat::test_path("..", "..", "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

toyMinHa <- 50                                  # the toy escape size: every toy fire is 50 ha or more

## A toy ELF in the runWave1() input format. `perfect = TRUE`: every simulated run is the observed year;
## otherwise runs are drawn from the whole record (the climatological forecast).
toyElf <- function(name = "toy", perfect = FALSE, seed = 1, scored = seq(1986L, 2021L, 5L), nRuns = 20L, years = 1986:2024) {
  set.seed(seed)
  cov <- rnorm(length(years))
  ign <- rnbinom(length(years), mu = exp(3 + 0.4 * cov), size = 5)
  esc <- rbinom(length(years), ign, plogis(-1.5 + 0.5 * cov))
  sizes <- lapply(seq_along(years), function(i) toyMinHa * exp(rexp(esc[i], 1 / 1.2)))
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

## A toy ELF with a known climate response: observed area = exp(5 + b z) in each scored year; run r is `f(r)` times a response `b * s` (slope factor s)
climToy <- function(b = 0.5, s = 1, f = function(r) 1, nRuns = 5L, years = 1986:2024) {
  z <- as.numeric(scale(sin(seq_along(years) * 1.7) + seq_along(years) / 40))   # distinct values: no ties at the quantiles
  obs <- exp(5 + b * z)
  sim <- do.call(rbind, lapply(seq_len(nRuns), function(r) data.frame(year = years, run = r, ignitions = 1, escapes = 1, area50 = f(r) * exp(5 + b * s * z))))
  list(elf = "clim", simYears = sim, simSizes = data.frame(year = years, run = 1L, ha = 100),
       obsYears = data.frame(year = years, ignitions = 1, escapes = 1, area50 = obs), obsSizes = data.frame(year = years, ha = 100),
       covariate = data.frame(year = years, value = 100 + 20 * z), z = z, obs = obs)
}

## Run the module from a temporary root: `params` are its parameters, outputs go to <root>/outputs. `.env` is where the options are undone.
runModule <- function(root, params, inputPath = root, .env = parent.frame()) {
  mp <- file.path(root, "modules"); dir.create(mp, showWarnings = FALSE)
  link <- file.path(mp, "fireSimulator_validation")
  if (!file.exists(link)) file.symlink(normalizePath(testthat::test_path("..", "..")), link) # the checkout's folder name may differ; a second link would be made inside the checkout
  withr::local_options(spades.useRequire = FALSE, spades.moduleCodeChecks = FALSE, .local_envir = .env)
  SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = 1), modules = "fireSimulator_validation", params = list(fireSimulator_validation = params),
    paths = list(modulePath = mp, outputPath = file.path(root, "outputs"), inputPath = inputPath,
                 cachePath = file.path(root, "cache"), scratchPath = file.path(root, "scratch")))
}
## the cache entries (one per cached event run) in <root>/cache
nCached <- function(root) length(unique(reproducible::showCache(file.path(root, "cache"), verbose = -1)$cacheId))
