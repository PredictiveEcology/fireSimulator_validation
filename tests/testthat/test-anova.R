## Analysis of variance of scores across arms, on synthetic designs where the answer is known.
## Arms = every combination of A (large effect), B (no effect) and C; ELFs differ a little.
mkAnovaToy <- function(confound = FALSE, seed = 1) {
  set.seed(seed)
  arms <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2"), C = c("c1", "c2"), stringsAsFactors = FALSE)
  if (confound) arms$C <- ifelse(arms$A == "a1", "c1", "c2")           # C varies with A in every arm
  arms$arm <- paste0("arm", seq_len(nrow(arms)))
  elfs <- c("e1", "e2", "e3")
  d <- merge(data.frame(ELF = elfs), arms)
  d$value <- ifelse(d$A == "a2", 1, 0) + c(e1 = 0, e2 = 0.1, e3 = -0.1)[d$ELF] + rnorm(nrow(d), sd = 0.02)
  d$metric <- "S1_N1_skill"
  d
}
anovaFactorCols <- c("A", "B", "C")

test_that("a factor with a large effect holds most of the sum of squares and a factor with none holds almost none", {
  t <- anovaTerms(mkAnovaToy(), anovaFactorCols)
  expect_gt(t$pctSS[t$term == "A"], 80)
  expect_lt(t$pctSS[t$term == "B"], 1)
  expect_lt(t$dropOnePctSS[t$term == "B"], 1)
  expect_gt(t$dropOnePctSS[t$term == "A"], 50)
  expect_lt(t$p[t$term == "A"], 1e-6)
  expect_gt(t$p[t$term == "B"], 0.01)
  expect_equal(sum(t$pctSS), 100)
  expect_true(grepl("not run-to-run noise", t$note[t$term == "Residual"]))
})

test_that("a confounded design is reported as confounded, not dropped silently", {
  t <- anovaTerms(mkAnovaToy(confound = TRUE), anovaFactorCols)
  expect_equal(t$df[t$term == "C"], 0)
  expect_match(t$note[t$term == "C"], "aliased with A")
  expect_true("C" %in% t$term)
  expect_equal(sum(t$pctSS), 100)
})

test_that("a factor with one level is reported as constant", {
  d <- mkAnovaToy(); d$C <- "c1"
  t <- anovaTerms(d, anovaFactorCols)
  expect_match(t$note[t$term == "C"], "constant")
})

test_that("an ELF missing from some arms is reported as an aliased interaction", {
  d <- mkAnovaToy(); d <- d[!(d$ELF == "e3" & d$A == "a2"), ]
  t <- anovaTerms(d, anovaFactorCols)
  expect_match(t$note[t$term == "ELF:A"], "aliased")
})

## score folders as mode "score" writes them: wave1_results.csv with level "ELF" rows
mkScoreFolders <- function(root, d) {
  dirs <- vapply(unique(d$arm), function(a) {
    x <- d[d$arm == a, ]; f <- file.path(root, paste0("score_", a)); dir.create(f)
    write.csv(data.frame(test = "S1", null = "N1", level = "ELF", ELF = x$ELF, measure = "skill", value = x$value, lo = NA, hi = NA, nELF = 1L, nUnits = 39L, verdict = NA),
              file.path(f, "wave1_results.csv"), row.names = FALSE); f
  }, "")
  stats::setNames(dirs, unique(d$arm))
}

test_that("score folders of one arm in several folders are combined, and a missing arm in the factor table stops", {
  d <- mkAnovaToy(); root <- withr::local_tempdir()
  one <- mkScoreFolders(root, d)
  fac <- unique(d[c("arm", "A", "B", "C")])
  long <- readAnovaScores(one, fac)
  expect_equal(nrow(long), nrow(d)); expect_equal(unique(long$metric), "S1_N1_skill")
  expect_error(readAnovaScores(one, fac[-1, ]), "no row in the factor table")
  r2 <- withr::local_tempdir(); halves <- c(mkScoreFolders(r2, d[d$ELF == "e1", ]), stats::setNames(mkScoreFolders(withr::local_tempdir(), d[d$ELF != "e1", ]), unique(d$arm)))
  expect_equal(nrow(readAnovaScores(halves, fac)), nrow(d))
})

test_that("module parameter mode = 'anova' writes both tables and the figure, without a hindcast", {
  skip_if_not_installed("SpaDES.core")
  d <- mkAnovaToy(); root <- withr::local_tempdir(); out <- file.path(root, "outputs"); dir.create(out)
  sc <- mkScoreFolders(root, d); write.csv(unique(d[c("arm", "A", "B", "C")]), file.path(root, "factors.csv"), row.names = FALSE)
  mp <- file.path(root, "modules"); dir.create(mp)
  file.symlink(normalizePath(testthat::test_path("..", "..")), file.path(mp, "fireSimulator_validation"))
  withr::local_options(spades.useRequire = FALSE, spades.moduleCodeChecks = FALSE)
  sim <- SpaDES.core::simInitAndSpades(times = list(start = 0, end = 1), modules = "fireSimulator_validation",
    params = list(fireSimulator_validation = list(mode = "anova", anovaScores = sc, anovaFactors = file.path(root, "factors.csv"), .plots = "png")),
    paths = list(modulePath = mp, outputPath = out, inputPath = root, cachePath = file.path(root, "cache"), scratchPath = file.path(root, "scratch")))
  tm <- read.csv(file.path(out, "anova_terms.csv")); cl <- read.csv(file.path(out, "anova_cells.csv"))
  expect_named(tm, c("metric", "term", "df", "pctSS", "dropOnePctSS", "F", "p", "note"))
  expect_gt(tm$pctSS[tm$term == "A"], 80)
  expect_true(all(c("metric", "ELF", "A", "B", "C", "n", "median") %in% names(cl)))
  expect_equal(sim$anovaResults$terms$pctSS, tm$pctSS, tolerance = 1e-8)
  expect_length(list.files(out, pattern = "^anova.*\\.png$", recursive = TRUE), 1)
  expect_false(file.exists(file.path(out, "wave1.rds")))
})
