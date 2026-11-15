## Weight tuning on a synthetic terms + bySet table where the right weights are known.
## Member k (1..10) has terms A = k and B = (8 - k)^2 and area ratio exp((k - 8) / 10): member 8 is the right one (ratio 1).
## Weight (1, 0) picks k = 1, (0, 1) picks k = 8, (1, 1) picks k = 7 (A + B = 7 + 1 = 8 = 8 + 0: a tie, the first).
mkWeightsToy <- function(elfs = c("a", "b", "c"), best = 8L) {
  k <- 1:10
  terms <- do.call(rbind, lapply(elfs, function(e) data.frame(ELF = e, member = k, A = k, B = (best - k)^2)))
  bySet <- do.call(rbind, lapply(elfs, function(e) data.frame(ELF = e, paramSet = k, area_ratio = exp((k - best) / 10),
                                                              area_runLo = 0.5, area_runHi = 2, annualArea_skillN1 = -abs(k - best))))
  list(terms = terms, bySet = bySet)
}
grid3 <- data.frame(A = c(1, 0, 1), B = c(0, 1, 1))

test_that("each weight vector picks the member that minimises the weighted terms, and reports that member's bySet row", {
  d <- mkWeightsToy()
  out <- tuneWeights(d$terms, d$bySet, grid3, tuneELFs = c("a", "b"), reportELFs = "c")
  expect_equal(out$tuning$member[out$tuning$ELF == "a"], c(1L, 8L, 7L))
  expect_equal(out$tuning$area_ratio[out$tuning$ELF == "a"], exp((c(1, 8, 7) - 8) / 10))
  expect_true(all(c("area_runLo", "area_runHi") %in% names(out$tuning)))
})

test_that("the weights that hindcast best on the tune ELFs are selected, and only report-ELF rows are the result", {
  d <- mkWeightsToy()
  out <- tuneWeights(d$terms, d$bySet, grid3, tuneELFs = c("a", "b"), reportELFs = "c")
  expect_equal(out$losses$weightId[out$losses$selected], 2L)
  expect_equal(out$report$ELF, "c"); expect_equal(out$report$weightId, 2L); expect_equal(out$report$member, 8L)
  expect_equal(unique(out$report$set), "report"); expect_equal(unique(out$tuning$set), "tune")
  expect_false(any(out$tuning$ELF == "c"))
})

test_that("leave-one-ELF-out selects without the held-out ELF and reports only it", {
  d <- mkWeightsToy()
  d$bySet$area_ratio[d$bySet$ELF == "c"] <- 1        # ELF c scores 1 whatever the member: the held-out ELF cannot influence the choice
  out <- tuneWeights(d$terms, d$bySet, grid3, leaveOneOut = TRUE)
  expect_equal(sort(out$report$ELF), c("a", "b", "c")); expect_equal(unique(out$report$weightId), 2L)
  expect_equal(nrow(out$losses), 9L)
  expect_true(all(out$report$set == "report"))
})

test_that("tuning and reporting the same ELF stops, unless the selection scores are other years'", {
  d <- mkWeightsToy()
  expect_error(tuneWeights(d$terms, d$bySet, grid3, tuneELFs = "a", reportELFs = "a"), "leaveOneOut")
  other <- d$bySet; other$area_ratio <- rev(other$area_ratio)
  expect_no_error(tuneWeights(d$terms, d$bySet, grid3, tuneELFs = "a", reportELFs = "a", tuneBySet = other))
  expect_error(tuneWeights(d$terms, d$bySet, grid3, tuneELFs = "a"), "reportELFs")
})

test_that("a grid of more than 50 vectors, a missing term and an unknown loss stop", {
  d <- mkWeightsToy()
  expect_error(tuneWeights(d$terms, d$bySet, data.frame(A = 1, B = 1:51), tuneELFs = "a", reportELFs = "b"), "at most 50")
  expect_error(tuneWeights(d$terms, d$bySet, data.frame(A = 1, Z = 1), tuneELFs = "a", reportELFs = "b"), "no column Z")
  expect_error(tuneWeights(d$terms, d$bySet, grid3, tuneELFs = "a", reportELFs = "b", loss = "nope"), "not in the bySet")
})

test_that("a skill loss is maximised", {
  d <- mkWeightsToy()
  out <- tuneWeights(d$terms, d$bySet, grid3, tuneELFs = c("a", "b"), reportELFs = "c", loss = "annualArea_skillN1")
  expect_equal(out$losses$weightId[out$losses$selected], 2L)
})

test_that("the module runs mode 'weights' from the three files and writes the report and the tuning rows", {
  d <- mkWeightsToy(); root <- withr::local_tempdir()
  write.csv(d$terms, file.path(root, "spreadFitTerms_x.csv"), row.names = FALSE)
  write.csv(d$bySet, file.path(root, "wave1_bySet.csv"), row.names = FALSE)
  write.csv(grid3, file.path(root, "grid.csv"), row.names = FALSE)
  out <- runWeights(file.path(root, "spreadFitTerms_x.csv"), file.path(root, "wave1_bySet.csv"), file.path(root, "grid.csv"),
                    outDir = root, tuneELFs = c("a", "b"), reportELFs = "c")
  expect_true(all(file.exists(file.path(root, c("weights_report.csv", "weights_tuning.csv", "weights_losses.csv")))))
  expect_equal(read.csv(file.path(root, "weights_report.csv"))$member, 8L)
})

test_that("module parameter mode = 'weights' reads the files and writes the tables, without a hindcast", {
  skip_if_not_installed("SpaDES.core")
  d <- mkWeightsToy(); root <- withr::local_tempdir(); out <- file.path(root, "outputs"); dir.create(out)
  write.csv(d$terms, file.path(root, "terms.csv"), row.names = FALSE); write.csv(d$bySet, file.path(out, "wave1_bySet.csv"), row.names = FALSE)
  write.csv(grid3, file.path(root, "grid.csv"), row.names = FALSE)
  mp <- file.path(root, "modules"); dir.create(mp)
  file.symlink(normalizePath(testthat::test_path("..", "..")), file.path(mp, "fireSimulator_validation"))
  withr::local_options(spades.useRequire = FALSE, spades.moduleCodeChecks = FALSE)
  SpaDES.core::simInitAndSpades(times = list(start = 0, end = 1), modules = "fireSimulator_validation",
    params = list(fireSimulator_validation = list(mode = "weights", weightsTerms = file.path(root, "terms.csv"), weightsGrid = file.path(root, "grid.csv"),
                                                  weightsTuneELFs = c("a", "b"), weightsReportELFs = "c", .plots = NA)),
    paths = list(modulePath = mp, outputPath = out, inputPath = root, cachePath = file.path(root, "cache"), scratchPath = file.path(root, "scratch")))
  expect_equal(read.csv(file.path(out, "weights_report.csv"))$member, 8L)
  expect_false(file.exists(file.path(out, "wave1.rds")))
})
