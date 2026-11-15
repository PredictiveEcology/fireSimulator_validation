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

## a terms file as fireSense_spreadFit writes it: member, rep, parameters, terms; no ELF column; the ELF is in the run name
mkTermsFile <- function(dir, elf, d, runName = paste0("fit_suffix_", elf, "_NA_NA_1")) {
  x <- d$terms[d$terms$ELF == elf, ]
  x <- do.call(rbind, lapply(1:2, function(r) data.frame(member = x$member, rep = r, youngAge = 0.1, A = x$A + c(0, 100)[r], B = x$B)))   # rep 2 of A is far off
  f <- file.path(dir, paste0("spreadFitTerms_", runName, ".csv")); write.csv(x, f, row.names = FALSE); f
}

test_that("terms files of several ELFs: the ELF comes from the run name, or from the names of the files", {
  d <- mkWeightsToy(c("5.2", "11.3")); root <- withr::local_tempdir()
  fs <- c(mkTermsFile(root, "5.2", d), mkTermsFile(root, "11.3", d))
  t <- readTermsFiles(fs)
  expect_setequal(unique(t$ELF), c("5.2", "11.3")); expect_true(all(c("member", "rep", "youngAge") %in% names(t)))
  expect_equal(unique(readTermsFiles(stats::setNames(fs[1], "x"))$ELF), "x")
  bad <- file.path(root, "spreadFitTerms_nothing.csv"); file.copy(fs[1], bad)
  expect_error(readTermsFiles(bad), "cannot tell the ELF")
})

test_that("replicates are reduced to their median before selection (the minimum over replicates would pick a lucky one)", {
  d <- mkWeightsToy(c("5.2", "11.3")); root <- withr::local_tempdir()
  t <- readTermsFiles(c(mkTermsFile(root, "5.2", d), mkTermsFile(root, "11.3", d)))
  g <- data.frame(A = 1, B = 0)
  # median of A over the 2 reps is k + 50, still increasing in k: member 1; by the single lucky rep 1 also 1 - so make rep 1 lucky for member 10
  t$A[t$rep == 1 & t$member == 10] <- -5
  out <- tuneWeights(t, d$bySet[d$bySet$ELF %in% c("5.2", "11.3"), ], g, tuneELFs = "5.2", reportELFs = "11.3")
  expect_equal(out$report$member, 1L)                                   # median of (-5, 110) = 52.5 > 51
  expect_equal(nrow(out$report), 1L)
})

test_that("a terms member with no bySet row stops, naming the count", {
  d <- mkWeightsToy(); d$bySet <- d$bySet[d$bySet$paramSet <= 8, ]
  expect_error(tuneWeights(d$terms, d$bySet, grid3, tuneELFs = c("a", "b"), reportELFs = "c"), "2 of 10 terms members have no bySet row")
})

test_that("NA terms stop; a NA loss is left out of the median with a warning naming ELF and weightId; all-NA stops", {
  d <- mkWeightsToy()
  t <- d$terms; t$A[2] <- NA
  expect_error(tuneWeights(t, d$bySet, grid3, tuneELFs = c("a", "b"), reportELFs = "c"), "NA or not finite: ELF a, member 2")
  t <- d$terms; t$B[3] <- Inf
  expect_error(tuneWeights(t, d$bySet, grid3, tuneELFs = c("a", "b"), reportELFs = "c"), "not finite")
  b <- d$bySet; b$area_ratio[b$ELF == "a" & b$paramSet == 8] <- NA          # weight (0, 1) picks member 8 on ELF a
  expect_warning(out <- tuneWeights(d$terms, b, grid3, tuneELFs = c("a", "b"), reportELFs = "c"), "ELF a, weightId 2")
  expect_equal(out$losses$weightId[out$losses$selected], 2L)               # ELF b still scores it perfectly
  b$area_ratio <- NA_real_
  expect_error(suppressWarnings(tuneWeights(d$terms, b, grid3, tuneELFs = c("a", "b"), reportELFs = "c")), "no finite loss")
})

test_that("ties go to the lowest member whatever the file order", {
  d <- mkWeightsToy(); g <- data.frame(A = 0, B = 0)
  out <- tuneWeights(d$terms[nrow(d$terms):1, ], d$bySet, g, tuneELFs = c("a", "b"), reportELFs = "c")
  expect_equal(out$report$member, 1L)
})
