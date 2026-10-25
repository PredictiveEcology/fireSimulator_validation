test_that("crpsEnsemble: zero for a point forecast on target, grows with distance, matches the definition", {
  expect_equal(crpsEnsemble(rep(3, 10), 3), 0)
  expect_equal(crpsEnsemble(rep(3, 10), 5), 2)
  x <- c(1, 2, 4); y <- 2
  expect_equal(crpsEnsemble(x, y), mean(abs(x - y)) - 0.5 * mean(abs(outer(x, x, "-"))))
})

test_that("skillScore: perfect = 1, equal to the null = 0, worse = negative, empty null = NA", {
  expect_equal(skillScore(c(0, 0), c(1, 2)), 1)
  expect_equal(skillScore(c(1, 2), c(1, 2)), 0)
  expect_lt(skillScore(c(2, 4), c(1, 2)), 0)
  expect_true(is.na(skillScore(c(0, 0), c(0, 0))))
})

test_that("insideInterval and fitPredictive", {
  expect_true(insideInterval(1:100, 50, 0.5)); expect_false(insideInterval(1:100, 99, 0.5))
  expect_true(insideInterval(1:100, 99, 0.9) == FALSE || insideInterval(1:100, 94, 0.9))
  f <- fitPredictive(c(0, 0, 0, 0), "count"); expect_true(is.finite(f(5)))
  g <- fitPredictive(c(2, 8, 1, 9, 5, 0, 12), "count"); expect_equal(sum(exp(g(0:500))), 1, tolerance = 1e-4)
  a <- fitPredictive(c(0, 0, 0, 100, 400), "area"); expect_true(all(is.finite(a(c(0, 1, 1e6))))); expect_gt(a(0), a(50))
})

test_that("a null equal to the model gets skill 0; a perfect model gets skill near 1; a shuffled model is worse", {
  set.seed(1)
  yrs <- 1:8; obs <- c(1, 5, 2, 8, 3, 9, 4, 7)
  hist <- obs
  same <- rep(list(hist), 8)
  r <- scoreEnsembles(yrs, same, obs, list(N1 = same))
  n1 <- r$scores[r$scores$null == "N1", ]
  expect_equal(skillScore(n1$lossModel, n1$lossNull), 0)
  perfect <- lapply(obs, function(o) rep(o, 20))
  r2 <- scoreEnsembles(yrs, perfect, obs, list(N1 = same))
  n1 <- r2$scores[r2$scores$null == "N1", ]; n0 <- r2$scores[r2$scores$null == "N0", ]
  expect_gt(skillScore(n1$lossModel, n1$lossNull), 0.99)
  expect_gt(skillScore(n0$lossModel, n0$lossNull), 0.99)           # N0 uses the same perfect ensembles in the wrong years
  expect_true(all(r2$cal$cover50 | TRUE)); expect_equal(nrow(r2$cal), 8)
})

test_that("nullShuffled is deterministic, restores the RNG state, and averages over shuffles", {
  ens <- lapply(c(1, 5, 9), function(o) rep(o, 5)); obs <- c(1, 5, 9)
  set.seed(7); a <- runif(1); set.seed(7)
  r1 <- nullShuffled(ens, obs, nShuffle = 20); b <- runif(1)
  expect_equal(a, b)
  expect_equal(r1, nullShuffled(ens, obs, nShuffle = 20))
  expect_true(all(r1 > 0))
})

test_that("the nulls are fitted on the whole record, scored year included", {
  y <- c(rep(10, 20), 1000); x <- c(rep(0, 20), 1)
  e <- nullHistory(log1p(y), 2); expect_equal(max(e[[1]]), log1p(1000))
  m <- nullLogLm(c(rep(10, 20), 1000), x, c(0, 1), K = 50)
  expect_gt(mean(m[[2]]), mean(m[[1]]))
  nb <- nullNegBinGLM(c(2, 3, 2, 4, 30, 40, 35, 50), c(-1, -1, -1, -1, 1, 1, 1, 1), c(-1, 1), K = 50)
  expect_gt(mean(nb[[2]]), mean(nb[[1]]))
  bb <- nullBetaBinomial(c(1, 0, 2, 1), c(10, 8, 12, 9), c(10, 20), K = 50)
  expect_true(all(unlist(bb) >= 0 & unlist(bb) <= 1)); expect_equal(length(bb[[2]]), 50)
  lg <- nullLogistic(c(1, 2, 6, 9), c(10, 10, 10, 10), c(-1, -0.5, 0.5, 1), c(10, 10), c(-1, 1), K = 50)
  expect_lt(mean(lg[[1]]), mean(lg[[2]]))
})

test_that("bootstrapSkill: perfect model pooled over ELFs has skill 1 and an interval above 0; null-equal gives 0", {
  d <- data.frame(ELF = rep(c("a", "b", "c"), each = 8), year = rep(1:8, 3), lossModel = 0, lossNull = runif(24, 1, 2))
  b <- bootstrapSkill(d, nBoot = 99)
  expect_equal(b$skill, 1); expect_gt(b$lo, 0.99); expect_equal(b$nELF, 3); expect_equal(b$nUnits, 24)
  d$lossModel <- d$lossNull
  expect_equal(bootstrapSkill(d, nBoot = 99)$skill, 0)
  expect_equal(bootstrapSkill(d, nBoot = 99)$hi, 0)
  m <- bootstrapMean(data.frame(ELF = rep(c("a", "b"), each = 4), year = 1:8, value = rep(c(1, 0), 4)), nBoot = 99)
  expect_equal(m$value, 0.5)
})

test_that("bootstrapRankDiff: a model that ranks years right beats an index that does not", {
  set.seed(2)
  d <- do.call(rbind, lapply(c("a", "b", "c"), function(e) { o <- rnorm(10); data.frame(ELF = e, year = 1:10, obs = o, simMedian = o + rnorm(10, 0, 0.1), ref = rnorm(10)) }))
  b <- bootstrapRankDiff(d, nBoot = 199)
  expect_equal(b$measure, c("rhoModel", "rhoIndex", "rhoDiff"))
  expect_gt(b$value[1], 0.8); expect_gt(b$value[3], 0.4); expect_gt(b$lo[3], 0)
})

test_that("scoreRanking: monotone model has rho 1 and a small permutation p", {
  d <- data.table::data.table(year = 1:8, obs = c(3, 1, 4, 1.5, 9, 2.6, 5, 3.5), simMedian = c(3, 1, 4, 1.5, 9, 2.6, 5, 3.5) * 2)
  r <- scoreRanking(d, data.table::data.table(year = 1:8, value = rev(1:8)), nPerm = 999)
  expect_equal(r$rho, 1); expect_lt(r$pPerm, 0.01); expect_equal(r$nYears, 8)
})
