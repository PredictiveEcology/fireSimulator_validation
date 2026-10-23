## The new code against numbers the hindcast-metrics prototype produced (needs the cluster's /mnt/fast and the session folder).
P <- "/home/emcintir/claudeSessions/2026-10-06-multiELF-predict/phase4"
H <- "/mnt/fast/outputs/hindcast_boundfix_rep81"

test_that("crpsEnsemble on log(1 + x) reproduces the prototype's per ELF-year CRPS of ignitions and fires of 50 ha or more", {
  skip_if_not(file.exists(file.path(P, "eval_boundfix/new/perELFyear.csv")) && dir.exists(H))
  old <- data.table::fread(file.path(P, "eval_boundfix/new/perELFyear.csv"))
  yrs <- seq(1986L, 2021L, 5L); e <- "6.6.1"
  s <- readHindcastSim(H, e, hindcastFinishedYears(H, e, yrs), 1:5)$simYears
  o <- do.call(rbind, lapply(yrs, function(y) as.data.frame(data.table::fread(file.path(H, e, y, 1, "hindcast_observed.csv")))))
  for (q in c("ignitions", "nFires50")) {
    col <- if (q == "nFires50") "escapes" else q
    new <- vapply(yrs, function(y) crpsEnsemble(log1p(s[[col]][s$year == y]), log1p(o[[q]][o$year == y])), 0)
    expect_equal(new, old$crps[old$ELF == e & old$quantity == q][order(old$hcYear[old$ELF == e & old$quantity == q])], tolerance = 1e-9)
  }
})

test_that("scoreSizesVsNull reproduces the prototype's pooled Anderson-Darling of the 6.6.1 power-law hindcast", {
  skip_if_not(file.exists(file.path(P, "metrics/powerLaw661/out/all_sizesPooled.csv")) && dir.exists("/mnt/fast/outputs/hindcast_powerLaw_rep87"))
  skip_if_not_installed("kSamples")
  HP <- "/mnt/fast/outputs/hindcast_powerLaw_rep87"; yrs <- seq(1986L, 2021L, 5L); e <- "6.6.1"
  pp <- data.table::fread(file.path(P, "metrics/powerLaw661/out/all_sizesPooled.csv"))
  pp <- pp[grepl("^Power law", pp$model) & pp$null == "history" & pp$ELF == e, ]
  s <- readHindcastSim(HP, e, hindcastFinishedYears(HP, e, yrs), 1:5)$simSizes
  ob <- do.call(rbind, lapply(yrs, function(y) { x <- data.table::fread(file.path(HP, e, y, 1, "hindcast_observedSizes.csv")); h <- x$ha[x$cause == "natural"]; data.frame(year = rep(y, length(h)), ha = h) }))
  r <- scoreSizesVsNull(s, ob, ob, nNull = 19, nulls = "history")$pooled
  expect_equal(r$adModel, pp$adModel, tolerance = 1e-6)
})
