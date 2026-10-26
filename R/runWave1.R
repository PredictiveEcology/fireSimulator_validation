## The runner: a set of ELFs in, a tidy results table and per-ELF tables out.

#' Wave 1 of the hindcast test plan
#'
#' Gates G1-G6 per ELF; skill tests S1-S5 against the nulls N0 (floor), N1 (the ELF's own history) and N2 (regression on
#' the ELF-median covariate), all fitted in-sample on the same years, scored year included. Verdicts are pooled over ELFs
#' with a bootstrap over ELFs and years (`verdict`: the interval of the skill, or of the Spearman difference to the index, against 0,
#' for N1 and N2 only; N0 is a floor); per-ELF numbers are reported, not judged.
#'
#' Each element of `elfs` is a plain list:
#' \itemize{
#'   \item `elf`: name.
#'   \item `simYears`: `data.frame(year, run, ignitions, escapes, area50)`, one row per simulated run-year (`area50` is
#'     the area, ha, in fires of 50 ha or more).
#'   \item `simSizes`: `data.frame(year, run, ha)`, simulated fire sizes.
#'   \item `obsYears`: `data.frame(year, ignitions, escapes, area50)`, the observed record, every year (scored years included).
#'   \item `obsSizes`: `data.frame(year, ha)`, observed fire sizes, every year.
#'   \item `covariate`: `data.frame(year, value)`, the ELF median of the fit's climate covariate, every year.
#'   \item `fit`: list with `objective` (ledger objective values), `pipelineObjective` (the same objective recomputed by the
#'     pipeline, per replicate; `NULL` if not available), `pc` (`p_c` of each kept set) and `p50`, `fitAreaHa` (area of the
#'     fires the fit used) and `totalAreaHa` (observed area).
#' }
#'
#' @param elfs List of ELF inputs. @param nBoot Bootstrap resamples. @param nNull Draws per year of the AD null.
#' @param nPerm Permutations for the S5 floor. @param seed Seed.
#' @param adFun The Anderson-Darling statistic (default fireSenseUtils').
#' @return A list. `results`: tidy table `(test, null, level, ELF, measure, value, lo, hi, nELF, nUnits, verdict)`, `level`
#'   being `"pooled"` or `"ELF"`; `gates`: `(ELF, gate, pass, value, detail, share100kha)`; `perELF`: named list (S1, S2, S3,
#'   S4ad, S4largest, S5) of per-ELF-year tables; `flagged`: ELFs for which G2 or G3 failed, whose skill tests are not meaningful.
runWave1 <- function(elfs, nBoot = 999, nNull = 199, nPerm = 9999, seed = 1, adFun = fireSenseUtils:::adStatistic) {
  names(elfs) <- vapply(elfs, `[[`, "", "elf")
  gl <- lapply(elfs, runGates)
  gates <- do.call(rbind, Map(function(g, nm) cbind(ELF = nm, g), gl, names(gl)))
  rownames(gates) <- NULL
  flagged <- names(gl)[!vapply(gl, attr, TRUE, "meaningful")]
  tag <- function(l) do.call(rbind, Map(function(x, nm) cbind(ELF = nm, x), l, names(l)))
  s <- lapply(elfs, function(e) list(S1 = scoreS1(e), S2 = scoreS2(e), S3 = scoreS3(e), S4 = scoreS4(e, nNull, adFun = adFun), S5 = scoreS5(e, nPerm)))
  perELF <- list(S1 = tag(lapply(s, function(x) x$S1$scores)), S2 = tag(lapply(s, function(x) x$S2$scores)),
                 S3 = tag(lapply(s, function(x) x$S3$scores)), S4ad = tag(lapply(s, function(x) x$S4$ad$scores)),
                 S4largest = tag(lapply(s, function(x) x$S4$largest$scores)), S5 = tag(lapply(s, `[[`, "S5")))
  cal <- list(S1 = tag(lapply(s, function(x) x$S1$cal)), S2 = tag(lapply(s, function(x) x$S2$cal)),
              S3 = tag(lapply(s, function(x) x$S3$cal)), S4largest = tag(lapply(s, function(x) x$S4$largest$cal)))
  perELF$calibration <- cal
  skillRows <- do.call(rbind, lapply(c("S1", "S2", "S3", "S4ad", "S4largest"), function(tst) {
    d <- perELF[[tst]]
    do.call(rbind, lapply(unique(d$null), function(nl) {
      x <- d[d$null == nl, ]
      perE <- do.call(rbind, lapply(split(x, x$ELF), function(z) data.frame(level = "ELF", ELF = z$ELF[1], measure = "skill",
                value = skillScore(z$lossModel, z$lossNull), lo = NA_real_, hi = NA_real_, nELF = 1L, nUnits = nrow(z))))
      p <- bootstrapSkill(x, nBoot, seed)
      rbind(data.frame(level = "pooled", ELF = "all", measure = "skill", value = p$skill, lo = p$lo, hi = p$hi, nELF = p$nELF, nUnits = p$nUnits), perE) |>
        transform(test = tst, null = nl)
    }))
  }))
  calRows <- do.call(rbind, lapply(names(cal), function(tst) {
    d <- cal[[tst]]
    do.call(rbind, lapply(intersect(c("cover50", "cover90", "logScore"), names(d)), function(m) {
      x <- data.frame(ELF = d$ELF, year = d$year, value = as.numeric(d[[m]]))
      p <- bootstrapMean(x, nBoot, seed)
      perE <- do.call(rbind, lapply(split(x, x$ELF), function(z) data.frame(level = "ELF", ELF = z$ELF[1], value = pooledValue(z$value),
                lo = NA_real_, hi = NA_real_, nELF = 1L, nUnits = nrow(z))))
      rbind(data.frame(level = "pooled", ELF = "all", p), perE) |> transform(test = tst, null = NA_character_, measure = m)
    }))
  }))
  rank <- tag(lapply(s, `[[`, "S5"))
  rankD <- do.call(rbind, lapply(names(s), function(nm) cbind(ELF = nm, attr(s[[nm]]$S5, "years"))))
  rb <- bootRankRows(rankD, nBoot, seed)
  rankRows <- rbind(rb, do.call(rbind, lapply(seq_len(nrow(rank)), function(i) data.frame(level = "ELF", ELF = rank$ELF[i],
              measure = c("rhoModel", "rhoIndex", "rhoDiff", "pPermN0"), value = unlist(rank[i, c("rho", "rhoIndex", "rhoDiff", "pPerm")]),
              lo = NA_real_, hi = NA_real_, nELF = 1L, nUnits = rank$nYears[i]))))
  rankRows$test <- "S5"; rankRows$null <- ifelse(rankRows$measure == "pPermN0", "N0", "N2")
  res <- rbind(skillRows, calRows, rankRows)
  judged <- res$level == "pooled" & res$measure %in% c("skill", "rhoDiff") & res$null %in% c("N1", "N2")   # N0 is a floor, not judged
  res$verdict <- ifelse(judged & is.finite(res$lo),
                        ifelse(res$lo > 0, "beats null", ifelse(res$hi < 0, "worse than null", "no difference")), NA_character_)
  res <- res[, c("test", "null", "level", "ELF", "measure", "value", "lo", "hi", "nELF", "nUnits", "verdict")]
  rownames(res) <- NULL
  list(results = res, gates = gates, perELF = perELF, flagged = flagged)
}

bootRankRows <- function(rankD, nBoot, seed) {
  b <- bootstrapRankDiff(rankD, nBoot, seed)
  data.frame(level = "pooled", ELF = "all", measure = b$measure, value = b$value, lo = b$lo, hi = b$hi,
             nELF = length(unique(rankD$ELF)), nUnits = nrow(rankD))
}
