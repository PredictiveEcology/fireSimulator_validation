## One self-contained HTML page of the validation: scores table, observed rank, overview, a figure per metric and per ELF.

## where the observed record ranks among the single runs (mid-rank percentile, 0-1), for the five ratios of the gates
observedRanks <- function(elf, minHa = 50) {
  sy <- elf$simYears; ob <- elf$obsYears; os <- elf$obsSizes[elf$obsSizes$ha >= minHa, ]
  rank <- function(v, o) mean(v < o) + 0.5 * mean(v == o)
  top10 <- function(x) if (length(x)) stats::median(utils::head(sort(x, TRUE), 10L)) else 0
  run <- stats::aggregate(cbind(ignitions, escapes, area50) ~ run, sy, sum)
  sz <- elf$simSizes[elf$simSizes$ha >= minHa, ]
  t10 <- vapply(run$run, function(r) top10(sz$ha[sz$run == r]), 0)
  c(area = rank(run$area50, sum(ob$area50)), top10 = rank(t10, top10(os$ha)), ign = rank(run$ignitions, sum(ob$ignitions)),
    esc = rank(run$escapes / run$ignitions, sum(ob$escapes) / sum(ob$ignitions)),
    ape = rank(run$area50 / pmax(run$escapes, 1), sum(ob$area50) / sum(ob$escapes)))
}

#' The validation page
#'
#' Draws the figures of [plotELFPanel()], [plotOverview()] and [plotMetricPanel()] and writes one HTML file with the images
#' inside it (no other files to publish): a sortable table of the five gate ratios with the spread of single runs and the
#' observed rank, a chart of the observed ranks, the overview, one figure per metric and one per ELF. The layout is
#' `pageTemplate.html` of the module.
#'
#' @param elfs List of ELF inputs ([runWave1()]); `elf` is the id shown on the page.
#' @param gates The `gates` table of [runWave1()] (its `ELF` values are the ids of `elfs`). Only G2, G4 and G6 are used.
#' @param file Where to write the page. @param template Path of `pageTemplate.html`.
#' @param subtitle One line under the title. @param label Name of the simulated series in the figures. @param minHa Escape size.
#' @param dpi Resolution of the inlined figures (80 keeps a 16-ELF page near 9 MB).
#' @return `file`, invisibly.
buildValidationPage <- function(elfs, gates, file, template, subtitle = "", label = "simulated", minHa = 50, dpi = 80) {
  ids <- vapply(elfs, `[[`, "", "elf"); names(elfs) <- ids
  summ <- summariseGates(gates); summ$ELF <- as.character(summ$ELF)
  pl <- lapply(elfs, elfPlots, label = label, minHa = minHa)
  tmp <- tempfile("page"); dir.create(tmp)
  save <- function(p, name, size) {
    size$dpi <- dpi; f <- file.path(tmp, paste0(name, ".png")); do.call(ggplot2::ggsave, c(list(f, p), size))
    paste0("data:image/png;base64,", gsub("\n", "", jsonlite::base64_enc(readBin(f, "raw", file.size(f))), fixed = TRUE))
  }
  img <- c(stats::setNames(lapply(ids, function(id) save(plotELFPanel(elfs[[id]], summ, label, minHa, plots = pl[[id]]), paste0("panel_", id), panelSize)), paste0("panel_", ids)),
           stats::setNames(lapply(names(metrics), function(m) save(plotMetricPanel(pl, summ, m), paste0("metric_", m), metricSize(length(pl)))), paste0("metric_", names(metrics))))
  overview <- save(plotOverview(summ), "overview", overviewSize)
  key <- c(area = "area", top10 = "top10_median", ign = "ignitions", esc = "escapeRate", ape = "areaPerEscape")
  g <- gates[gates$gate %in% c("G2", "G4", "G6"), ]
  rows <- lapply(ids, function(id) {
    z <- summ[summ$ELF == id, ]; gi <- g[g$ELF == id, ]; rk <- observedRanks(elfs[[id]], minHa)
    r <- list(elf = id, gates = sum(gi$pass %in% TRUE), failed = paste(gi$gate[gi$pass %in% FALSE], collapse = " "))
    for (k in names(key)) r[[k]] <- round(c(z[[paste0(key[k], "_ratio")]], z[[paste0(key[k], "_runLo")]], z[[paste0(key[k], "_runHi")]], rk[[k]]), 3)
    r })
  tpl <- paste(readLines(template, warn = FALSE), collapse = "\n")
  fill <- function(pat, val) tpl <<- sub(pat, val, tpl, fixed = TRUE)
  fill("/*ROWS*/[]", as.character(jsonlite::toJSON(rows, auto_unbox = TRUE)))
  fill("/*IMG*/{}", as.character(jsonlite::toJSON(img, auto_unbox = TRUE)))
  tpl <- gsub("__OVERVIEW__", overview, tpl, fixed = TRUE)
  tpl <- gsub("__SUBTITLE__", subtitle, tpl, fixed = TRUE); tpl <- gsub("__NGATES__", "3", tpl, fixed = TRUE)
  tpl <- gsub("__UPDATED__", format(Sys.time(), "%b %d, %H:%M"), tpl, fixed = TRUE)
  writeLines(tpl, file)
  invisible(file)
}
