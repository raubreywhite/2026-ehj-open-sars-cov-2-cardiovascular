# Illustrates, with 100 persons per group, how each source of bias moves a
# ratio of risks. The numbers are illustrative: they show the mechanism, and are
# not simulation output. In every panel infection multiplies each person's
# risk of a CVD diagnosis by 1.2. Over the follow-up, a high-risk person has a
# 50% risk and a low-risk person a 5% risk. Diagnoses per group are those
# risks summed and rounded, halves up.
#
# It draws figures/illustration_selection.png, _depletion, _unrecorded and
# _deaths.
#
# Needs R 4.6 with data.table, ggplot2 and knitr.
# Run it from the repository root with: Rscript Illustrations.R

library(data.table)
library(ggplot2)

# CFG ====
ILL <- list()
ILL$hr <- 1.2
ILL$risk <- c(high = 0.5, low = 0.05)

# PART 1 -- DATA CREATION ====
# A group: n_high and n_low persons, with or without infection. unrec of the
# persons are infected without a record, split over high and low risk in
# proportion. died high-risk persons die at infection before any diagnosis.
group <- function(panel, label, n_high, n_low, infected, unrec = 0L, died = 0L) {
  d <- data.table(risk = rep(c("high", "low"), c(n_high, n_low)))
  d[, id := .I]
  u_high <- round(unrec * n_high / (n_high + n_low))
  d[, unrec := (risk == "high" & id <= u_high) | (risk == "low" & id > n_high & id <= n_high + unrec - u_high)]
  d[, died := risk == "high" & id > n_high - died]
  d[, p := ILL$risk[risk] * fifelse(infected | unrec, ILL$hr, 1)]
  d[died == TRUE, p := 0]
  # Diagnoses: the expected number in each subgroup of risk and unrecorded
  # infection, rounded, given to the first living persons of that subgroup.
  dx_n <- d[died == FALSE, .(n = floor(sum(p) + 0.5)), by = .(risk, unrec)]
  d[, dx := FALSE]
  for (i in seq_len(nrow(dx_n))) {
    ids <- d[risk == dx_n$risk[i] & unrec == dx_n$unrec[i] & died == FALSE, head(id, dx_n$n[i])]
    d[ids, dx := TRUE]
  }
  d[, `:=`(
    panel = panel,
    label = label,
    row = 10L - (id - 1L) %/% 10L,
    col = (id - 1L) %% 10L + 1L
  )]
  return(d[])
}
P <- c(
  sel = "Selection into early infection",
  dep = "Depletion of susceptibles",
  unr = "Unrecorded infections",
  dai = "Deaths at infection"
)
d <- rbind(
  group(P["sel"], "Infected before Omicron", 10L, 90L, TRUE),
  group(P["sel"], "Comparison group", 20L, 80L, FALSE),
  # Year 2 after a year in which the high-risk persons were diagnosed first:
  # the persons still free of CVD at the start of year 2.
  group(P["dep"], "Infected, CVD-free after a year", 8L, 75L, TRUE),
  group(P["dep"], "Comparison, CVD-free after a year", 10L, 76L, FALSE),
  group(P["unr"], "Infected, recorded", 20L, 80L, TRUE),
  group(P["unr"], "Comparison group", 20L, 80L, FALSE, unrec = 50L),
  group(P["dai"], "Infected", 20L, 80L, TRUE, died = 3L),
  group(P["dai"], "Comparison group", 20L, 80L, FALSE)
)
d[, panel := factor(panel, levels = P)]

# PART 2 -- ANALYSIS ====
# The ratio of risks a study would see: diagnoses per person at risk,
# infected against comparison. Persons who died at infection are not at risk.
tab <- d[, .(n = sum(!died), dx = sum(dx)), keyby = .(panel, label)]
tab[, infected := grepl("^Infected", label)]
ratio <- tab[, .(seen = (dx[infected] / n[infected]) / (dx[!infected] / n[!infected])), keyby = panel]

# PART 3 -- OUTPUT ====
print(knitr::kable(tab, format = "pipe"))
print(knitr::kable(ratio, format = "pipe", digits = 2))
dir.create("results", showWarnings = FALSE)
saveRDS(list(tab = tab, ratio = ratio, cfg = ILL), "results/illustrations.rds")

# One figure per panel: the infected group on the left, the comparison group
# on the right.
d[tab, on = .(panel, label), strip := sprintf("%s\n%d of %d diagnosed", label, i.dx, i.n)]
d[, status := fcase(
  died, "Died at infection",
  risk == "high", "High underlying risk",
  default = "Low underlying risk"
)]
FILL <- c("High underlying risk" = "#4d4d4d", "Low underlying risk" = "#d9d9d9", "Died at infection" = "white")
EDGE_COL <- c("Infected without a record" = "#2b83ba", "Died at infection" = "grey45")
FILES <- c(sel = "selection", dep = "depletion", unr = "unrecorded", dai = "deaths")
for (k in names(P)) {
  x <- d[panel == P[k]]
  x[, strip := factor(strip, levels = unique(x[order(-grepl("^Infected", label)), strip]))]
  x[, edge := fcase(unrec, "Infected without a record", died, "Died at infection", default = NA_character_)]
  q <- ggplot(x, aes(col, row))
  q <- q + geom_tile(aes(fill = status), colour = "white", linewidth = 0.8, width = 0.92, height = 0.92)
  if (any(!is.na(x$edge))) {
    q <- q + geom_tile(data = x[!is.na(edge)], aes(colour = edge), fill = NA, linewidth = 0.9, width = 0.8, height = 0.8)
    q <- q + scale_colour_manual(NULL, values = EDGE_COL[unique(x[!is.na(edge), edge])])
  }
  q <- q + geom_point(data = x[dx == TRUE], aes(shape = "Diagnosed with CVD"), size = 2.6, stroke = 1.2, colour = "#d7191c")
  q <- q + scale_fill_manual(NULL, values = FILL, breaks = setdiff(names(FILL)[names(FILL) %in% x$status], "Died at infection"))
  q <- q + scale_shape_manual(NULL, values = 4)
  q <- q + facet_wrap(~strip, nrow = 1)
  q <- q + coord_equal()
  q <- q +
    theme_void(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 11, margin = margin(b = 4)),
      plot.margin = margin(6, 6, 6, 6)
    )
  ggsave(
    sprintf("figures/illustration_%s.png", FILES[[k]]),
    q,
    width = 195,
    height = 110,
    units = "mm",
    dpi = 180,
    bg = "white"
  )
}
