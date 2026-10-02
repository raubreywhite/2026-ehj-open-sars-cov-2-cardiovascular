# Illustrates, with 100 persons per group, how each source of bias moves a
# ratio of risks. The numbers are illustrative: they show the mechanism, and
# are not simulation output. In every example infection multiplies each
# person's risk of a CVD diagnosis by 1.2. Over the follow-up (each year, in
# the depletion example), a high-risk person has a 50% risk and a low-risk
# person a 5% risk. Diagnoses per subgroup are those risks summed and rounded,
# halves up. Each figure shows the persons as squares, and the working under
# them.
#
# It draws figures/illustration_selection.png, _depletion, _unrecorded and
# _deaths.
#
# Needs R 4.6 with data.table, ggplot2, patchwork and knitr.
# Run it from the repository root with: Rscript Illustrations.R

library(data.table)
library(ggplot2)

# CFG ====
ILL <- list()
ILL$hr <- 1.2
ILL$risk <- c(high = 0.5, low = 0.05)
P <- c(
  sel = "Selection into early infection",
  dep = "Depletion of susceptibles",
  unr = "Unrecorded infections",
  dai = "Deaths at infection"
)

# PART 1 -- DATA CREATION ====
# One group at one stage: n_high and n_low persons, with or without
# infection. unrec of them are infected without a record, split over high and
# low risk in proportion. died high-risk persons die at infection before any
# diagnosis. Diagnoses: the expected number per subgroup of risk and
# unrecorded infection, rounded halves up, given to the first living persons
# of the subgroup.
group <- function(
  panel,
  stage,
  label,
  infected,
  n_high,
  n_low,
  unrec = 0L,
  died = 0L
) {
  d <- data.table(risk = rep(c("high", "low"), c(n_high, n_low)))
  d[, id := .I]
  u_high <- round(unrec * n_high / (n_high + n_low))
  d[,
    unrec := (risk == "high" & id <= u_high) |
      (risk == "low" & id > n_high & id <= n_high + unrec - u_high)
  ]
  d[, died := risk == "high" & id > n_high - died]
  d[, infected := infected | unrec]
  d[, p := ILL$risk[risk] * fifelse(infected, ILL$hr, 1)]
  d[died == TRUE, p := 0]
  dx_n <- d[died == FALSE, .(n = floor(sum(p) + 0.5)), by = .(risk, unrec)]
  d[, dx := FALSE]
  for (i in seq_len(nrow(dx_n))) {
    ids <- d[
      risk == dx_n$risk[i] & unrec == dx_n$unrec[i] & died == FALSE,
      head(id, dx_n$n[i])
    ]
    d[ids, dx := TRUE]
  }
  d[, removed := FALSE]
  d[, `:=`(
    panel = panel,
    stage = stage,
    label = label,
    row = 10L - (id - 1L) %/% 10L,
    col = (id - 1L) %% 10L + 1L
  )]
  return(d[])
}
# Depletion: year 1 with all 100 persons, then year 2 in the same positions.
# The persons diagnosed in year 1 are removed: they are no longer free of
# CVD. The year-2 diagnoses go to the first persons still followed in each
# risk group.
depletion <- function(label, infected) {
  y1 <- group(P["dep"], "Year 1", label, infected, 20L, 80L)
  y2 <- copy(y1)[, `:=`(stage = "Year 2", removed = dx, dx = FALSE)]
  for (r in c("high", "low")) {
    ids <- y2[risk == r & removed == FALSE, id]
    y2[head(ids, floor(sum(y2[ids, p]) + 0.5)), dx := TRUE]
  }
  return(rbind(y1, y2))
}
d <- rbind(
  group(P["sel"], "", "Infected before Omicron", TRUE, 10L, 90L),
  group(P["sel"], "", "Comparison group", FALSE, 20L, 80L),
  depletion("Infected", TRUE),
  depletion("Comparison group", FALSE),
  group(P["unr"], "", "Infected, recorded", TRUE, 20L, 80L),
  group(P["unr"], "", "Comparison group", FALSE, 20L, 80L, unrec = 50L),
  group(P["dai"], "", "Infected", TRUE, 20L, 80L, died = 3L),
  group(P["dai"], "", "Comparison group", FALSE, 20L, 80L)
)
d[, exposed := label %like% "^Infected"]

# PART 2 -- ANALYSIS ====
# Per group and stage: persons at risk (those who did not die at infection
# and were not removed after a diagnosis in year 1) and diagnoses. The ratio
# of risks, infected against comparison, per stage.
# ex: the expected diagnoses, unrounded. The squares show the rounded
# diagnoses dx; the working and the ratio use ex. In the depletion example
# the persons removed after year 1 are whole persons, so rounding still moves
# the year-2 ratio slightly.
tab <- d[,
  .(n = sum(!died & !removed), dx = sum(dx), ex = sum(p[!died & !removed]), died = sum(died)),
  keyby = .(panel, stage, label, exposed)
]
ratio <- tab[,
  .(seen = (ex[exposed] / n[exposed]) / (ex[!exposed] / n[!exposed])),
  keyby = .(panel, stage)
]
# The working per subgroup of risk and unrecorded infection.
cells <- d[
  died == FALSE & removed == FALSE,
  .(n = .N, p = p[1], expected = sum(p), dx = sum(dx), infected = infected[1]),
  keyby = .(panel, stage, label, exposed, risk, unrec)
]
pc <- function(x) sprintf("%.0f%%", 100 * x)
# A count without decimals when it is whole, else with one.
num <- function(x) if (abs(x - round(x)) < 1e-9) sprintf("%d", round(x)) else sprintf("%.1f", x)
term <- function(n, r, infected, expected, dx) {
  return(sprintf(
    "%d × %s%s = %s",
    n,
    pc(ILL$risk[[r]]),
    if (infected) sprintf(" × %.1f", ILL$hr) else "",
    num(expected)
  ))
}
# The lines of working under each figure: one per group and stage, then the
# ratio of each stage.
working <- function(k) {
  out <- character(0)
  for (sg in unique(tab[panel == P[k], stage])) {
    for (g in tab[panel == P[k] & stage == sg][order(-exposed), label]) {
      cl <- cells[panel == P[k] & stage == sg & label == g]
      parts <- vapply(
        seq_len(nrow(cl)),
        function(i) {
          sprintf(
            "%s %s%s",
            term(cl$n[i], cl$risk[i], cl$infected[i], cl$expected[i], cl$dx[i]),
            if (cl$risk[i] == "high") "high-risk" else "low-risk",
            if (cl$unrec[i]) " without a record" else ""
          )
        },
        character(1)
      )
      tb <- tab[panel == P[k] & stage == sg & label == g]
      # More than two terms go on a second, indented line.
      first <- head(parts, 2L)
      rest <- tail(parts, -2L)
      tail_txt <- sprintf(
        "   →   %s / %d = %.1f%%%s",
        num(tb$ex),
        tb$n,
        100 * tb$ex / tb$n,
        if (tb$died > 0) sprintf("   (%d died at infection)", tb$died) else ""
      )
      head_txt <- sprintf(
        "%s%s:   %s",
        if (sg == "") "" else paste0(sg, ", "),
        g,
        paste(first, collapse = "  +  ")
      )
      if (length(rest) == 0L) {
        out <- c(out, paste0(head_txt, tail_txt))
      } else {
        out <- c(
          out,
          head_txt,
          paste0("        +  ", paste(rest, collapse = "  +  "), tail_txt)
        )
      }
    }
    ti <- tab[panel == P[k] & stage == sg & exposed == TRUE]
    tc <- tab[panel == P[k] & stage == sg & exposed == FALSE]
    out <- c(
      out,
      sprintf(
        "%s:   %.1f%% ÷ %.1f%% = %.2f   (true %.2f)",
        if (sg == "") {
          "The study would see"
        } else {
          paste0(sg, ", the study would see")
        },
        100 * ti$ex / ti$n,
        100 * tc$ex / tc$n,
        ratio[panel == P[k] & stage == sg, seen],
        ILL$hr
      ),
      ""
    )
  }
  return(c(head(out, -1L), "", "The squares show whole persons; the working uses the unrounded numbers."))
}

# PART 3 -- OUTPUT ====
print(knitr::kable(tab, format = "pipe"))
print(knitr::kable(ratio, format = "pipe", digits = 2))
dir.create("results", showWarnings = FALSE)
saveRDS(
  list(tab = tab, ratio = ratio, cells = cells, cfg = ILL),
  "results/illustrations.rds"
)

FILL <- c(
  "High underlying risk" = "#4d4d4d",
  "Low underlying risk" = "#d9d9d9",
  "Died at infection" = "white",
  "Removed" = "white"
)
EDGE_COL <- c(
  "Infected without a record" = "#2b83ba",
  "Died at infection" = "grey45",
  "Diagnosed with CVD in year 1, no longer followed" = "#f4a6a6"
)
FILES <- c(
  sel = "selection",
  dep = "depletion",
  unr = "unrecorded",
  dai = "deaths"
)
d[
  tab,
  on = .(panel, stage, label),
  strip := sprintf(
    "%s%s\n%d of %d diagnosed",
    fifelse(stage == "", "", paste0(stage, ": ")),
    label,
    i.dx,
    i.n
  )
]
d[,
  status := fcase(
    died , "Died at infection" , removed , "Removed" , risk == "high" , "High underlying risk" ,
    default = "Low underlying risk"
  )
]
d[,
  edge := fcase(
    unrec , "Infected without a record" , died , "Died at infection" , removed , "Diagnosed with CVD in year 1, no longer followed" ,
    default = NA_character_
  )
]
for (k in names(P)) {
  x <- d[panel == P[k]]
  x[, strip := factor(strip, levels = unique(x[order(stage, -exposed), strip]))]
  q <- ggplot(x, aes(col, row))
  q <- q +
    geom_tile(
      aes(fill = status),
      colour = "white",
      linewidth = 0.8,
      width = 0.92,
      height = 0.92
    )
  if (any(!is.na(x$edge))) {
    q <- q +
      geom_tile(
        data = x[!is.na(edge)],
        aes(colour = edge),
        fill = NA,
        linewidth = 0.9,
        width = 0.8,
        height = 0.8
      )
    q <- q +
      scale_colour_manual(
        NULL,
        values = EDGE_COL[unique(x[!is.na(edge), edge])]
      )
  }
  q <- q +
    geom_point(
      data = x[dx == TRUE],
      aes(shape = "Diagnosed with CVD"),
      size = 2.4,
      stroke = 1.1,
      colour = "#d7191c"
    )
  q <- q +
    scale_fill_manual(
      NULL,
      values = FILL,
      breaks = setdiff(
        names(FILL)[names(FILL) %in% x$status],
        c("Died at infection", "Removed")
      )
    )
  q <- q + scale_shape_manual(NULL, values = 4)
  q <- q + scale_y_continuous(limits = c(0.5, 10.5))
  q <- q + facet_wrap(~strip, ncol = 2)
  q <- q + coord_equal()
  q <- q +
    theme_void(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 11, margin = margin(b = 4)),
      plot.margin = margin(6, 6, 6, 6)
    )
  # The working, as lines of text in a white space under the squares.
  w <- working(k)
  tx <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = rev(seq_along(w)),
      label = w,
      hjust = 0,
      size = 3.4
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, length(w) + 0.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(6, 12, 10, 12))
  n_rows <- length(unique(x$stage))
  qq <- patchwork::wrap_plots(
    q,
    tx,
    ncol = 1,
    heights = c(n_rows * 105, 7 * length(w))
  )
  ggsave(
    sprintf("figures/illustration_%s.png", FILES[[k]]),
    qq,
    width = 240,
    height = 105 * n_rows + 30 + 7 * length(w),
    units = "mm",
    dpi = 170,
    bg = "white"
  )
}
