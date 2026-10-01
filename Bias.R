# Which sources of bias move the estimates at 12 months or more below the true
# hazard ratios, and how much selection into early infection is needed? It
# uses the settings of Run.R, so run Run.R first.
#
# Selection: the model of Run.R with pre_frailty 0, 0.025, 0.05, 0.075, 0.10
# and 0.15 (the value of Run.R). Per value, the geometric mean of the 12
# estimates at 12 months or more, analysed as the paper does, and the ratio of
# mean frailty of persons with a recorded infection before the Omicron wave to
# that of persons infected later.
#
# Decomposition: a Shapley decomposition of the bias, log(true geometric mean)
# - log(estimated geometric mean) of the 12 outcomes, per window from day 2 to
# 12 months or more, over four sources, each switched on or off. Day 0-1 is
# left out: some outcomes have fewer than 5 events there in some cohorts.
#   risk:  unmeasured differences in risk with selection into early infection
#          (off: frailty SD 1e-6 and pre_frailty 0; an SD of 0 divides by 0);
#   unrec: unrecorded infections (off: contam 1e-6, about 2 persons; a share
#          of 0 draws no dates and stops sample());
#   dai:   death at infection (off: p80 0);
#   cens:  censoring at the first CVD diagnosis of any type (off: each outcome
#          followed to its own first diagnosis, death or fu).
# A source's share is its mean marginal effect over all orders in which the
# sources can be switched on. The shares add up to the bias with all sources
# on less the bias with all off.
#
# Each cohort is saved to results/bias/ as it finishes, with the md5 of Run.R
# and R/functions.R, so a stopped run resumes. A saved cohort from other code
# stops the script. With all sources on and censoring, the estimates MUST
# equal those of Run.R.
#
# It draws figures/bias_steps.png: per window, the estimate as the sources
# are added one at a time.
#
# Needs R 4.6 with data.table, ggplot2 and knitr.
# Run it from the repository root with: Rscript Bias.R

library(data.table)
library(ggplot2)

# CFG ====
# The settings, functions and SSI checks of Run.R: every top-level expression
# before its simulation.
ex <- parse("Run.R")
stop_at <- which(vapply(
  ex,
  function(e) startsWith(deparse(e)[1], "res <- parallel::mclapply"),
  logical(1)
))[1]
for (e in ex[seq_len(stop_at - 1L)]) {
  eval(e)
}
BIAS <- list()
BIAS$pf <- c(0, 0.025, 0.05, 0.075, 0.10, CFG$pre_frailty)
BIAS$src <- c("risk", "unrec", "dai", "cens")
BIAS$dir <- "results/bias"
BIAS$windows <- 2:5
BIAS$code <- unname(tools::md5sum(c("Run.R", "R/functions.R")))
cfg_of <- function(frail, unrec, pf, dai) {
  return(modifyList(
    CFG,
    list(
      fsd = if (frail) CFG$fsd else 1e-6,
      contam = if (unrec) CFG$contam else 1e-6,
      pre_frailty = pf,
      p80 = if (dai) CFG$p80 else 0
    )
  ))
}

# PART 1 -- DATA CREATION ====
# The settings: the selection values with all other sources on, and the 8
# combinations of risk, unrec and dai. Risk off means no frailty and no
# selection. Two settings appear in both lists once only.
set_sel <- data.table(frail = 1L, unrec = 1L, pf = BIAS$pf, dai = 1L)
set_shp <- CJ(risk = 0:1, unrec = 0:1, dai = 0:1)[, .(
  frail = risk,
  unrec,
  pf = risk * CFG$pre_frailty,
  dai
)]
grid <- unique(rbind(set_sel, set_shp))
grid <- grid[, .(seed = CFG$seeds), by = names(grid)]
grid[,
  file := sprintf(
    "%s/f%d_u%d_pf%.3f_d%d_seed%02d.rds",
    BIAS$dir,
    frail,
    unrec,
    pf,
    dai,
    seed
  )
]
dir.create(BIAS$dir, showWarnings = FALSE, recursive = TRUE)
todo <- grid[!file.exists(file)]
cat(sprintf("%d of %d cohorts to run\n", nrow(todo), nrow(grid)))
res <- parallel::mclapply(
  seq_len(nrow(todo)),
  function(i) {
    g <- todo[i]
    cfg <- cfg_of(g$frail, g$unrec, g$pf, g$dai)
    z <- sim(g$seed, cfg)
    saveRDS(
      list(
        code = BIAS$code,
        est = rbind(
          analyse(z, g$seed, cfg)[, cens := 1L],
          analyse_uncensored(z, g$seed, cfg)[, cens := 0L]
        ),
        sel = describe(z, g$seed, cfg)$sel
      ),
      g$file
    )
    return(TRUE)
  },
  mc.cores = CFG$n_core,
  mc.preschedule = FALSE
)
stopifnot(vapply(res, isTRUE, logical(1)))
x <- lapply(grid$file, readRDS)
stopifnot(vapply(x, function(r) identical(r$code, BIAS$code), logical(1)))
est <- rbindlist(Map(
  function(r, i) cbind(r$est, grid[i, .(frail, unrec, pf, dai)]),
  x,
  seq_len(nrow(grid))
))
# The frailty of the pre-Omicron and later infected, for the selection
# settings only.
i_sel <- which(grid$frail == 1L & grid$unrec == 1L & grid$dai == 1L)
sel <- rbindlist(lapply(i_sel, function(i) cbind(x[[i]]$sel, grid[i, .(pf)])))

# PART 2 -- ANALYSIS ====
true_w <- data.table(
  window = 1:5,
  true_lg = colMeans(log(CFG$tm[CFG$outcomes, ]))
)
true_lg <- true_w[window == 5L, true_lg]
full <- est[
  frail == 1L & unrec == 1L & pf == CFG$pre_frailty & dai == 1L & cens == 1L,
  .(window, outcome, hr, se, ew, seed)
]
stopifnot(isTRUE(all.equal(
  full[order(seed, outcome, window)],
  readRDS("results/run.rds")$est[order(seed, outcome, window)],
  check.attributes = FALSE
)))
# lg: the mean log estimate of the 12 outcomes, per window and seed.
g <- est[
  window %in% BIAS$windows,
  .(n = .N, lg = mean(log(hr))),
  keyby = .(window, frail, unrec, pf, dai, cens, seed)
]
stopifnot(all(g$n == length(CFG$outcomes)))
g5 <- g[window == 5L]

## Selection ----
fr <- sel[, .(u = sum(su) / sum(n)), keyby = .(pf, grp)]
fr <- dcast(fr, pf ~ grp, value.var = "u")
sel_tab <- g5[
  frail == 1L & unrec == 1L & dai == 1L & cens == 1L,
  .(gm = exp(mean(lg)), gm_se = sd(lg) / sqrt(.N)),
  keyby = pf
]
sel_tab[fr, on = "pf", lower := 1 - i.pre / i.later]

## Decomposition ----
v <- g[frail == as.integer(pf > 0) & (pf == 0 | pf == CFG$pre_frailty)]
v[true_w, on = "window", bias := i.true_lg - lg]
v[, risk := frail]
v[, k := paste0(risk, unrec, dai, cens)]
key <- function(on) paste(as.integer(BIAS$src %in% on), collapse = "")
n <- length(BIAS$src)
shp <- rbindlist(lapply(BIAS$windows, function(w) {
  rbindlist(lapply(CFG$seeds, function(s) {
    b <- setNames(
      v[window == w & seed == s, bias],
      v[window == w & seed == s, k]
    )
    stopifnot(length(b) == 2L^n)
    return(rbindlist(lapply(BIAS$src, function(i) {
      others <- setdiff(BIAS$src, i)
      phi <- 0
      for (m in 0:(n - 1L)) {
        subsets <- if (m == 0L) {
          list(character(0))
        } else {
          utils::combn(others, m, simplify = FALSE)
        }
        for (S in subsets) {
          wt <- factorial(m) * factorial(n - m - 1L) / factorial(n)
          phi <- phi + wt * (b[[key(c(S, i))]] - b[[key(S)]])
        }
      }
      return(data.table(window = w, seed = s, source = i, phi = phi))
    })))
  }))
}))
tot <- v[,
  .(on = bias[k == "1111"], off = bias[k == "0000"]),
  keyby = .(window, seed)
]
stopifnot(isTRUE(all.equal(
  shp[, sum(phi), keyby = .(window, seed)]$V1,
  tot[, on - off]
)))
shp_tab <- shp[,
  .(phi = mean(phi), phi_se = sd(phi) / sqrt(.N)),
  keyby = .(window, source)
]
shp_tab[, share := phi / sum(phi), by = window]
gm_tab <- tot[, .(on = exp(-mean(on)), off = exp(-mean(off))), keyby = window]
gm_tab[
  true_w,
  on = "window",
  `:=`(
    true = exp(i.true_lg),
    on = exp(i.true_lg) * on,
    off = exp(i.true_lg) * off
  )
]
# The gap log(true) - log(all on) is the residual log(true) - log(all off),
# from the analysis itself, plus the four Shapley contributions.
gm_tab[, `:=`(gap = log(true / on), residual = log(true / off))]

# PART 3 -- OUTPUT ====
cat(
  "\nGeometric mean of the 12 outcomes per window: true, all sources on, all off; the log gap true - on, and its residual true - off\n"
)
print(knitr::kable(gm_tab, format = "pipe", digits = 3))
cat("\nShapley contributions per window (phi, log scale); share is of the change from all sources off to all on, not of the gap from the true value:\n")
print(knitr::kable(shp_tab[order(window, -phi)], format = "pipe", digits = 3))
cat(
  "\nSelection: geometric mean at 12 months or more, and how much lower the mean frailty of the pre-Omicron infected is than that of the later infected:\n"
)
print(knitr::kable(sel_tab, format = "pipe", digits = 3))
## Steps ----
# Real simulated settings, adding one source at a time in a fixed order: the
# paper's analysis alone (all sources off, not censored), then censoring,
# unrecorded infections, death at infection, unmeasured risk without
# selection, and selection of increasing strength. The order is a choice: it
# changes the size of each step, not the first or last point.
STEPS <- list(
  list(lab = "Analysis alone", frail = 0L, unrec = 0L, pf = 0, dai = 0L, cens = 0L),
  list(lab = "+ stop at first CVD", frail = 0L, unrec = 0L, pf = 0, dai = 0L, cens = 1L),
  list(lab = "+ unrecorded infections", frail = 0L, unrec = 1L, pf = 0, dai = 0L, cens = 1L),
  list(lab = "+ deaths at infection", frail = 0L, unrec = 1L, pf = 0, dai = 1L, cens = 1L),
  list(lab = "+ unmeasured risk", frail = 1L, unrec = 1L, pf = 0, dai = 1L, cens = 1L)
)
for (p in BIAS$pf[BIAS$pf > 0]) {
  STEPS[[length(STEPS) + 1L]] <- list(
    lab = sprintf("+ selection, %.0f%%", 100 * sel_tab[pf == p, lower]),
    frail = 1L,
    unrec = 1L,
    pf = p,
    dai = 1L,
    cens = 1L
  )
}
steps <- rbindlist(lapply(seq_along(STEPS), function(k) {
  s <- STEPS[[k]]
  x <- g[
    frail == s$frail & unrec == s$unrec & abs(pf - s$pf) < 1e-9 & dai == s$dai & cens == s$cens,
    .(gm = exp(mean(lg)), gm_se = sd(lg) / sqrt(.N), n = .N),
    keyby = window
  ]
  stopifnot(all(x$n == length(CFG$seeds)))
  return(x[, `:=`(step = s$lab, k = k, n = NULL)])
}))
steps <- rbind(
  true_w[window %in% BIAS$windows, .(window, gm = exp(true_lg), gm_se = 0, step = "True", k = 0L)],
  steps
)
STEP_LEV <- c("True", vapply(STEPS, `[[`, character(1), "lab"))
cat("\nSteps: geometric mean of the 12 outcomes per window, adding one source at a time:\n")
print(knitr::kable(dcast(steps, k + step ~ window, value.var = "gm")[order(k)], format = "pipe", digits = 3))

saveRDS(
  list(
    sel_tab = sel_tab,
    shp_tab = shp_tab,
    gm_tab = gm_tab,
    shp = shp,
    tot = tot,
    true_gm = exp(true_lg),
    steps = steps,
    step_lev = STEP_LEV
  ),
  "results/bias.rds"
)

## Figure: from the true hazard ratios to the estimates, step by step ----
WLAB <- c(
  "2" = "Day 2 to <1 month",
  "3" = "1 to 5 months",
  "4" = "6 to 11 months",
  "5" = "12 months or more"
)
pd <- copy(steps)
pd[, `:=`(
  lo = gm * exp(-1.96 * gm_se),
  hi = gm * exp(1.96 * gm_se),
  step = factor(step, levels = STEP_LEV),
  wl = factor(WLAB[as.character(window)], levels = WLAB),
  what = fifelse(k == 0L, "True value", "Simulated estimate (95% interval)")
)]
pub <- readRDS("results/run.rds")$pooled[window %in% BIAS$windows, .(window, pub)]
pub[, wl := factor(WLAB[as.character(window)], levels = WLAB)]
q <- ggplot(pd, aes(x = step, y = gm))
q <- q + geom_hline(yintercept = 1, colour = "grey40")
q <- q + geom_hline(data = pub, aes(yintercept = pub, linetype = "Published estimate"))
q <- q + geom_line(data = pd[k > 0L], aes(group = 1), colour = "#1b9e77", linewidth = 0.7)
q <- q + geom_pointrange(aes(ymin = lo, ymax = hi, colour = what), size = 0.4)
q <- q + scale_colour_manual(NULL, values = c("True value" = "#d7191c", "Simulated estimate (95% interval)" = "#1b9e77"))
q <- q + scale_linetype_manual(NULL, values = c("Published estimate" = "dotted"))
q <- q + scale_x_discrete(limits = STEP_LEV)
q <- q + scale_y_continuous(trans = "log2", breaks = function(l) pretty(l, n = 5))
q <- q + facet_wrap(~wl, ncol = 1, scales = "free_y", axes = "all", axis.labels = "all_y")
q <- q + labs(x = NULL, y = "Hazard ratio (geometric mean of 12 outcomes)")
q <- q +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(3.5, "pt"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )
ggsave(
  "figures/bias_steps.png",
  q,
  width = 190,
  height = 270,
  units = "mm",
  dpi = 200,
  bg = "white"
)
