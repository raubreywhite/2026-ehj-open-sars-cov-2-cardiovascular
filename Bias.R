# Which sources of bias move the estimates below the true hazard ratios, and
# how much avoidance of infection by persons at higher underlying risk is
# needed? It uses the settings of Run.R, so run Run.R first.
#
# Avoidance: the model of Run.R with avoid 0, 0.02, 0.04, 0.06 and 0.08 (the
# value of Run.R). Per value, the geometric mean of the 12 estimates per
# window, analysed as the paper does, and the mean frailty of test-negative
# person-time in the second half of 2022 relative to the whole cohort.
#
# Decomposition: a Shapley decomposition of the bias, log(true geometric mean)
# - log(estimated geometric mean), per window, over four sources, each switched
# on or off. In day 0-1 some outcomes have fewer than 5 events in some cohorts
# and no estimate, so each window uses the outcomes with an estimate in every
# cohort and setting: some of the 12 in day 0-1, all 12 in the other windows.
#   risk:  unmeasured differences in risk, with avoidance of infection by
#          persons at higher risk (off: frailty SD 1e-6 and avoid 0; an SD of
#          0 divides by 0);
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
# It draws figures/bias_steps_w1.png to _w5.png: per window, the estimate as
# the sources are added one at a time.
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
BIAS$av <- c(0, 0.02, 0.04, 0.06, CFG$avoid)
BIAS$src <- c("risk", "unrec", "dai", "cens")
BIAS$dir <- "results/bias"
BIAS$windows <- 1:5
# The code that a saved cohort comes from: the md5 of the parsed code of Run.R
# and R/functions.R, without comments, so a change to a comment keeps the
# saved cohorts. The variant cohorts also depend on the variant analysis in
# this file.
code_md5 <- function(files) {
  tf <- tempfile()
  writeLines(unlist(lapply(files, function(f) {
    vapply(parse(f, keep.source = FALSE), function(e) paste(deparse(e, width.cutoff = 500L), collapse = "\n"), "")
  })), tf)
  return(unname(tools::md5sum(tf)))
}
BIAS$code <- code_md5(c("Run.R", "R/functions.R"))
BIAS$vcode <- code_md5(c("Run.R", "R/functions.R", "Bias.R"))
# Variant periods by the date of the positive test (paper, Methods), in days
# since 1 March 2020: original strain 1 February-31 December 2020, alpha
# 15 March-30 June 2021, delta 15 July-15 November 2021, omicron 28 December
# 2021-31 January 2022.
BIAS$variant <- data.table(
  variant = c("Original", "Alpha", "Delta", "Omicron"),
  lo = as.numeric(as.Date(c("2020-02-01", "2021-03-15", "2021-07-15", "2021-12-28")) - as.Date("2020-03-01")),
  hi = as.numeric(as.Date(c("2020-12-31", "2021-06-30", "2021-11-15", "2022-01-31")) - as.Date("2020-03-01"))
)
# Supplementary Table 5: hazard ratios for 1-12 months after the start of
# follow-up, by variant, for the 12 outcomes.
# fmt: skip
BIAS$pub_variant <- fread(header = TRUE, text = "
outcome Original Alpha Delta Omicron
ischemic_heart_disease 0.86 0.94 0.89 0.95
cerebral_infarction 0.99 0.91 0.99 0.88
cerebrovascular_hemorrhage 0.90 0.85 1.24 0.92
other_cerebrovascular 1.09 0.86 1.03 1.02
venous_embolism 1.06 1.22 0.99 1.08
pulmonary_embolism 1.51 1.69 1.43 0.89
aneurysm_dissection 1.08 1.04 0.79 1.04
heart_failure 0.61 0.46 1.30 0.81
inflammatory_heart_disease 1.52 1.46 1.22 0.92
arrhythmias 1.10 0.92 1.19 1.04
conduction_disorders 0.79 1.38 0.66 1.18
valve_disorders 0.89 1.12 1.10 1.01
")
cfg_of <- function(frail, unrec, av, dai) {
  return(modifyList(
    CFG,
    list(
      fsd = if (frail) CFG$fsd else 1e-6,
      contam = if (unrec) CFG$contam else 1e-6,
      avoid = av,
      p80 = if (dai) CFG$p80 else 0
    )
  ))
}

# PART 1 -- DATA CREATION ====
# The settings: the avoidance values with all other sources on, and the 8
# combinations of risk, unrec and dai. Risk off means no frailty and no
# avoidance. A setting in both lists appears once.
set_av <- data.table(frail = 1L, unrec = 1L, av = BIAS$av, dai = 1L)
set_shp <- CJ(risk = 0:1, unrec = 0:1, dai = 0:1)[, .(
  frail = risk,
  unrec,
  av = risk * CFG$avoid,
  dai
)]
grid <- unique(rbind(set_av, set_shp))
grid <- grid[, .(seed = CFG$seeds), by = names(grid)]
grid[,
  file := sprintf(
    "%s/f%d_u%d_av%.3f_d%d_seed%02d.rds",
    BIAS$dir,
    frail,
    unrec,
    av,
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
    cfg <- cfg_of(g$frail, g$unrec, g$av, g$dai)
    z <- sim(g$seed, cfg)
    saveRDS(
      list(
        code = BIAS$code,
        est = rbind(
          analyse(z, g$seed, cfg)[, cens := 1L],
          analyse_uncensored(z, g$seed, cfg)[, cens := 0L]
        ),
        negu = describe(z, g$seed, cfg)$negu
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
  function(r, i) cbind(r$est, grid[i, .(frail, unrec, av, dai)]),
  x,
  seq_len(nrow(grid))
))
# The frailty of test-negative person-time by half-year, for the avoidance
# settings only.
i_av <- which(grid$frail == 1L & grid$unrec == 1L & grid$dai == 1L)
negu <- rbindlist(lapply(i_av, function(i) cbind(x[[i]]$negu, grid[i, .(av)])))

# PART 2 -- ANALYSIS ====
# The outcomes of each window: those with an estimate in every cohort and
# setting, analysed both ways.
oset <- est[, .N, by = .(window, outcome)][N == 2L * nrow(grid), .(window, outcome)]
oset[, true := CFG$tm[cbind(outcome, c("day1", "m_lt1", "m_1_5", "m_6_11", "m_ge12")[window])]]
oset[CFG$pub, on = .(outcome, window), pub := i.pub]
true_w <- oset[, .(true_lg = mean(log(true)), pub = exp(mean(log(pub))), n_out = .N), keyby = window]
stopifnot(true_w[window > 1L, all(n_out == length(CFG$outcomes))])
true_lg <- true_w[window == 5L, true_lg]
full <- est[
  frail == 1L & unrec == 1L & av == CFG$avoid & dai == 1L & cens == 1L,
  .(window, outcome, hr, se, ew, seed)
]
stopifnot(isTRUE(all.equal(
  full[order(seed, outcome, window)],
  readRDS("results/run.rds")$est[order(seed, outcome, window)],
  check.attributes = FALSE
)))
# lg: the mean log estimate of the window's outcomes, per window and seed.
g <- est[oset, on = .(window, outcome)][,
  .(n = .N, lg = mean(log(hr))),
  keyby = .(window, frail, unrec, av, dai, cens, seed)
]
g[true_w, on = "window", n_out := i.n_out]
stopifnot(all(g$n == g$n_out))

## Avoidance ----
# In the last half-year, 28 August to 31 December 2022 (days 910-1036), and
# as in analyse() up to the first CVD diagnosis: higher, how much higher the
# mean frailty of test-negative person-time is than the cohort mean; vs_inf,
# its ratio to the mean frailty of person-time after a recorded infection.
fr <- negu[half == max(half), .(higher = sum(su) / sum(st) / mean(u_all) - 1, vs_inf = (sum(su) / sum(st)) / (sum(su_inf) / sum(st_inf))), keyby = av]
av_tab <- g[
  frail == 1L & unrec == 1L & dai == 1L & cens == 1L,
  .(gm = exp(mean(lg))),
  keyby = .(av, window)
]
av_tab <- dcast(av_tab, av ~ window, value.var = "gm")
setnames(av_tab, as.character(BIAS$windows), paste0("w", BIAS$windows))
av_tab[fr, on = "av", `:=`(higher = i.higher, vs_inf = i.vs_inf)]

## Decomposition ----
v <- g[frail == as.integer(av > 0) & (av == 0 | av == CFG$avoid)]
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
  "\nGeometric mean of the outcomes of each window: true, all sources on, all off; the log gap true - on, and its residual true - off\n"
)
print(knitr::kable(gm_tab, format = "pipe", digits = 3))
cat("\nShapley contributions per window (phi, log scale); share is of the change from all sources off to all on, not of the gap from the true value:\n")
print(knitr::kable(shp_tab[order(window, -phi)], format = "pipe", digits = 3))
cat(
  "\nAvoidance: geometric mean per window; in the last half-year, the frailty of the comparison group against the cohort mean (higher) and against person-time after a recorded infection (vs_inf):\n"
)
print(knitr::kable(av_tab, format = "pipe", digits = 3))
## Steps ----
# Real simulated settings, adding one source at a time in a fixed order: the
# four sources off (not censored), then censoring,
# unrecorded infections, death at infection, differences in underlying risk
# without avoidance, and avoidance of increasing strength. The order is a choice: it
# changes the size of each step, not the first or last point.
STEPS <- list(
  list(lab = "All biases off", frail = 0L, unrec = 0L, av = 0, dai = 0L, cens = 0L),
  list(lab = "+ stop at first diagnosis of any outcome", frail = 0L, unrec = 0L, av = 0, dai = 0L, cens = 1L),
  list(lab = "+ unrecorded infections", frail = 0L, unrec = 1L, av = 0, dai = 0L, cens = 1L),
  list(lab = "+ deaths caused by infection", frail = 0L, unrec = 1L, av = 0, dai = 1L, cens = 1L),
  list(lab = "+ differences in underlying risk", frail = 1L, unrec = 1L, av = 0, dai = 1L, cens = 1L)
)
for (p in BIAS$av[BIAS$av > 0]) {
  STEPS[[length(STEPS) + 1L]] <- list(
    lab = sprintf("with avoidance: underlying risk of comparison group %.2f x infected", av_tab[av == p, vs_inf]),
    frail = 1L,
    unrec = 1L,
    av = p,
    dai = 1L,
    cens = 1L
  )
}
steps <- rbindlist(lapply(seq_along(STEPS), function(k) {
  s <- STEPS[[k]]
  x <- g[
    frail == s$frail & unrec == s$unrec & abs(av - s$av) < 1e-9 & dai == s$dai & cens == s$cens,
    .(gm = exp(mean(lg)), h = stats::qt(0.975, .N - 1) * sd(lg) * sqrt(1 + 1 / .N), n = .N),
    keyby = window
  ]
  stopifnot(all(x$n == length(CFG$seeds)))
  return(x[, `:=`(step = s$lab, k = k, n = NULL)])
}))
steps <- rbind(
  true_w[window %in% BIAS$windows, .(window, gm = exp(true_lg), h = 0, step = "True", k = 0L)],
  steps
)
STEP_LEV <- c("True", vapply(STEPS, `[[`, character(1), "lab"))
cat(sprintf("\nSteps: geometric mean of the window's outcomes (%s per window), adding one source at a time:\n", paste(true_w$n_out, collapse = ", ")))
print(knitr::kable(dcast(steps, k + step ~ window, value.var = "gm")[order(k)], format = "pipe", digits = 3))

## Variant strata ----
# The main model of Run.R, analysed as the paper's Supplementary Table 5: per
# variant, the hazard ratio for 1-12 months after the start of follow-up
# (windows 1-5 months and 6-11 months together) of persons whose positive
# test falls in the variant's period, against test-negative time. Other
# exposed time is left out. The cohorts are simulated again, and saved to
# results/bias/ with the code md5.
by_variant <- function(z, seed, cfg) {
  # The test date: entry plus bio, which places a positive first test 30 days
  # before entry, as in the paper's design.
  test_day <- fifelse(z$d$seen, z$d$entry + z$d$bio, NA_real_)
  return(rbindlist(lapply(seq_len(nrow(BIAS$variant)), function(v) {
    inv <- !is.na(test_day) & test_day >= BIAS$variant$lo[v] & test_day < BIAS$variant$hi[v] + 1
    w <- z$win
    w[w > 0L] <- NA_integer_
    w[inv[row(z$win)] & z$win %in% 3:4] <- 1L
    zz <- z
    zz$win <- w
    return(analyse(zz, seed, cfg)[window == 1L][, variant := BIAS$variant$variant[v]])
  })))
}
vfile <- sprintf("%s/variant_seed%02d.rds", BIAS$dir, CFG$seeds)
vres <- parallel::mclapply(
  seq_along(CFG$seeds),
  function(i) {
    if (file.exists(vfile[i])) {
      return(readRDS(vfile[i]))
    }
    z <- sim(CFG$seeds[i])
    r <- list(code = BIAS$vcode, est = by_variant(z, CFG$seeds[i], CFG))
    saveRDS(r, vfile[i])
    return(r)
  },
  mc.cores = CFG$n_core,
  mc.preschedule = FALSE
)
stopifnot(vapply(vres, function(r) identical(r$code, BIAS$vcode), logical(1)))
vest <- rbindlist(Map(function(r, s) r$est[, seed := s], vres, CFG$seeds))
# Per variant, the outcomes with an estimate in every cohort (a small
# stratum can have fewer than 5 events); the published geometric mean is
# taken over the same outcomes.
vset <- vest[outcome %in% CFG$outcomes, .N, by = .(variant, outcome)][N == length(CFG$seeds), .(variant, outcome)]
var_tab <- vest[vset, on = .(variant, outcome)][, .(lg = mean(log(hr))), keyby = .(variant, seed)]
var_tab <- var_tab[, .(sim = exp(mean(lg))), keyby = variant]
pv <- melt(BIAS$pub_variant, id.vars = "outcome", variable.name = "variant", value.name = "pub", variable.factor = FALSE)
var_tab[pv[vset, on = .(variant, outcome)][, .(published = exp(mean(log(pub))), n_out = .N), keyby = variant], on = "variant", `:=`(published = i.published, n_out = i.n_out)]
var_tab[pv[, .(published_12 = exp(mean(log(pub)))), keyby = variant], on = "variant", published_12 := i.published_12]
var_tab <- var_tab[match(BIAS$variant$variant, variant)]
cat("\nVariant strata, 1-12 months, geometric mean of the 12 outcomes, simulated and published:\n")
print(knitr::kable(var_tab, format = "pipe", digits = 3))

saveRDS(
  list(
    var_tab = var_tab,
    av_tab = av_tab,
    shp_tab = shp_tab,
    gm_tab = gm_tab,
    shp = shp,
    tot = tot,
    true_gm = exp(true_lg),
    steps = steps,
    true_w = true_w,
    step_lev = STEP_LEV
  ),
  "results/bias.rds"
)

## Figures: from the true hazard ratios to the estimates, step by step ----
# One figure per window. Each estimate carries the 95% prediction interval of
# one study's geometric mean, from the spread over the seeds, as in
# compare_pub().
pd <- copy(steps)
pd[, `:=`(
  lo = gm * exp(-h),
  hi = gm * exp(h),
  step = factor(step, levels = STEP_LEV),
  what = fifelse(k == 0L, "True value", "Simulated estimate (95% prediction interval, one study)")
)]
for (w in BIAS$windows) {
  x <- pd[window == w]
  q <- ggplot(x, aes(x = step, y = gm))
  # The line at 1 only where the estimates come near it.
  if (min(x$lo) < 1.5) {
    q <- q + geom_hline(yintercept = 1, colour = "grey40")
  }
  q <- q + geom_hline(aes(yintercept = true_w[window == w, pub], linetype = "Published (geometric mean)"))
  q <- q + geom_line(data = x[k > 0L], aes(group = 1), colour = "#1b9e77", linewidth = 0.7)
  q <- q + geom_pointrange(aes(ymin = lo, ymax = hi, colour = what), size = 0.45)
  q <- q + scale_colour_manual(NULL, values = c("True value" = "#d7191c", "Simulated estimate (95% prediction interval, one study)" = "#1b9e77"))
  q <- q + scale_linetype_manual(NULL, values = c("Published (geometric mean)" = "dotted"))
  q <- q + scale_x_discrete(limits = STEP_LEV)
  q <- q + scale_y_continuous(trans = "log2", breaks = function(l) pretty(l, n = 6))
  q <- q + labs(x = NULL, y = "Hazard ratio, geometric mean")
  q <- q +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      plot.margin = margin(12, 12, 8, 12),
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
    sprintf("figures/bias_steps_w%d.png", w),
    q,
    width = 190,
    height = 130,
    units = "mm",
    dpi = 200,
    bg = "white"
  )
}
