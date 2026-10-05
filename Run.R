# Reproduces the letter on Boyd et al. 2026 (Eur Heart J Open,
# doi:10.1093/ehjopen/oeag121). Each of 100 seeds simulates a cohort of
# 4,508,489 synthetic persons under assumed true hazard ratios. Poisson models
# with 5-year age bands and calendar year analyse each cohort, as an
# approximation to the paper's Cox analysis.
#
# Per outcome and window, it compares the estimates with the published hazard
# ratios. It also analyses the cohorts by virus variant, as Supp Table 5 does.
# It simulates the same cohorts without infection, and gives the extra CVD
# diagnoses and the true risk and rate differences after a recorded infection.
# It adds the biases one at a time: unrecorded infections, depletion of
# susceptibles and avoidance of infection. It compares the seed-1 cohort with
# the paper, its supplement, Statens Serum Institut (SSI) data in data/ssi/
# and Erikstrup et al. 2022. It saves the results to results/run.rds. It draws
# figures/forest_12m.png, figures/forest_outcome_truths.png,
# figures/bias_steps.png and three figures/cohort_*.png.
#
# Needs R 4.6 with data.table, ggplot2 and patchwork.
# On Windows the run uses 1 worker.
# Run it from the repository root with: Rscript Run.R

library(data.table)
library(ggplot2)

# CFG ====
CFG <- list()

## Cohort (Boyd et al. 2026, Results and Table 1) ----
CFG$n <- 4508489L
CFG$seeds <- 1:100
CFG$n_core <- if (.Platform$OS.type == "windows") 1L else 2L
CFG$p_infected <- 2698261 / 4508489
CFG$age_median <- 34.2 # the reference age of the hazards
DAY_YR <- 365.25
STUDY <- 1036 # days from 1 March 2020 to 31 December 2022

# Age at the start of follow-up: a monotone spline through quantile knots. The
# knots are the quartiles 18.5, 34.2 and 52.3, and the Table 1 band edges at 40
# and 65 years. Follow-up starts 30 days after the first test, so a band edge
# at age a is a knot at a + 30 / 365.25. The knots at 0, 72, 97 and 105 years
# are assumed.
CFG$age_knots <- list(
  p = c(
    0,
    0.25,
    0.5,
    2603206 / 4508489,
    0.75,
    4053204 / 4508489,
    0.95,
    0.999,
    1
  ),
  age = c(0, 18.5, 34.2, 40 + 30 / 365.25, 52.3, 65 + 30 / 365.25, 72, 97, 105)
)
# Entry day: a monotone spline through the quartiles that the paper's
# follow-up quartiles give (27.5, 25.2 and 21.7 months), and a fitted tail.
CFG$entry_knots <- list(
  p = c(0, 0.25, 0.5, 0.75, 0.95, 1),
  day = c(
    0,
    1036 - 27.5 * 30.4375,
    1036 - 25.2 * 30.4375,
    1036 - 21.7 * 30.4375,
    759.73,
    1035
  )
)
# Recorded infection by age: SSI confirmed cases per resident, 2020-W10 to
# 2022-W52, per age group. edges are the lower age edges of the groups.
CFG$inf_age <- list(
  edges = c(0, 3, 6, 12, 16, 20, 40, 65, 80),
  rr = c(
    0.300438,
    0.585972,
    0.784431,
    0.787726,
    0.769325,
    0.679085,
    0.579198,
    0.343298,
    0.280519
  )
)

# Infection dates are piecewise uniform, in days since 1 March 2020. A share
# f_pre falls on (90, o0), before the Omicron wave, and f_late on (740, 1036)
# by calendar month. The rest is the wave, on (o0, 740). Day 740 ends 10 March
# 2022, when widespread testing ended. f_pre, o0 and f_late are fitted jointly
# with the entry tail to the person-years of Supp Tables 1 and 7.
CFG$f_pre <- 0.163861
CFG$o0 <- 660.5955
CFG$f_late <- 0.038642
# 11 March 2022, then the first day of each month, April 2022 to January 2023.
CFG$late_edges <- c(740, 761, 791, 822, 852, 883, 914, 944, 975, 1005, 1036)
# late_n: SSI first infections per month, from 11 March 2022. Source:
# data/ssi/24_reinfektioner_daglig_region.csv, 2022, all 5 regions.
# fmt: skip
CFG$late_n <- c(
  117294, 51046, 18319, 30252, 46303, 29218, 19935, 24852, 10590, 19262
)
# Hidden infections after day 740 are dated from the wastewater index instead.
# ww_edges: ISO weeks 10-52 of 2022, clipped to (740, 1036). Week 1 starts on
# day 673.
CFG$ww_edges <- pmin(pmax(673 + 7 * (9:52), 740), 1036)
# ww_w: SSI rna_mean_faeces, national, per week. Source:
# data/ssi/2026-09-23_dk_wastewater_data.csv, 2022-W10 to 2022-W52.
# fmt: skip
CFG$ww_w <- c(
  19380.62, 16403.12, 13474.51, 8960.46, 5383.24, 6189.51, 3864.95, 2177.80,
  886.66, 665.00, 789.17, 703.63, 844.51, 861.59, 2104.56, 3278.83, 1654.76,
  4295.19, 5494.91, 5100.89, 4191.00, 5012.95, 5279.03, 3557.26, 3003.98,
  2185.26, 2370.14, 1946.78, 2571.02, 2910.75, 5318.46, 6146.41, 5009.72,
  3684.16, 2137.53, 1881.64, 1808.10, 2734.82, 4134.47, 9195.75, 12283.69,
  11186.86, 13819.55
)
# Cumulative shares at the edges. Linear interpolation on them gives a date.
CFG$inf_x <- c(90, CFG$o0, CFG$late_edges)
CFG$inf_cm <- c(
  0,
  CFG$f_pre,
  1 - CFG$f_late,
  1 - CFG$f_late + CFG$f_late * cumsum(CFG$late_n) / sum(CFG$late_n)
)
CFG$ww_cm <- cumsum(c(0, CFG$ww_w * diff(CFG$ww_edges)))
CFG$ww_cm <- CFG$ww_cm / CFG$ww_cm[length(CFG$ww_cm)]

## Baseline rates per 1000 person-years, test-negative group (Supp Table 1) ----
# The first events per outcome in the row "No positive test", over 6,188,401
# person-years. The events sum to 54,247.
CFG$rate <- 1000 *
  c(
    ischemic_heart_disease = 8415,
    cerebral_infarction = 6064,
    cerebrovascular_hemorrhage = 1503,
    other_cerebrovascular = 4287,
    venous_embolism = 4926,
    arterial_embolism = 201,
    pulmonary_embolism = 2425,
    aneurysm_dissection = 2370,
    cardiac_arrest = 833,
    heart_failure = 1902,
    cardiomyopathy = 424,
    inflammatory_heart_disease = 1037,
    arrhythmias = 14039,
    conduction_disorders = 1769,
    valve_disorders = 4052
  ) /
  6188401
CFG$rate_true <- 1000 * 54247 / 6188401
# sc: the baseline scale of every outcome hazard. It is fitted jointly with the
# true hazard ratios, so that the full model gives the published test-negative
# first-event rate, rate_true (checks_v05/refit3.R in the internal
# repository).
CFG$sc <- 0.1579235837242475

## The simulated world ----
# fsd: SD of log frailty, a shared unrecorded cardiovascular risk.
# p80: the probability that infection kills an 80-year-old of average frailty
# (u = 1); it is multiplied by the person's frailty u.
# detect_early, detect_late: timing weights p of unrecorded infections, before
# and after widespread testing ended on day 740. Unrecorded infections follow
# the dates of recorded ones weighted by (1 - p) / p. contam sets their number.
# avoid: avoidance of infection by persons at higher underlying risk. A person
# 1 SD of log frailty more frail has exp(-avoid) times the weight in the draws
# of recorded and unrecorded infections.
# All five are assumed. No source measures them.
CFG$fsd <- 1.65
CFG$p80 <- 0.03
CFG$detect_early <- 0.40
CFG$detect_late <- 0.05
CFG$avoid <- 0.08

# contam: the share of never-recorded persons with an unrecorded infection. One
# third of infections were unrecorded (Erikstrup et al. 2022).
CFG$contam <- (CFG$p_infected / (1 - 1 / 3) - CFG$p_infected) /
  (1 - CFG$p_infected)

# The assumed true hazard ratios, one row per outcome, fitted so the simulated
# estimates reproduce the published ones. The windows start at 0, 1, 30,
# 182.62 and 365.25 days since infection. Arterial embolism, cardiac arrest
# and cardiomyopathy have the 3 lowest rates and are not among the 12 fitted
# outcomes. Their row is the geometric mean of the 12 fitted rows, per window.
# The fit is checks_v05/refit3.R in the internal repository.
EDGE <- c(0, 1, 30, 182.62, 365.25)
# The values are the fitted ones, to full precision.
# fmt: skip
CFG$tm <- matrix(c(
  8.0000382240960626, 1.2554860317102041, 1.2554860317102041, 1.2554860317102041, 1.2159019205112402,
  12.827617110464306, 1.3833851934048165, 1.146402217020321, 1.146402217020321, 1.01,
  7.8248079211829378, 1.4739622875993277, 1.1734746167983776, 1.1734746167983776, 1.0696503385599083,
  12.103577112908974, 1.5595265397061977, 1.4104638079726644, 1.4104638079726644, 1.1205924349453964,
  4.6360232770273893, 2.0019057864785603, 1.5974484219621812, 1.5974484219621812, 1.3362842926791776,
  8.6477865519139456, 1.6178171884012815, 1.3429671482158905, 1.3162919373516657, 1.1875353829080393,
  27.793038008044036, 7.4577626291016159, 2.0306040987397789, 1.6347830623133339, 1.4710666493484521,
  3.9688444176142879, 1.2750801485041992, 1.2750801485041992, 1.2750801485041992, 1.01,
  8.6477865519139456, 1.6178171884012815, 1.3429671482158905, 1.3162919373516657, 1.1875353829080393,
  15.141859000422764, 1.211065032683613, 1.0344613278947299, 1.01, 1.01,
  8.6477865519139456, 1.6178171884012815, 1.3429671482158905, 1.3162919373516657, 1.1875353829080393,
  5.5155226996984394, 1.2494151323788241, 1.2494151323788241, 1.2494151323788241, 1.1989293359197164,
  11.37036760861678, 1.5482729786936489, 1.4968619678690311, 1.4968619678690311, 1.3105672251880738,
  12.008318463790717, 1.3900738307862621, 1.3900738307862621, 1.3900738307862621, 1.3900738307862621,
  3.086496450834844, 1.2988825704687255, 1.2988825704687255, 1.2988825704687255, 1.2200694295047818
), 15, byrow = TRUE, dimnames = list(names(CFG$rate), c("day1", "m_lt1", "m_1_5", "m_6_11", "m_ge12")))

## Analysis ----
# 5-year age bands, and calendar year split at 1 January 2021, 2022 and 2023.
CFG$ageband <- 5
CAL_EDGES <- c(-1, 306, 671, 1036, 1e6)

## Published hazard ratios (Supp Table 1) ----
# The 12 outcomes with the highest published baseline rate, highest first. Per
# window from day 1 to 12 months or more, a row holds the HR and its 95% CI.
pub <- fread(
  header = FALSE,
  text = "
arrhythmias                7.67 6.17 9.53  1.19 1.07 1.32  1.07 1.02 1.12  1.05 1.00 1.11  0.90 0.82 0.99
ischemic_heart_disease     5.84 4.21 8.11  1.01 0.87 1.17  0.93 0.87 1.00  1.03 0.93 1.14  0.89 0.78 1.00
cerebral_infarction        8.29 5.97 11.5  1.16 0.99 1.37  0.88 0.81 0.96  0.91 0.83 0.99  0.74 0.62 0.88
venous_embolism            2.67 1.48 4.84  1.42 1.22 1.66  1.12 1.03 1.21  1.09 1.00 1.18  0.90 0.77 1.05
other_cerebrovascular      8.75 5.95 12.9  1.15 0.94 1.40  1.01 0.92 1.11  1.02 0.93 1.12  0.77 0.64 0.94
valve_disorders            2.51 1.20 5.28  1.03 0.83 1.28  0.98 0.89 1.08  1.03 0.93 1.14  0.86 0.71 1.04
pulmonary_embolism         16.1 11.0 23.5  4.14 3.57 4.79  1.21 1.08 1.36  0.89 0.78 1.02  0.77 0.61 0.98
aneurysm_dissection        3.85 1.83 8.09  0.90 0.68 1.21  1.02 0.90 1.15  0.98 0.86 1.11  0.73 0.56 0.95
heart_failure              10.8 6.51 18.0  1.00 0.73 1.37  0.87 0.73 1.01  0.73 0.62 0.87  0.57 0.41 0.80
conduction_disorders       8.85 5.01 15.6  0.76 0.53 1.09  1.08 0.94 1.24  1.13 0.98 1.30  1.13 0.88 1.44
cerebrovascular_hemorrhage 6.93 3.45 13.9  1.19 0.87 1.63  0.86 0.73 1.01  0.96 0.82 1.13  0.85 0.63 1.14
inflammatory_heart_disease 4.90 2.04 11.8  1.08 0.76 1.54  0.99 0.84 1.18  1.05 0.89 1.23  0.91 0.67 1.23
"
)
CFG$pub <- melt(
  pub,
  id.vars = 1L,
  measure.vars = list(seq(2, 14, 3), seq(3, 15, 3), seq(4, 16, 3)),
  variable.name = "window",
  value.name = c("pub", "lo", "hi")
)
setnames(CFG$pub, 1L, "outcome")
CFG$pub[, window := as.integer(window)]
CFG$outcomes <- pub$V1
rm(pub)

## Targets of the cohort comparison ----
# Paper, Results: age median 34.2 (IQR 18.5-52.3) years, 8,909,627
# person-years, follow-up median 25.2 (IQR 21.7-27.5) months, 70,885 persons
# with a first CVD diagnosis, 16,475 of them after a positive test. Paper,
# Table 1: persons by age at the first test (<40, 40-64, 65 or more), first
# test positive plus first test negative.
PAPER <- list(
  age_q = c(18.5, 34.2, 52.3),
  py = 8909627,
  fu_q = c(21.7, 25.2, 27.5),
  cvd_n = 70885,
  cvd_pos_n = 16475,
  age_n = c(131489 + 2471717, 61057 + 1388941, 24352 + 430933)
)
# Supp Tables 1 and 7, person-years per window: test-negative, day 1, <1 month,
# 1-5, 6-11 and 12 months or more. Table 1 runs to 31 December 2022, and Table
# 7 to 10 March 2022. Table 7 has 38,419 first events in test-negative time.
PUB_PY <- cbind(
  t1 = c(6188401, 6814, 201905, 1101947, 1115292, 295266),
  t7 = c(4732238, 6161, 163013, 242882, 118442, 50715)
)
PUB_EV7 <- 38419
WIN_LAB <- c(
  "Test-negative",
  "Day 1",
  "<1 month",
  "1-5 months",
  "6-11 months",
  "12+ months"
)
# Erikstrup et al. 2022: 66% (95% CI 63-70%) of all healthy blood donors aged
# 17-72 were infected from 1 November 2021 to 15 March 2022, days 610 to 744.
ERIK <- c(lo = 610, hi = 745)
# The SSI source files. data/ssi/README.md gives their origin.
SSI <- c(
  age = "data/ssi/18_fnkt_alder_uge_testede_positive_nyindlagte.csv",
  month = "data/ssi/24_reinfektioner_daglig_region.csv",
  ww = "data/ssi/2026-09-23_dk_wastewater_data.csv"
)
SSI_MD5 <- c(
  age = "58245bd78ff14fda0831f04009329575",
  month = "c3b069badfdaaf40b0d4f903d9abed26",
  ww = "f1ef3aff3ec7a4987c1f1a31c78fcc50"
)
# The first day of each month, March 2020 to December 2022.
MON <- seq(as.Date("2020-03-01"), as.Date("2022-12-01"), by = "month")
MON_DAY <- as.numeric(MON - MON[1])

## Figure legend ----
LAB_TRUE <- "True hazard ratio (assumed)"
LAB_PUB <- "Published estimate, Boyd et al. (95% CI)"
LAB_SIM <- "Biased estimate expected under the true hazard ratio (95% prediction interval)"
WL <- c(
  "Day 0-1",
  "Day 2 to <1 month",
  "1 to 5 months",
  "6 to 11 months",
  "12 months or more"
)

## Variant strata (paper, Methods, and Supp Table 5) ----
# Variant periods by the date of the positive test, in days since 1 March 2020:
#   original strain: 1 February-31 December 2020
#   alpha:           15 March-30 June 2021
#   delta:           15 July-15 November 2021
#   omicron:         28 December 2021-31 January 2022
CFG$variant <- data.table(
  variant = c("Original", "Alpha", "Delta", "Omicron"),
  lo = as.numeric(
    as.Date(c("2020-02-01", "2021-03-15", "2021-07-15", "2021-12-28")) -
      as.Date("2020-03-01")
  ),
  hi = as.numeric(
    as.Date(c("2020-12-31", "2021-06-30", "2021-11-15", "2022-01-31")) -
      as.Date("2020-03-01")
  )
)
# Supp Table 5: hazard ratios for 1-12 months after the start of follow-up,
# by variant, for the 12 outcomes.
# fmt: skip
CFG$pub_variant <- fread(header = TRUE, text = "
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

## Bias steps ----
# The steps add one bias at a time to a model with all three off. The order is
# unrecorded infections, depletion of susceptibles, then avoidance of infection
# at the strength of the main model. The last step is the main model. The order
# is a choice: it changes the size of each step, not the first or last point.
# Death at infection (p80) is on in every step. A bias is switched off as
# follows:
#   unrec:     contam 1e-6, about 2 persons. A share of 0 draws no dates and
#              stops sample().
#   depletion: frailty SD 1e-6, because an SD of 0 divides by 0. Each outcome
#              is followed to its own first diagnosis, death or fu (cens 0).
#   avoidance: avoid 0.
# CFG$bias holds the settings that need cohorts of their own. nodai is the
# main model without death at infection.
CFG$steps <- c(
  off = "All biases off",
  unrec = "+ unrecorded infections",
  depl = "+ depletion of susceptibles",
  full = "+ avoidance of infection"
)
CFG$bias <- data.table(
  setting = c("off", "unrec", "depl", "nodai"),
  frail = c(0L, 0L, 1L, 1L),
  unrec = c(0L, 1L, 1L, 1L),
  av = c(0, 0, 0, CFG$avoid),
  dai = c(1L, 1L, 1L, 0L),
  cens = c(0L, 0L, 1L, 1L)
)

# PART 1 -- DATA CREATION ====
## Functions ----
# The date at cumulative share p, on a curve with cumulative shares cm at x.
qcurve <- function(p, cm, x) stats::approx(cm, x, xout = p, rule = 2)$y

# Exactly k of length(w) units, drawn by systematic sampling over a random
# order. Unit i is drawn with probability min(c * w[i], 1), where c makes the
# probabilities sum to k. It returns the drawn indices.
pps_draw <- function(w, k) {
  pr <- numeric(length(w))
  free <- rep(TRUE, length(w))
  repeat {
    pr[free] <- w[free] * (k - sum(!free)) / sum(w[free])
    over <- free & pr > 1
    if (!any(over)) {
      break
    }
    pr[over] <- 1
    free[over] <- FALSE
  }
  o <- sample.int(length(w))
  cs <- cumsum(pr[o])
  cs <- cs * (k / cs[length(cs)])
  u <- runif(1)
  return(o[floor(cs - u) > floor(c(0, cs[-length(cs)]) - u)])
}

# The new persons for infections that must move, one per date t. For each
# date in turn, it draws a person at random from cand[[key]] who is free and
# alive on that date (last >= t). Where there is none, it draws from all free
# persons alive on that date, with probability proportional to w. A drawn
# person is no longer free. Rejection over 64 draws, then the full set, keeps
# the draw uniform.
move_to_living <- function(t, key, cand, free, last, w) {
  to <- integer(length(t))
  for (i in seq_along(t)) {
    cc <- cand[[key[i]]]
    j <- cc[sample.int(length(cc), min(length(cc), 64L), replace = TRUE)]
    j <- j[free[j] & last[j] >= t[i]]
    if (length(j) == 0L) {
      j <- cc[free[cc] & last[cc] >= t[i]]
      pw <- NULL
      if (length(j) == 0L) {
        j <- which(free & last >= t[i])
        pw <- w[j]
      }
      j <- j[sample.int(length(j), 1L, prob = pw)]
    }
    to[i] <- j[1L]
    free[to[i]] <- FALSE
  }
  return(to)
}

# One simulated cohort. Follow-up is cut into 10 segments per person: before
# infection, and one per exposure window on the analysis clock and on the
# biological clock. The hazard of each outcome is Gompertz in age, times
# frailty, times the baseline scale sc, times the outcome's hazard ratio in the
# segment's window. One Exp(1) threshold per person and outcome inverts in
# closed form to the event time.
sim <- function(seed, cfg = CFG, infect = TRUE, sc = cfg$sc) {
  n <- cfg$n
  set.seed(seed)
  lam <- cfg$rate / 1000 / DAY_YR # per day at the median age
  b <- log(2) / 7.5 # age slope of the outcome hazard, per year
  bd <- log(2) / 8 # age slope of the death hazard, per year

  ak <- cfg$age_knots
  ek <- cfg$entry_knots
  d <- data.table(
    age0 = stats::splinefun(ak$p, ak$age, method = "hyman")(runif(n))
  )
  d[, entry := stats::splinefun(ek$p, ek$day, method = "hyman")(runif(.N))]
  d[, fu := STUDY - entry]
  d[, u := exp(rnorm(.N, -cfg$fsd^2 / 2, cfg$fsd))]

  # Natural death: Gompertz in age, proportional to frailty, and independent
  # of infection. It comes before infection, so that every infection can fall
  # at or before it. last is its calendar day, Inf without a death in
  # follow-up. ki is the draw for death at infection below. The draws use
  # their own stream, and the main stream is restored after them.
  rs <- .GlobalEnv$.Random.seed
  Ad <- (DAY_YR / bd) * exp(bd * (d$age0 - cfg$age_median))
  kd <- (0.0008 / DAY_YR) * d$u
  Hd <- kd * (Ad * exp(bd * d$fu / DAY_YR) - Ad)
  set.seed(seed + 2000000L)
  ed <- rexp(n)
  ki <- runif(n)
  .GlobalEnv$.Random.seed <- rs
  tdn <- rep(NA_real_, n)
  die <- ed < Hd
  tdn[die] <- (DAY_YR / bd) * log((Ad[die] + ed[die] / kd[die]) / Ad[die])
  last <- d$entry + fcoalesce(tdn, Inf)

  # Recorded infections. Exactly round(n * p_infected) persons are infected.
  # pps_draw() draws each with a probability proportional to the rr of their
  # age group, times exp(-avoid * zf): persons at higher underlying risk are
  # less likely to be infected. zf is the standardised log frailty.
  ag <- findInterval(d$age0, cfg$inf_age$edges)
  rr <- cfg$inf_age$rr[ag]
  zf <- (log(d$u) + cfg$fsd^2 / 2) / cfg$fsd
  tilt <- exp(-cfg$avoid * zf)
  hit <- rep(FALSE, n)
  hit[pps_draw(rr * tilt, round(n * cfg$p_infected))] <- TRUE
  # A share f_pre of persons has a date before o0, at random. The other dates
  # come from the rest of the curve.
  pre <- runif(n) < cfg$f_pre
  rec <- rep(NA_real_, n)
  rec[hit & pre] <- runif(sum(hit & pre), 90, cfg$o0)
  rec[hit & !pre] <- qcurve(
    runif(sum(hit & !pre), cfg$f_pre, 1),
    cfg$inf_cm,
    cfg$inf_x
  )

  # Unrecorded infections. Their number is set by contam. Their dates follow
  # the fitted curve of recorded infections, weighted by (1 - p) / p, with p
  # the timing weight detect_early before day 740 and detect_late after. They
  # go to never-recorded persons by pps_draw(), with the same rr by age group
  # and the same tilt by underlying risk.
  cand <- which(is.na(rec))
  k <- round(length(cand) * cfg$contam)
  m <- 4L * k
  dt0 <- qcurve(runif(m), cfg$inf_cm, cfg$inf_x)
  pdet <- ifelse(dt0 < 740, cfg$detect_early, cfg$detect_late)
  dt0 <- dt0[sample(seq_len(m), k, replace = TRUE, prob = (1 - pdet) / pdet)]
  pick <- cand[pps_draw(rr[cand] * tilt[cand], k)]
  # A hidden date at or after day 740 is redrawn from the wastewater shape, on
  # its own stream. The main stream is restored, so later draws do not change.
  rs <- .GlobalEnv$.Random.seed
  set.seed(seed + 6000000L)
  late <- which(dt0 >= 740)
  dt0[late] <- qcurve(runif(length(late)), cfg$ww_cm, cfg$ww_edges)
  unr <- rep(NA_real_, n)
  unr[pick] <- dt0

  # Only a living person is infected. A recorded infection dated after natural
  # death moves, with its date, to a person of the same age group and the same
  # pre draw who is alive on that date and has no recorded infection. An
  # unrecorded infection moves in the same way, to a never-recorded person of
  # the same age group without one, when it falls after natural death or its
  # person took over a recorded infection. move_to_living() draws the new
  # persons, on their own stream.
  set.seed(seed + 3000000L)
  nag <- length(cfg$inf_age$rr)
  bad <- which(rec > last)
  moved <- c(rec = length(bad))
  tb <- rec[bad]
  rec[bad] <- NA_real_
  key <- ag + nag * pre
  free <- is.na(rec)
  cand <- split(which(free), factor(key[free], levels = seq_len(2L * nag)))
  rec[move_to_living(tb, key[bad], cand, free, last, rr)] <- tb
  bad <- which(!is.na(unr) & (!is.na(rec) | unr > last))
  moved <- c(moved, unr = length(bad))
  tb <- unr[bad]
  unr[bad] <- NA_real_
  free <- is.na(rec) & is.na(unr)
  cand <- split(which(free), factor(ag[free], levels = seq_len(nag)))
  unr[move_to_living(tb, ag[bad], cand, free, last, rr)] <- tb
  .GlobalEnv$.Random.seed <- rs

  # bio: infection on the biological clock, in days since entry. ana: the
  # analysis clock, which starts a recorded infection before entry at entry.
  d[, inf_cal := fcoalesce(rec, unr)]
  d[, bio := inf_cal - entry]
  d[!is.na(rec), bio := pmax(bio, -30)]
  d[!is.na(bio) & bio >= fu, bio := NA_real_]
  d[, ana := fifelse(is.na(rec), NA_real_, pmax(bio, 0))]
  d[, seen := !is.na(ana)]
  # infect = FALSE keeps every draw but removes every infection. Death and the
  # event thresholds come from their own streams, so the persons, their natural
  # death times and their thresholds stay the same as with infection.
  if (!infect) {
    d[, c("inf_cal", "bio", "ana") := NA_real_]
    d[, seen := FALSE]
  }

  # Segment edges in days since entry. A recorded infection splits follow-up
  # on the analysis clock, and an unrecorded one on the biological clock. A
  # recorded infection before entry also splits it on the biological clock,
  # so each segment lies in one window of each clock. A bubble sort puts the
  # edges in order.
  brk <- outer(d$ana, EDGE, "+")
  hid <- which(is.na(d$ana) & !is.na(d$bio))
  brk[hid, ] <- pmax(outer(d$bio[hid], EDGE, "+"), 0)
  s0 <- cbind(0, pmin(brk, d$fu))
  rm(brk)
  s0[is.na(s0)] <- 0
  s0 <- cbind(s0, d$fu)
  for (j in 2:ncol(s0)) {
    s0[, j] <- pmax(s0[, j], s0[, j - 1L])
  }
  # The edges bio + EDGE[-1] of a recorded infection before entry, clamped to
  # (0, fu). The edge bio + EDGE[1] clamps to 0, which is already column 1.
  # Other rows get 0, which only adds segments of length 0.
  neg <- which(!is.na(d$ana) & d$bio < 0)
  bb <- matrix(0, n, length(EDGE) - 1L)
  bb[neg, ] <- pmin(pmax(outer(d$bio[neg], EDGE[-1L], "+"), 0), d$fu[neg])
  s0 <- cbind(s0, bb)
  rm(bb)
  for (p in seq_len(ncol(s0) - 1L)) {
    for (j in seq_len(ncol(s0) - 1L)) {
      sw <- s0[, j] > s0[, j + 1L]
      tmp <- s0[sw, j]
      s0[sw, j] <- s0[sw, j + 1L]
      s0[sw, j + 1L] <- tmp
    }
  }
  mid <- (s0[, -1L] + s0[, -ncol(s0)]) / 2

  # Window 0 to 5 of each segment: the true hazard ratio reads the biological
  # clock, and the analysis reads the analysis clock. 0 is unexposed.
  wb <- findInterval(mid - d$bio, EDGE)
  wb[is.na(wb)] <- 0L
  win <- matrix(findInterval(mid - d$ana, EDGE), nrow(mid))
  win[is.na(win)] <- 0L

  # The cumulative Gompertz hazard over a segment [t0, t1] is
  # u * lam * sc * hr * (G(t1) - G(t0)), with G(t) = Acon * exp(b * t / DAY_YR).
  Acon <- (DAY_YR / b) * exp(b * (d$age0 - cfg$age_median))
  Glo <- Acon * exp(b * s0[, -ncol(s0)] / DAY_YR)
  Ghi <- Acon * exp(b * s0[, -1L] / DAY_YR)

  # Natural death, then death at infection with a probability that rises with
  # age and frailty. Every person was alive at entry, so death at infection
  # applies only to an infection during follow-up (bio >= 0).
  d[, tdeath := tdn]
  pcd <- cfg$p80 *
    exp(bd * (d$age0 - cfg$age_median)) /
    exp(bd * (80 - cfg$age_median)) *
    d$u
  kill <- !is.na(d$bio) & d$bio >= 0 & ki < pmin(pcd, 0.95)
  d[kill, tdeath := pmin(tdeath, bio, na.rm = TRUE)]

  # Event time of each outcome: where the cumulative hazard reaches E. Row k
  # of cfg$tm gives outcome k its hazard ratio per window.
  set.seed(seed + 1000000L)
  E <- matrix(rexp(n * length(lam)), nrow = n)
  for (k in seq_along(lam)) {
    kk <- d$u *
      lam[k] *
      sc *
      matrix(c(1, cfg$tm[k, ])[wb + 1L], nrow(mid))
    Hc <- kk * (Ghi - Glo)
    for (j in 2:ncol(Hc)) {
      Hc[, j] <- Hc[, j] + Hc[, j - 1L]
    }
    jj <- rowSums(Hc < E[, k]) + 1L # first segment that reaches E
    ii <- which(jj <= ncol(Hc))
    jj <- jj[ii]
    Hlo <- Hc[cbind(ii, pmax(jj - 1L, 1L))]
    Hlo[jj == 1L] <- 0
    tev <- rep(NA_real_, n)
    tev[ii] <- (DAY_YR / b) *
      log(
        (Glo[cbind(ii, jj)] + (E[ii, k] - Hlo) / kk[cbind(ii, jj)]) / Acon[ii]
      )
    d[, (names(lam)[k]) := tev]
  }
  return(list(d = d, s0 = s0, win = win, sc = sc, moved = moved))
}

# The first CVD diagnosis of each person, and whether it comes before death
# and the end of follow-up.
first_cvd <- function(d, cfg = CFG) {
  f <- do.call(pmin, c(d[, names(cfg$rate), with = FALSE], na.rm = TRUE))
  return(data.table(
    seen = d$seen,
    bio = d$bio,
    ana = d$ana,
    fu = d$fu,
    f = f,
    ev = !is.na(f) & f <= pmin(d$fu, d$tdeath, na.rm = TRUE)
  ))
}

# The excess of one seed: x1 is first_cvd() with infection, and x0 the same
# persons without it. It returns the excess persons with a CVD diagnosis by
# infection status with infection, and the risk of a first diagnosis in the
# 12 months after a recorded infection, among recorded persons followed that
# long, with and without infection.
excess <- function(x1, x0) {
  nv <- which(is.na(x1$bio))
  stopifnot(identical(x1[nv, .(f, ev)], x0[nv, .(f, ev)]))
  grp <- fifelse(
    x1$seen,
    "recorded",
    fifelse(is.na(x1$bio), "never", "unrecorded")
  )
  tot <- data.table(grp = grp, e = x1$ev - x0$ev)[,
    .(excess = sum(e)),
    keyby = grp
  ]
  tt <- DAY_YR
  j <- which(x1$seen & x1$ana + tt <= x1$fu)
  a <- x1$ana[j]
  risk <- function(x) mean(x$ev[j] & x$f[j] > a & x$f[j] <= a + tt)
  return(list(
    tot = tot,
    risk = c(n = length(j), with = risk(x1), without = risk(x0))
  ))
}

# The risk of a first diagnosis of each fitted outcome in one cohort d. j and a
# are the recorded persons and the start of their analysis clock, both from the
# cohort with infection. Window w runs from a + lo[w] to a + hi[w], with hi 730.5
# days (24 months) for 12 months or more. Its base is the persons in j with
# a + hi[w] <= fu. A diagnosis counts when it comes after a and at or before
# death and fu. Other outcomes do not censor it. Per outcome and window, it
# returns the base n and, on that base, the risk in (a, a + hi] (cum).
risk_after <- function(d, j, a, cfg = CFG) {
  hi <- c(EDGE[-1L], 2 * DAY_YR)
  fu <- d$fu[j]
  end <- pmin(fu, d$tdeath[j], na.rm = TRUE)
  retval <- list()
  for (o in cfg$outcomes) {
    tk <- d[[o]][j]
    ok <- !is.na(tk) & tk > a & tk <= end
    for (w in seq_along(hi)) {
      b <- a + hi[w] <= fu
      retval[[length(retval) + 1L]] <- data.table(
        outcome = o,
        window = w,
        n = sum(b),
        cum = mean(ok[b] & tk[b] <= a[b] + hi[w])
      )
    }
  }
  return(rbindlist(retval))
}

# The rate of a first diagnosis of each fitted outcome in one cohort d, per
# window since the analysis clock. j and a are as in risk_after(). Window w
# runs from a + lo[w] to a + hi[w], and 12 months or more runs to the end of
# follow-up. A person is at risk in a window from its start to the first
# diagnosis of the outcome, death, fu or the end of the window. Other outcomes
# do not censor it. Per outcome and window, it returns the first diagnoses ev
# and the person-years py.
rate_after <- function(d, j, a, cfg = CFG) {
  lo <- EDGE
  hi <- c(EDGE[-1L], Inf)
  end <- pmin(d$fu[j], d$tdeath[j], na.rm = TRUE)
  retval <- list()
  for (o in cfg$outcomes) {
    tk <- d[[o]][j]
    ok <- !is.na(tk) & tk <= end
    stop_t <- pmin(end, tk, na.rm = TRUE)
    for (w in seq_along(lo)) {
      s <- a + lo[w]
      retval[[length(retval) + 1L]] <- data.table(
        outcome = o,
        window = w,
        ev = sum(ok & tk > s & tk <= a + hi[w]),
        py = sum(pmax(pmin(stop_t, a + hi[w]) - s, 0)) / DAY_YR
      )
    }
  }
  return(rbindlist(retval))
}

# The rate differences per 100,000 person-years: k1 is rate_after() of the
# cohort with infection, and k0 of the same persons without it.
rate_diff <- function(k1, k0) {
  stopifnot(identical(k1[, .(outcome, window)], k0[, .(outcome, window)]))
  return(k1[, .(
    outcome,
    window,
    rd = 1e5 * (ev / py - k0$ev / k0$py)
  )])
}

# The risk differences in percentage points: k1 is risk_after() of the cohort
# with infection, and k0 of the same persons without it.
rd_pp <- function(k1, k0) {
  stopifnot(identical(k1[, .(outcome, window, n)], k0[, .(outcome, window, n)]))
  return(k1[, .(
    outcome,
    window,
    n,
    cum = 100 * (cum - k0$cum)
  )])
}

# The paper's analysis. Follow-up stops at the first cardiovascular event of
# any type. A Poisson model per outcome adjusts the exposure window for 5-year
# age band and calendar year. With cens = FALSE, the person-time of each
# outcome stops at that outcome's own first diagnosis, death or fu instead.
analyse <- function(z, seed, cfg = CFG, cens = TRUE) {
  d <- z$d
  tfirst <- do.call(pmin, c(d[, names(cfg$rate), with = FALSE], na.rm = TRUE))
  sadm <- pmin(d$fu, d$tdeath, na.rm = TRUE)
  stop_t <- if (cens) pmin(tfirst, sadm, na.rm = TRUE) else sadm

  # Person-time split at the window, age-band and calendar-year edges. A
  # stretch assigned whole to the band that holds its midpoint biases the
  # estimates, because stretches that end in an event are short.
  lo <- pmin(z$s0[, -ncol(z$s0)], stop_t)
  hi <- pmin(z$s0[, -1L], stop_t)
  keep <- as.vector(hi - lo > 0)
  p <- data.table(
    lo = as.vector(lo)[keep],
    hi = as.vector(hi)[keep],
    win = as.vector(z$win)[keep],
    row = rep(seq_len(nrow(d)), times = ncol(lo))[keep]
  )
  a0 <- (d$age0[p$row] + p$lo / DAY_YR) / cfg$ageband
  a1 <- (d$age0[p$row] + p$hi / DAY_YR) / cfg$ageband
  nb <- floor(a1 - 1e-9) - floor(a0) + 1L
  p <- p[rep(seq_len(.N), times = nb)]
  p[, ab := as.integer(floor(rep(a0, times = nb)) + sequence(nb) - 1L)]
  p[, s0 := pmax(lo, (ab * cfg$ageband - d$age0[row]) * DAY_YR)]
  p[, s1 := pmin(hi, ((ab + 1L) * cfg$ageband - d$age0[row]) * DAY_YR)]
  p <- p[s1 > s0]
  b0 <- findInterval(d$entry[p$row] + p$s0, CAL_EDGES)
  b1 <- findInterval(d$entry[p$row] + p$s1 - 1e-9, CAL_EDGES)
  nc <- as.integer(b1 - b0 + 1L)
  p <- p[rep(seq_len(.N), times = nc)]
  p[, cp := as.integer(rep(b0, times = nc) + sequence(nc) - 1L)]
  p[, q0 := pmax(s0, CAL_EDGES[cp] - d$entry[row])]
  p[, q1 := pmin(s1, CAL_EDGES[cp + 1L] - d$entry[row])]
  p <- p[q1 > q0]
  if (cens) {
    p[, py := (q1 - q0) / DAY_YR]
    p[, last := !duplicated(row, fromLast = TRUE)]
  }

  retval <- list()
  for (o in cfg$outcomes) {
    tk <- d[[o]]
    if (cens) {
      isev <- !is.na(tk) & !is.na(tfirst) & tk == tfirst & tk <= sadm
      p[, ev := as.integer(isev[row] & abs(q1 - stop_t[row]) < 1e-8 & last)]
      a <- p[, .(py = sum(py), ev = sum(ev)), keyby = .(win, ab, cp)]
    } else {
      st <- pmin(tk, sadm, na.rm = TRUE)
      isev <- !is.na(tk) & tk <= sadm
      po <- p[q0 < st[row]]
      po[, q1o := pmin(q1, st[row])]
      po[, py := (q1o - q0) / DAY_YR]
      po[, ev := as.integer(isev[row] & abs(q1o - st[row]) < 1e-8)]
      a <- po[, .(py = sum(py), ev = sum(ev)), keyby = .(win, ab, cp)]
    }
    ew <- a[, .(ew = sum(ev)), keyby = .(window = win)]
    a[, `:=`(win = factor(win), ab = factor(ab), cp = factor(cp))]
    fit <- stats::glm(
      ev ~ win + ab + cp + offset(log(py)),
      family = stats::poisson(),
      data = a
    )
    cf <- stats::coef(summary(fit))
    rn <- grep("^win", rownames(cf), value = TRUE)
    r <- data.table(
      outcome = o,
      window = as.integer(sub("^win", "", rn)),
      hr = exp(cf[rn, "Estimate"]),
      se = cf[rn, "Std. Error"]
    )
    # A window with fewer than 5 events gives no estimate.
    retval[[o]] <- merge(r, ew, by = "window")[ew >= 5L]
  }
  return(rbindlist(retval)[, seed := seed][])
}

# The analysis of Supp Table 5. Per variant, it estimates one hazard ratio for
# 1-12 months after the start of follow-up, the windows 1-5 months and 6-11
# months together. The exposed persons are those whose positive test falls in
# the variant's period, and the reference is test-negative time. Other exposed
# time is left out. The test date is entry plus bio, which places a positive
# first test 30 days before entry, as in the paper's design.
by_variant <- function(z, seed, cfg = CFG) {
  test_day <- fifelse(z$d$seen, z$d$entry + z$d$bio, NA_real_)
  return(rbindlist(lapply(seq_len(nrow(cfg$variant)), function(v) {
    inv <- !is.na(test_day) &
      test_day >= cfg$variant$lo[v] &
      test_day < cfg$variant$hi[v] + 1
    w <- z$win
    w[w > 0L] <- NA_integer_
    w[inv[row(z$win)] & z$win %in% 3:4] <- 1L
    zz <- z
    zz$win <- w
    return(analyse(zz, seed, cfg)[window == 1L][,
      variant := cfg$variant$variant[v]
    ])
  })))
}

# The estimates of all seeds against the published ones. Per outcome and
# window: mean m and SD s of the log HR over the k seeds. The 95% prediction
# interval of one study's estimate is exp(m +- qt(0.975, k - 1) * s * sqrt(1 +
# 1/k)). An estimate in fewer than 60% of the seeds gives no interval. Per seed
# and window, lg is the mean log HR over the outcomes of oset. gm is its mean
# over the seeds, exponentiated. sig counts, per seed, the outcomes whose 95% CI at 12 months or more
# lies below 1.
compare_pub <- function(est, oset, cfg = CFG) {
  ow <- est[,
    .(m = mean(log(hr)), s = sd(log(hr)), k = .N),
    keyby = .(outcome, window)
  ]
  ow <- ow[cfg$pub, on = .(outcome, window)]
  ow[, true := cfg$tm[cbind(outcome, colnames(cfg$tm)[window])]]
  ow[, h := qt(0.975, k - 1) * s * sqrt(1 + 1 / k)]
  ow[, `:=`(l95 = exp(m - h), u95 = exp(m + h))]
  ow[k < 0.6 * length(cfg$seeds), c("l95", "u95") := NA_real_]
  ow[, pub_in_pi := pub >= l95 & pub <= u95]
  lg <- est[oset, on = .(window, outcome), nomatch = NULL][,
    .(lg = mean(log(hr))),
    keyby = .(window, seed)
  ]
  pl <- cfg$pub[oset, on = .(outcome, window), nomatch = NULL][,
    .(pub_lg = mean(log(pub))),
    keyby = window
  ]
  lg <- lg[pl, on = "window"]
  pooled <- lg[,
    .(
      gm = exp(mean(lg)),
      pub = exp(pub_lg[1])
    ),
    keyby = window
  ]
  sig <- est[
    window == 5L,
    .(n = sum(exp(log(hr) + 1.96 * se) < 1)),
    keyby = seed
  ]
  return(list(ow = ow, pooled = pooled, sig = sig))
}

# The cohort statistics of one seed. Every seed gives the frailty of
# test-negative person-time. The first seed also returns the comparisons with
# the paper, the supplement and SSI.
describe <- function(z, seed, cfg = CFG) {
  d <- z$d
  n <- nrow(d)
  tfirst <- do.call(pmin, c(d[, names(cfg$rate), with = FALSE], na.rm = TRUE))
  # Mean frailty of test-negative person-time per half-year since 1 March
  # 2020, as sums: su is person-days times frailty, st person-days. The
  # test-negative time of a person runs from entry to the recorded infection,
  # the first CVD diagnosis, death or fu, as in analyse().
  stop_neg <- pmin(
    fcoalesce(d$ana, d$fu),
    d$fu,
    fcoalesce(d$tdeath, Inf),
    fcoalesce(tfirst, Inf)
  )
  # The same for person-time after a recorded infection: from the recorded
  # infection to the first CVD diagnosis, death or fu (st_inf, su_inf).
  stop_all <- pmin(d$fu, fcoalesce(d$tdeath, Inf), fcoalesce(tfirst, Inf))
  start_inf <- fcoalesce(d$ana, Inf)
  edges <- c(seq(0, 1036, by = 182), 1036)
  negu <- rbindlist(lapply(seq_len(length(edges) - 1L), function(j) {
    lo <- pmax(edges[j] - d$entry, 0)
    hi <- edges[j + 1L] - d$entry
    t <- pmax(pmin(hi, stop_neg) - lo, 0)
    ti <- pmax(pmin(hi, stop_all) - pmax(lo, start_inf), 0)
    return(data.table(
      half = j,
      lo = edges[j],
      st = sum(t),
      su = sum(t * d$u),
      st_inf = sum(ti),
      su_inf = sum(ti * d$u)
    ))
  }))[, `:=`(seed = seed, u_all = mean(d$u))][]
  if (seed != cfg$seeds[1]) {
    return(list(negu = negu))
  }

  # Stop times and first events, as analyse() sets them. wev is the window of
  # the first event on the analysis clock, 0 for test-negative time.
  sadm <- pmin(d$fu, d$tdeath, na.rm = TRUE)
  stop_t <- pmin(tfirst, sadm, na.rm = TRUE)
  ev <- !is.na(tfirst) & tfirst <= sadm
  ko <- rep(NA_integer_, n)
  for (k in seq_along(cfg$rate)) {
    tk <- d[[names(cfg$rate)[k]]]
    ko[ev & !is.na(tk) & tk == tfirst] <- k
  }
  wev <- fifelse(
    !d$seen | tfirst < d$ana,
    0L,
    findInterval(tfirst - d$ana, EDGE)
  )[ev]
  cal <- (d$entry + tfirst)[ev]

  # Person-years per analysis window, to the stop time (Supp Table 1) and to
  # the stop time or day 740, whichever is first (Supp Table 7).
  lo <- z$s0[, -ncol(z$s0)]
  hi <- z$s0[, -1L]
  pyw <- function(s) {
    w <- pmin(hi, s) - pmin(lo, s)
    return(vapply(0:5, function(k) sum(w[z$win == k]), numeric(1)) / DAY_YR)
  }
  py <- cbind(t1 = pyw(stop_t), t7 = pyw(pmax(pmin(stop_t, 740 - d$entry), 0)))
  rm(lo, hi)
  neg <- wev == 0L
  sim_ev <- tabulate(ko[ev][neg], length(cfg$rate))
  n7 <- sum(neg & cal < 740)
  tn <- c(
    whole = sum(neg) / py[[1, "t1"]],
    before = n7 / py[[1, "t7"]],
    after = (sum(neg) - n7) / (py[[1, "t1"]] - py[[1, "t7"]])
  )

  # Recorded infection by SSI age group, against the mean pps_draw()
  # probability. Unrecorded infection among the never-recorded.
  ag <- findInterval(d$age0, cfg$inf_age$edges)
  rr <- cfg$inf_age$rr[ag]
  pr <- rr * round(n * cfg$p_infected) / sum(rr)
  hid <- !d$seen & !is.na(d$bio)
  inf_age <- data.table(ag = ag, seen = d$seen, pr = pr, hid = hid)[,
    .(n = .N, rec = mean(seen), target = mean(pr), unrec = mean(hid[!seen])),
    keyby = ag
  ]

  # Infections by calendar month, recorded or not. An infection counts only
  # if it comes before death, or at death when it kills.
  ok <- !is.na(d$inf_cal) &
    (is.na(d$tdeath) | (!is.na(d$bio) & d$bio <= d$tdeath))
  i <- which(ok)
  x <- d$inf_cal[i]
  s <- d$seen[i]
  m <- findInterval(x, MON_DAY)
  erik <- x >= ERIK[["lo"]] & x < ERIK[["hi"]]
  return(list(
    negu = negu,
    coh = c(
      n = n,
      age_q = stats::quantile(d$age0, c(0.25, 0.5, 0.75), names = FALSE),
      py = sum(stop_t) / DAY_YR,
      fu_q = stats::quantile(
        stop_t / (DAY_YR / 12),
        c(0.25, 0.5, 0.75),
        names = FALSE
      ),
      rec_n = sum(d$seen & ok),
      cvd_n = sum(ev),
      cvd_pos_n = sum(wev >= 1L),
      erik_n = sum(erik),
      erik_unrec = mean(!s[erik])
    ),
    age_n = as.vector(table(cut(
      d$age0 - 30 / DAY_YR,
      c(-Inf, 40, 65, Inf),
      right = FALSE
    ))),
    py = py,
    sim_ev = sim_ev,
    tn = tn,
    inf_age = inf_age,
    mon = cbind(
      recorded = tabulate(m[s], length(MON)),
      unrecorded = tabulate(m[!s], length(MON))
    )
  ))
}

# The layers both forest plots share: the prediction band, the published CI
# and square, and the true HR, all on the outcome's row y.
forest_layers <- function(labels, bands) {
  ny <- length(labels)
  return(list(
    geom_rect(
      data = bands,
      aes(xmin = 0, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
      fill = "grey93",
      inherit.aes = FALSE
    ),
    geom_vline(xintercept = 1, colour = "grey40", linewidth = 0.5),
    geom_rect(
      aes(
        xmin = l95,
        xmax = u95,
        ymin = y - 0.2,
        ymax = y + 0.2,
        fill = LAB_SIM
      ),
      alpha = 1,
      na.rm = TRUE
    ),
    geom_errorbar(
      aes(y = y, xmin = lo, xmax = hi, linetype = LAB_PUB),
      orientation = "y",
      width = 0,
      linewidth = 0.6,
      colour = "black"
    ),
    geom_point(
      aes(x = pub, y = y, shape = LAB_PUB),
      size = 2,
      colour = "black"
    ),
    geom_point(aes(x = true, y = y, colour = LAB_TRUE), size = 2.3, shape = 16),
    scale_colour_manual(NULL, values = "#d7191c"),
    scale_shape_manual(NULL, values = 15),
    scale_linetype_manual(NULL, values = "solid"),
    guides(
      colour = guide_legend(order = 1),
      fill = guide_legend(order = 2),
      shape = guide_legend(order = 3),
      linetype = guide_legend(order = 3)
    ),
    scale_y_continuous(
      breaks = seq_len(ny),
      labels = labels,
      expand = expansion(add = c(0.5, 0.5))
    )
  ))
}
# The shared theme, with the legend at the bottom.
theme_forest <- function(base_size = 9) {
  return(
    theme_minimal(base_size = base_size) +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(0, "pt"),
        legend.key.spacing.y = unit(1, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.3),
        axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.ticks = element_line(colour = "black", linewidth = 0.6),
        axis.ticks.length = unit(3.5, "pt"),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black")
      )
  )
}
pretty_name <- function(x) {
  return(sub("^(.)", "\\U\\1", gsub("_", " ", x), perl = TRUE))
}

## SSI cases per resident by age group ----
# The md5 of each SSI file MUST equal SSI_MD5.
stopifnot(unname(tools::md5sum(SSI)) == SSI_MD5)
# The sum of "Antal positive" over ISO weeks 2020-W10 to 2022-W52, over the
# median across those weeks of 1e5 * "Antal testede" / "Testede pr. 100.000
# borgere". The cases include reinfections. CFG$inf_age$rr MUST equal them.
x <- fread(SSI[["age"]])
setnames(x, c("week", "age", "test_r", "pos_r", "adm_r", "test_n", "pos_n"))
x <- x[week >= "2020-U10" & week <= "2022-U52"]
ssi_age <- x[,
  .(
    pos_n = sum(pos_n),
    pop = median(1e5 * test_n[test_r > 0] / test_r[test_r > 0])
  ),
  keyby = age
]
ssi_age[, rr := pos_n / pop]
stopifnot(nrow(ssi_age) == 9L, max(abs(ssi_age$rr - CFG$inf_age$rr)) <= 1e-6)

## SSI first infections per month, all 5 regions ----
# CFG$late_n MUST equal them from 11 March 2022 on.
x <- fread(SSI[["month"]], encoding = "Latin-1")
setnames(x, c("date", "region", "infected", "pop", "type", "type_count"))
x <- x[iconv(type, "latin1", "UTF-8") == "2.Bekræftede tilfælde"]
x[, day := as.numeric(as.Date(date) - MON[1])]
x <- x[day >= 0 & day < STUDY]
ssi_mon <- x[, .(n = sum(infected)), keyby = .(m = findInterval(day, MON_DAY))]
late <- x[
  day >= 740,
  sum(infected),
  keyby = .(findInterval(day, CFG$late_edges))
]
stopifnot(
  identical(ssi_mon$m, seq_along(MON)),
  identical(as.numeric(late$V1), CFG$late_n)
)

## SSI wastewater index, national, per week ----
# CFG$ww_w MUST equal it for 2022-W10 to 2022-W52.
x <- fread(SSI[["ww"]])
stopifnot(identical(
  x[week >= "2022-W10" & week <= "2022-W52", rna_mean_faeces],
  CFG$ww_w
))
rm(x, late)


## The simulated cohorts ----
# Each cohort is simulated, analysed and described in the worker that holds
# it, because one cohort needs about 11 GiB. est holds one row per seed,
# outcome and window, and vest the same per variant. The worker then
# simulates the same cohort without infection, with the same baseline scale
# sc, and compares the two. The event times with infection do not depend on
# death at infection. The cohort without infection holds the natural death
# times of the same persons. So e1, the event times of the recorded persons,
# with the death times of d0 gives the risk with infection without death at
# infection. rdn holds the true risk differences after a recorded infection
# and rr the rate differences, per outcome and window, both without death at
# infection.
res <- parallel::mclapply(
  CFG$seeds,
  function(s) {
    z <- sim(s)
    r <- list(
      est = analyse(z, s),
      desc = describe(z, s),
      vest = by_variant(z, s),
      moved = z$moved
    )
    x1 <- first_cvd(z$d)
    j <- which(z$d$seen)
    a <- z$d$ana[j]
    e1 <- z$d[j, c("fu", CFG$outcomes), with = FALSE]
    sc <- z$sc
    rm(z)
    invisible(gc())
    d0 <- sim(s, infect = FALSE, sc = sc)$d
    r$excess <- excess(x1, first_cvd(d0))
    k0 <- risk_after(d0, j, a)
    e1[, tdeath := d0$tdeath[j]]
    r$rdn <- rd_pp(risk_after(e1, seq_along(j), a), k0)[, seed := s][]
    r$rr <- rate_diff(
      rate_after(e1, seq_along(j), a),
      rate_after(d0, j, a)
    )[, seed := s][]
    return(r)
  },
  mc.cores = CFG$n_core,
  mc.preschedule = FALSE
)
stopifnot(!vapply(res, inherits, logical(1), "try-error"))
est <- rbindlist(lapply(res, `[[`, "est"))
vest <- rbindlist(lapply(res, `[[`, "vest"))
negu <- rbindlist(lapply(res, function(r) r$desc$negu))
# Infections moved to another living person because their first person died
# before them, per cohort: recorded and unrecorded.
moved <- do.call(rbind, lapply(res, `[[`, "moved"))
v <- res[[1]]$desc
exc <- lapply(res, `[[`, "excess")
rdn <- rbindlist(lapply(res, `[[`, "rdn"))
rr <- rbindlist(lapply(res, `[[`, "rr"))
rm(res)

## The bias-step cohorts ----
# Each setting of CFG$bias, per seed, analysed with or without censoring at
# the first CVD diagnosis, as cens sets. The main model is the last step, so
# est above gives it.
bgrid <- CFG$bias[, .(seed = CFG$seeds), by = names(CFG$bias)]
bres <- parallel::mclapply(
  seq_len(nrow(bgrid)),
  function(i) {
    g <- bgrid[i]
    cfg <- modifyList(
      CFG,
      list(
        fsd = if (g$frail == 1L) CFG$fsd else 1e-6,
        contam = if (g$unrec == 1L) CFG$contam else 1e-6,
        avoid = g$av,
        p80 = if (g$dai == 1L) CFG$p80 else 0
      )
    )
    r <- analyse(sim(g$seed, cfg), g$seed, cfg, cens = g$cens == 1L)
    return(r[, setting := g$setting][])
  },
  mc.cores = CFG$n_core,
  mc.preschedule = FALSE
)
stopifnot(!vapply(bres, inherits, logical(1), "try-error"))
best <- rbindlist(bres)
rm(bres)

# PART 2 -- ANALYSIS ====
## The outcome set of each window ----
# The outcomes with an estimate in every analysed cohort: the main cohorts and
# the bias-step cohorts. At day 0-1 some outcomes have fewer than 5 events in
# some cohorts and no estimate, so the set holds some of the 12. In the other
# windows it MUST hold all 12. Every average over the outcomes below uses this
# set, except the variant averages, which use vset. true_w holds, per window, the mean log true
# HR and the published geometric mean.
oset <- rbind(
  est[, .(setting = "full", window, outcome, seed)],
  best[, .(setting, window, outcome, seed)]
)[, .N, keyby = .(window, outcome)][
  N == length(CFG$seeds) * (1L + uniqueN(best$setting)),
  .(window, outcome)
]
stopifnot(
  oset[window > 1L, .N] == 4L * length(CFG$outcomes),
  oset[window == 1L, .N] > 0L
)
true_w <- CFG$pub[oset, on = .(outcome, window)][,
  .(
    true_lg = mean(log(CFG$tm[cbind(outcome, colnames(CFG$tm)[window])])),
    pub = exp(mean(log(pub))),
    n_out = .N
  ),
  keyby = window
]

## The estimates against the published ones ----
# The prediction intervals, the pooled rows and the count of outcomes with a
# 95% CI below 1.
cp <- compare_pub(est, oset)
ow <- cp$ow
pooled <- cp$pooled
sig <- cp$sig
n_sig_pub <- CFG$pub[window == 5L, sum(hi < 1)]

## The absolute excess ----
# The counterfactual, per seed: persons with a CVD diagnosis, with infection
# less without it, by infection status with infection. And the risk of a first
# CVD diagnosis in the 12 months after a recorded infection, with and without.
xt <- rbindlist(Map(function(x, s) x$tot[, seed := s], exc, CFG$seeds))
xt <- rbindlist(
  list(xt, xt[, .(grp = "all", excess = sum(excess)), keyby = seed]),
  use.names = TRUE
)
xt <- xt[,
  .(m = mean(excess)),
  keyby = grp
]
xk <- do.call(rbind, lapply(exc, `[[`, "risk"))
# The true risk differences without death at infection, in percentage points,
# per outcome and window: the mean over the seeds. cum runs to the window's
# end, 24 months for 12 months or more.
rdpn <- rdn[, .(n = mean(n), cum = mean(cum)), keyby = .(outcome, window)]
# The rate differences per 100,000 person-years, per outcome and window: the
# mean over the seeds.
rrp <- rr[, .(rd = mean(rd)), keyby = .(outcome, window)]

## The bias steps ----
# lg: the mean log estimate over the outcomes of oset, per setting, window and
# seed. full is the main model.
g <- rbind(
  est[, .(setting = "full", window, outcome, hr, seed)],
  best[, .(setting, window, outcome, hr, seed)]
)[oset, on = .(window, outcome), nomatch = NULL][,
  .(n = .N, lg = mean(log(hr))),
  keyby = .(setting, window, seed)
]
g[true_w, on = "window", n_out := i.n_out]
stopifnot(g[, all(n == n_out)])
# Per step and window, gm is the geometric mean over the seeds. h is the
# half-width of the 95% prediction interval of one study, on the log scale, as
# in compare_pub().
steps <- g[
  setting %in% names(CFG$steps),
  .(
    gm = exp(mean(lg)),
    h = stats::qt(0.975, .N - 1) * sd(lg) * sqrt(1 + 1 / .N),
    n = .N
  ),
  keyby = .(setting, window)
]
stopifnot(all(steps$n == length(CFG$seeds)))
steps[, `:=`(
  step = unname(CFG$steps[setting]),
  k = match(setting, names(CFG$steps))
)]
steps <- rbind(
  true_w[, .(window, gm = exp(true_lg), h = 0, step = "True", k = 0L)],
  steps[order(k, window), .(window, gm, h, step, k)]
)
step_lev <- c("True", unname(CFG$steps))
# The main model with and without death at infection, per window.
dai_tab <- dcast(
  g[
    setting %in% c("full", "nodai"),
    .(gm = exp(mean(lg))),
    keyby = .(window, setting)
  ],
  window ~ setting,
  value.var = "gm"
)
setnames(dai_tab, c("nodai", "full"), c("without", "with"))
setcolorder(dai_tab, c("window", "without", "with"))

## Variant strata ----
# Per variant, the outcomes with an estimate in every cohort (a small stratum
# can have fewer than 5 events). The published geometric mean is taken over
# the same outcomes.
vset <- vest[outcome %in% CFG$outcomes, .N, by = .(variant, outcome)][
  N == length(CFG$seeds),
  .(variant, outcome)
]
var_tab <- vest[vset, on = .(variant, outcome)][,
  .(lg = mean(log(hr))),
  keyby = .(variant, seed)
]
var_tab <- var_tab[, .(sim = exp(mean(lg))), keyby = variant]
pv <- melt(
  CFG$pub_variant,
  id.vars = "outcome",
  variable.name = "variant",
  value.name = "pub",
  variable.factor = FALSE
)
var_tab[
  pv[vset, on = .(variant, outcome)][,
    .(published = exp(mean(log(pub))), n_out = .N),
    keyby = variant
  ],
  on = "variant",
  `:=`(published = i.published, n_out = i.n_out)
]
var_tab[
  pv[, .(published_12 = exp(mean(log(pub)))), keyby = variant],
  on = "variant",
  published_12 := i.published_12
]
var_tab <- var_tab[match(CFG$variant$variant, variant)]

## The seed-1 cohort against the published data ----
N <- sum(PAPER$age_n)
f0 <- function(x) format(round(x), big.mark = ",", scientific = FALSE)
pct <- function(x) sprintf("%.1f%%", 100 * x)
q3 <- function(x) sprintf("%.1f (%.1f-%.1f)", x[2], x[1], x[3])
cv <- as.list(v$coh)
tn_pub <- 1000 *
  c(
    whole = 54247 / PUB_PY[1, "t1"],
    before = PUB_EV7 / PUB_PY[1, "t7"],
    after = (54247 - PUB_EV7) / (PUB_PY[1, "t1"] - PUB_PY[1, "t7"])
  )
cmp <- data.table(
  quantity = c(
    "Persons",
    "Age at start of follow-up, median (IQR), years",
    "Age at first test <40 / 40-64 / 65+",
    "Person-years",
    "Follow-up, median (IQR), months",
    "Recorded positives",
    "Persons with a first CVD diagnosis",
    "... after a recorded positive",
    "Test-negative CVD rate, whole period",
    "Test-negative CVD rate, to 10 March 2022",
    "Test-negative CVD rate, after 10 March 2022",
    "Infected, 1 November 2021 to 15 March 2022",
    "Unrecorded share of those infections"
  ),
  simulated = c(
    f0(cv$n),
    q3(unlist(cv[c("age_q1", "age_q2", "age_q3")])),
    paste(pct(v$age_n / cv$n), collapse = " / "),
    f0(cv$py),
    q3(unlist(cv[c("fu_q1", "fu_q2", "fu_q3")])),
    sprintf("%s (%s)", f0(cv$rec_n), pct(cv$rec_n / cv$n)),
    f0(cv$cvd_n),
    f0(cv$cvd_pos_n),
    sprintf("%.3f", 1000 * v$tn),
    sprintf("%s of all persons", pct(cv$erik_n / cv$n)),
    pct(cv$erik_unrec)
  ),
  published = c(
    f0(N),
    q3(PAPER$age_q),
    paste(pct(PAPER$age_n / N), collapse = " / "),
    f0(PAPER$py),
    q3(PAPER$fu_q),
    sprintf("%s (%s)", f0(CFG$p_infected * N), pct(CFG$p_infected)),
    f0(PAPER$cvd_n),
    f0(PAPER$cvd_pos_n),
    sprintf("%.3f", tn_pub),
    "66% (63-70%) of all blood donors aged 17-72",
    "one third"
  ),
  source = c(
    rep("Paper, Results", 2),
    "Paper, Table 1",
    rep("Paper, Results", 5),
    "Supp Table 1",
    "Supp Table 7",
    "Supp Table 1 minus Table 7",
    rep("Erikstrup et al. 2022", 2)
  )
)
py_cmp <- data.table(
  window = WIN_LAB,
  sim_t1 = v$py[, "t1"],
  pub_t1 = PUB_PY[, "t1"],
  sim_t7 = v$py[, "t7"],
  pub_t7 = PUB_PY[, "t7"]
)
py_cmp[, `:=`(ratio_t1 = sim_t1 / pub_t1, ratio_t7 = sim_t7 / pub_t7)]
rt <- data.table(
  outcome = pretty_name(names(CFG$rate)),
  sim_ev = v$sim_ev,
  pub_ev = round(CFG$rate * 6188401 / 1000),
  sim_rate = 1000 * v$sim_ev / v$py[1, "t1"],
  pub_rate = CFG$rate
)
rt[, ratio := sim_rate / pub_rate]
rt <- rt[order(-pub_rate)]
inf_age <- v$inf_age[, .(
  age = ssi_age$age[ag],
  ssi_rr = ssi_age$rr[ag],
  n,
  recorded = rec,
  target = target,
  z = (rec - target) / sqrt(target * (1 - target) / n),
  unrecorded = unrec
)]
# SSI first infections, scaled so their total over the study period equals
# the paper's 2,698,261 recorded positives.
mon <- data.table(month = format(MON, "%Y-%m"), v$mon, ssi = ssi_mon$n)
mon[, ssi_scaled := ssi * CFG$p_infected * cv$n / sum(ssi)]
mon[, ratio := recorded / ssi_scaled]

# PART 3 -- OUTPUT ====
d12 <- ow[window == 5L]
d12[, y := match(outcome, d12[order(pub), outcome])]
dir.create("results", showWarnings = FALSE)
saveRDS(
  list(
    est = est,
    ow = ow,
    pooled = pooled,
    sig = sig,
    n_sig_pub = n_sig_pub,
    xt = xt,
    xk = xk,
    rdpn = rdpn,
    rrp = rrp,
    negu = negu,
    moved = moved,
    cmp = cmp,
    py_cmp = py_cmp,
    rt = rt,
    inf_age = inf_age,
    mon = mon,
    true_w = true_w,
    steps = steps,
    step_lev = step_lev,
    dai_tab = dai_tab,
    var_tab = var_tab,
    cfg = CFG,
    paper = PAPER,
    pub_py = PUB_PY
  ),
  "results/run.rds"
)

## Figure: 12 months or more, with true hazard ratio and risk difference ----
# Rows sorted by the published HR, the highest at the top. The figures give the
# cumulative risk difference without death at infection, as extra diagnoses
# per 100,000 persons with a recorded infection (pp x 1000). rd_lab() rounds it
# to 1 decimal. At 12 months or more it runs over 0-24 months.
rd_lab <- function(x) {
  x <- round(1000 * x, 1)
  x[x == 0] <- 0
  return(sprintf("%.1f", x))
}
d12[rdpn[window == 5L], on = "outcome", rd := i.cum]
ny <- nrow(d12)
bands <- data.table(y = seq(2L, ny, by = 2L))
q <- ggplot(d12)
q <- q + forest_layers(pretty_name(d12[order(y), outcome]), bands)
q <- q + scale_fill_manual(NULL, values = "#82cab4")
q <- q + scale_x_log10(breaks = c(0.5, 0.71, 1, 1.41))
q <- q +
  coord_cartesian(xlim = c(0.38, 1.6), ylim = c(-1.5, ny + 0.5))
# The direction of the effect: an arrow from 1 to each side, inside the panel
# under the last outcome, with its label beneath.
q <- q +
  annotate(
    "segment",
    x = c(1 / 1.03, 1.03),
    xend = c(0.55, 1.5),
    y = 0.25,
    yend = 0.25,
    arrow = grid::arrow(length = unit(1.5, "mm"), type = "closed"),
    linewidth = 0.5
  )
q <- q +
  annotate(
    "label",
    x = c(1 / 1.03, 1.03),
    hjust = c(1, 0),
    y = -0.8,
    label = c("SARS-CoV-2\nprotective", "SARS-CoV-2\nharmful"),
    lineheight = 0.9,
    size = 2.8,
    fill = "white",
    linewidth = 0,
    label.padding = unit(1, "pt")
  )
q <- q + labs(x = "Hazard ratio at 12 months or more", y = NULL)
q <- q +
  theme_forest() +
  theme(plot.margin = margin(5, 2, 5, 5))
t <- ggplot(
  rbind(
    d12[, .(y, value = sprintf("%.2f", true), x = 0.16)],
    d12[, .(y, value = rd_lab(rd), x = 0.69)]
  ),
  aes(x = x, y = y, label = value)
)
t <- t +
  geom_rect(
    data = bands,
    aes(xmin = 0, xmax = 1, ymin = y - 0.5, ymax = y + 0.5),
    inherit.aes = FALSE,
    fill = "grey93"
  )
t <- t + geom_text(size = 2.9)
t <- t +
  scale_x_continuous(
    position = "top",
    breaks = c(0.16, 0.69),
    labels = c(
      "True HR,\n12 months\nor more",
      "Extra diagnoses\nper 100,000,\n0-24 months,\nwithout death\nat infection"
    ),
    limits = c(0, 1),
    expand = c(0, 0)
  )
t <- t + scale_y_continuous(expand = expansion(add = c(0.5, 0.5)))
t <- t + coord_cartesian(ylim = c(-1.5, ny + 0.5))
t <- t +
  theme_void(base_size = 9) +
  theme(
    axis.text.x.top = element_text(
      face = "bold",
      size = 8,
      margin = margin(b = 3)
    ),
    plot.margin = margin(5, 5, 5, 0)
  )
qq <- patchwork::wrap_plots(q, t, widths = c(1, 0.6)) +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = unit(0, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.margin = margin(0, 0, 0, 0)
  )
dir.create("figures", showWarnings = FALSE)
ggsave(
  "figures/forest_12m.png",
  qq,
  width = 160,
  height = 115,
  units = "mm",
  dpi = 300,
  bg = "white"
)

## Figure: all 5 windows, A4 portrait ----
# Rows sorted by the true HR at 12 months or more, the highest at the top.
# Ties are broken by the published HR. The month panels share one x range.
ord <- ow[window == 5L][order(true, pub), outcome]
pd <- copy(ow)
pd[, y := match(outcome, ord)]
pd[, window_pretty := factor(WL[window], levels = WL)]
blank <- CJ(window = 3:5, x = c(0.45, 1.6))[,
  window_pretty := factor(WL[window], levels = WL)
]
q <- ggplot(pd)
q <- q + forest_layers(pretty_name(ord), bands)
q <- q + geom_blank(data = blank, aes(x = x))
q <- q +
  scale_fill_manual(
    NULL,
    values = "#82cab4",
    labels = paste(strwrap(LAB_SIM, 40), collapse = "\n")
  )
# The month panels keep the breaks 0.71 and 1.41 while their range stays
# below 3.
q <- q +
  scale_x_log10(
    breaks = function(l) {
      if (max(l) > 3) c(0.5, 1, 2, 4, 8, 16, 32) else c(0.5, 0.71, 1, 1.41)
    },
    expand = expansion(mult = 0.04)
  )
q <- q + coord_cartesian(ylim = c(0.5, ny + 0.5))
q <- q +
  facet_wrap(
    ~window_pretty,
    ncol = 2,
    scales = "free_x",
    axes = "all",
    axis.labels = "margins"
  )
q <- q +
  labs(
    x = "Hazard ratio",
    y = NULL,
    title = "Assumed true, biased and published hazard ratios"
  )
q <- q +
  theme_forest() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.605, 0.16),
    legend.justification.inside = c(0, 0.5)
  )
ggsave(
  "figures/forest_outcome_truths.png",
  q,
  width = 210,
  height = 297,
  units = "mm",
  dpi = 200,
  bg = "white"
)

## Figure: from the true hazard ratios to the estimates, step by step ----
# One panel per window. Each estimate carries the 95% prediction interval of
# one study's geometric mean, from the spread over the seeds, as in
# compare_pub().
pd <- copy(steps)
pd[, `:=`(
  lo = gm * exp(-h),
  hi = gm * exp(h),
  step = factor(step, levels = rev(step_lev)),
  win = factor(WL[window], levels = WL),
  what = fifelse(
    k == 0L,
    "True value",
    "Simulated estimate (95% prediction interval, one study)"
  )
)]
pp <- true_w[, .(win = factor(WL[window], levels = WL), pub)]
q <- ggplot(pd, aes(x = gm, y = step))
q <- q + geom_vline(data = pp[pub < 2], aes(xintercept = 1), colour = "grey60")
q <- q +
  geom_vline(
    data = pp,
    aes(xintercept = pub, linetype = "Published (geometric mean)")
  )
q <- q +
  geom_path(
    data = pd[k > 0L],
    aes(group = 1),
    colour = "#1b9e77",
    linewidth = 0.6
  )
q <- q +
  geom_pointrange(
    aes(xmin = lo, xmax = hi, colour = what),
    size = 0.4,
    linewidth = 0.7
  )
q <- q +
  scale_colour_manual(
    NULL,
    values = c(
      "True value" = "#d7191c",
      "Simulated estimate (95% prediction interval, one study)" = "#1b9e77"
    )
  )
q <- q +
  scale_linetype_manual(
    NULL,
    values = c("Published (geometric mean)" = "dashed")
  )
q <- q +
  scale_x_continuous(trans = "log2", breaks = function(l) pretty(l, n = 4))
q <- q +
  facet_wrap(
    ~win,
    nrow = 1,
    scales = "free_x",
    axes = "all",
    axis.labels = "all_x"
  )
q <- q + labs(x = "Hazard ratio, geometric mean over the outcomes", y = NULL)
q <- q +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing.x = unit(14, "pt"),
    strip.text = element_text(colour = "black", face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(3.5, "pt"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )
# What each step shows, as a text column to the right of the panels.
WHAT <- c(
  "True" = "The effect of infection",
  "All biases off" = "Close to the true value",
  "+ unrecorded infections" = "Dilutes toward 1",
  "+ depletion of susceptibles" = "Can push below 1",
  "+ avoidance of infection" = "Can push below 1"
)
tx <- data.table(
  step = factor(names(WHAT), levels = rev(step_lev)),
  what = WHAT
)
qt <- ggplot(tx, aes(x = 0, y = step, label = what))
qt <- qt + geom_text(hjust = 0, size = 3.9)
qt <- qt + scale_x_continuous(limits = c(0, 1), expand = c(0, 0))
qt <- qt + facet_wrap(~"What it shows")
qt <- qt +
  theme_void(base_size = 12) +
  theme(
    strip.text = element_text(
      colour = "black",
      face = "bold",
      hjust = 0,
      margin = margin(b = 5.5)
    )
  )
qq <- patchwork::wrap_plots(q, qt, widths = c(5, 1.25)) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(
  "figures/bias_steps.png",
  qq,
  width = 340,
  height = 95,
  units = "mm",
  dpi = 200,
  bg = "white"
)

## Cohort figures ----
th <- theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(3.5, "pt"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )
# Infections per month, against SSI. The line marks 10 March 2022.
pd <- melt(
  mon[, .(
    date = MON,
    recorded = as.numeric(recorded),
    unrecorded = as.numeric(unrecorded),
    ssi_scaled
  )],
  id.vars = "date"
)
pd[,
  series := factor(
    variable,
    c("recorded", "unrecorded", "ssi_scaled"),
    c(
      "Recorded, simulated",
      "Unrecorded, simulated",
      "SSI first infections, scaled to the cohort"
    )
  )
]
q <- ggplot(
  pd,
  aes(x = date, y = value / 1000, colour = series, linetype = series)
)
q <- q +
  geom_vline(xintercept = MON[1] + 740, colour = "grey50", linewidth = 0.4)
q <- q + geom_line(linewidth = 0.8)
q <- q + geom_point(size = 1.3)
q <- q + scale_colour_manual(NULL, values = c("#1b9e77", "#d95f02", "black"))
q <- q + scale_linetype_manual(NULL, values = c("solid", "solid", "22"))
q <- q + scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y")
q <- q + scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
q <- q + expand_limits(y = 0)
q <- q + labs(x = "Month of infection", y = "Infections per month (thousands)")
q <- q + th
ggsave(
  "figures/cohort_infection_timing.png",
  q,
  width = 200,
  height = 110,
  units = "mm",
  dpi = 300,
  bg = "white"
)
# Person-years per window, simulated against published.
pd <- melt(
  py_cmp[, .(window, sim_t1, pub_t1, sim_t7, pub_t7)],
  id.vars = "window"
)
pd[, `:=`(
  source = factor(
    fifelse(substr(variable, 1, 3) == "sim", "Simulated", "Published"),
    c("Simulated", "Published")
  ),
  table = fifelse(
    grepl("t1$", variable),
    "Supp Table 1: to 31 December 2022",
    "Supp Table 7: to 10 March 2022"
  ),
  window = factor(window, WIN_LAB),
  lab = fifelse(
    value >= 1e5,
    sprintf("%.2fM", value / 1e6),
    sprintf("%.1fk", value / 1e3)
  )
)]
pp <- list()
for (tb in unique(pd$table)) {
  for (neg in c(TRUE, FALSE)) {
    q <- ggplot(
      pd[table == tb & (window == WIN_LAB[1]) == neg],
      aes(x = window, y = value / 1e6, fill = source)
    )
    q <- q + geom_col(position = position_dodge(width = 0.8), width = 0.75)
    q <- q +
      geom_text(
        aes(label = lab),
        position = position_dodge(width = 0.8),
        vjust = -0.4,
        size = 2.6
      )
    q <- q + scale_fill_manual(NULL, values = c("#1b9e77", "grey45"))
    q <- q + scale_y_continuous(expand = expansion(mult = c(0, 0.12)))
    q <- q + expand_limits(y = 0)
    q <- q + labs(x = NULL, y = "Person-years (millions)", title = if (neg) tb)
    q <- q + th + theme(plot.title = element_text(size = 11))
    pp[[length(pp) + 1L]] <- q
  }
}
q <- patchwork::wrap_plots(pp, ncol = 2, widths = c(1, 4)) +
  patchwork::plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = "bottom")
ggsave(
  "figures/cohort_person_time.png",
  q,
  width = 200,
  height = 150,
  units = "mm",
  dpi = 300,
  bg = "white"
)
# Baseline rate per outcome, simulated against published.
pd <- melt(rt[, .(outcome, sim_rate, pub_rate)], id.vars = "outcome")
pd[, `:=`(
  outcome = factor(outcome, rev(rt$outcome)),
  source = factor(
    fifelse(variable == "sim_rate", "Simulated", "Published"),
    c("Simulated", "Published")
  )
)]
q <- ggplot(pd, aes(x = value, y = outcome, colour = source, shape = source))
q <- q + geom_point(size = 2.4, position = position_dodge(width = 0.5))
q <- q + scale_colour_manual(NULL, values = c("#1b9e77", "black"))
q <- q + scale_shape_manual(NULL, values = c(16, 15))
q <- q + scale_x_continuous(expand = expansion(mult = c(0, 0.03)))
q <- q + expand_limits(x = 0)
q <- q +
  labs(x = "First events per 1000 person-years, test-negative time", y = NULL)
q <- q + th
ggsave(
  "figures/cohort_baseline_rates.png",
  q,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 300,
  bg = "white"
)
