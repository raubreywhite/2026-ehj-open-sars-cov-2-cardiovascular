# Reproduces the letter on Boyd et al. 2026 (Eur Heart J Open,
# doi:10.1093/ehjopen/oeag121). Each of 10 seeds simulates a cohort of
# 4,508,489 synthetic persons under assumed true hazard ratios, and analyses it
# with Poisson models (5-year age bands, calendar year) as an approximation to
# the paper's Cox analysis.
#
# It prints, per outcome at 12 months or more, the assumed true hazard ratio,
# one study's 95% prediction interval and the published estimate. Over all 5
# windows it counts the published hazard ratios inside that interval. It then
# compares the seed-1 cohort with the paper, its supplement, Statens Serum
# Institut (SSI) data in data/ssi/ and Erikstrup et al. 2022. It draws
# figures/forest_12m.png, figures/forest_outcome_truths.png and three
# figures/cohort_*.png, and saves the estimates to results/run.rds.
#
# Needs R 4.6 with data.table, ggplot2, patchwork and knitr.
# With 2 workers the run took 299 seconds on a 20-core Linux machine. The
# largest R process peaked at 8.66 GiB, and all its R processes together at
# 17.23 GiB. Memory was sampled every 2 seconds, so these are lower bounds.
# On Windows the run uses 1 worker.
# Run it from the repository root with: Rscript Run.R

library(data.table)
library(ggplot2)

# CFG ====
CFG <- list()

## Cohort (Boyd et al. 2026, Results and Table 1) ----
CFG$n <- 4508489L
CFG$seeds <- 1:10
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
# The calendar step: from day 740 the hazard of every outcome is step_mult
# times higher. step_mult was solved against the test-negative rates 8.119
# (Supp Table 7) and 10.870 per 1000 person-years after day 740.
CFG$step_day <- 740
CFG$step_mult <- 1.1750

## The simulated world ----
# fsd: SD of log frailty, a shared unrecorded cardiovascular risk.
# p80: the probability that infection kills an 80-year-old of average frailty
# (u = 1); it is multiplied by the person's frailty u.
# detect_early, detect_late: the probability that a test records an infection,
# before and after widespread testing ended on day 740.
# pre_frailty: selection into early infection. A person 1 SD of log frailty
# less frail has exp(pre_frailty) times the odds of an infection before o0.
# All five are assumed. No source measures them.
CFG$fsd <- 1.65
CFG$p80 <- 0.03
CFG$detect_early <- 0.40
CFG$detect_late <- 0.05
CFG$pre_frailty <- 0.15

# contam: the share of never-recorded persons with an unrecorded infection. One
# third of infections were unrecorded (Erikstrup et al. 2022).
CFG$contam <- (CFG$p_infected / (1 - 1 / 3) - CFG$p_infected) /
  (1 - CFG$p_infected)

# The assumed true hazard ratios, one row per outcome, fitted so the simulated
# estimates reproduce the published ones. The windows start at 0, 1, 30,
# 182.62 and 365.25 days since infection.
EDGE <- c(0, 1, 30, 182.62, 365.25)
# The values are the fitted ones, to full precision.
# fmt: skip
CFG$tm <- matrix(c(
  7.262658810455609, 1.1527392137677965, 1.0833576125225581, 1.0833576125225581, 1.0833576125225581,
  10.347612936183138, 1.2719701857808678, 1.01, 1.01, 1.01,
  7.496398730724567, 1.281590790373247, 1.01, 1.01, 1.01,
  11.378584828357422, 1.2690264530984252, 1.0768075162635082, 1.0768075162635082, 1.01,
  4.2683453330071535, 1.6852939405635237, 1.226454289676099, 1.226454289676099, 1.2044506868494067,
  7.291893823922938, 1.2701457122335622, 1.113293310288862, 1.113293310288862, 1.0177552600987922,
  22.753908147597887, 6.475188817032425, 1.5739623304842238, 1.2632101510915585, 1.2632101510915585,
  2.9804472349724365, 1.01, 1.01, 1.01, 1.01,
  7.291893823922938, 1.2701457122335622, 1.113293310288862, 1.113293310288862, 1.0177552600987922,
  16.212648646447782, 1.1054619786967321, 1.01, 1.01, 1.01,
  7.291893823922938, 1.2701457122335622, 1.113293310288862, 1.113293310288862, 1.0177552600987922,
  5.5155226996984394, 1.3702123309931373, 1.1633588707814906, 1.1633588707814906, 1.1633588707814906,
  9.166584483991715, 1.3823345308307158, 1.127829612468253, 1.127829612468253, 1.1184495424019756,
  11.133502702879513, 1.1576971449486315, 1.1576971449486315, 1.1576971449486315, 1.1576971449486315,
  3.217387769877089, 1.0388171766645773, 1.0388171766645773, 1.0388171766645773, 1.031029189424833
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
# test positive plus first test negative, and the first-test positives alone.
PAPER <- list(
  age_q = c(18.5, 34.2, 52.3),
  py = 8909627,
  fu_q = c(21.7, 25.2, 27.5),
  cvd_n = 70885,
  cvd_pos_n = 16475,
  age_n = c(131489 + 2471717, 61057 + 1388941, 24352 + 430933),
  first_pos_n = 131489 + 61057 + 24352
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
# Erikstrup et al. 2022: 66% (95% CI 63-70%) of healthy blood donors aged
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

# FUNCTIONS ====
source("R/functions.R")

# PART 1 -- DATA CREATION ====
stopifnot(unname(tools::md5sum(SSI)) == SSI_MD5)

## SSI cases per resident by age group ----
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
# it, because one cohort needs about 8.6 GiB. est holds one row per seed,
# outcome and window.
res <- parallel::mclapply(
  CFG$seeds,
  function(s) {
    z <- sim(s)
    return(list(est = analyse(z, s), desc = describe(z, s)))
  },
  mc.cores = CFG$n_core,
  mc.preschedule = FALSE
)
stopifnot(!vapply(res, inherits, logical(1), "try-error"))
est <- rbindlist(lapply(res, `[[`, "est"))
sel <- rbindlist(lapply(res, function(r) r$desc$sel))
v <- res[[1]]$desc
rm(res)

# PART 2 -- ANALYSIS ====
# Per outcome and window: mean m and SD s of the log HR over the k seeds. The
# 95% prediction interval of one study's estimate is
# exp(m +- qt(0.975, k - 1) * s * sqrt(1 + 1/k)). Fewer than 6 seeds give no
# interval.
ow <- est[,
  .(m = mean(log(hr)), s = sd(log(hr)), k = .N),
  keyby = .(outcome, window)
]
ow <- ow[CFG$pub, on = .(outcome, window)]
ow[, true := CFG$tm[cbind(outcome, colnames(CFG$tm)[window])]]
ow[, h := qt(0.975, k - 1) * s * sqrt(1 + 1 / k)]
ow[, `:=`(l95 = exp(m - h), u95 = exp(m + h))]
ow[k < 6L, c("l95", "u95") := NA_real_]
ow[, pub_in_pi := pub >= l95 & pub <= u95]

## Pooled over the 12 outcomes ----
# Per seed and window, the mean log HR over the outcomes with an estimate. gm
# is its mean over the seeds, exponentiated, against the geometric mean of the
# published HRs. p is a one-sample t-test over the seeds, so it reflects Monte
# Carlo error only.
lg <- est[, .(lg = mean(log(hr))), keyby = .(window, seed)]
lg <- lg[CFG$pub[, .(pub_lg = mean(log(pub))), keyby = window], on = "window"]
pooled <- lg[,
  .(
    gm = exp(mean(lg)),
    pub = exp(pub_lg[1]),
    p = stats::t.test(lg, mu = pub_lg[1])$p.value
  ),
  keyby = window
]
# Per seed, the outcomes whose 95% CI at 12 months or more lies below 1.
sig <- est[window == 5L, .(n = sum(exp(log(hr) + 1.96 * se) < 1)), keyby = seed]
n_sig_pub <- CFG$pub[window == 5L, sum(hi < 1)]

## Selection into early infection ----
# The recorded pre-wave infected against the later infected, pooled over the
# seeds: the ratio of mean frailty and of modelled baseline CVD rate.
sp <- sel[, .(u = sum(su) / sum(n), rate = sum(sh) / sum(spy)), keyby = grp]
s1 <- sel[seed == CFG$seeds[1]][, .(grp, u = su / n, rate = sh / spy)]
ratio <- function(x) {
  return(x[grp == "pre", c(u = u, rate = rate)] / x[grp == "later", c(u, rate)])
}
sel_ratio <- rbind(all = ratio(sp), seed1 = ratio(s1))

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
    "Recorded infection before the start of follow-up",
    "Persons with a first CVD diagnosis",
    "... after a recorded positive",
    "Test-negative CVD rate, whole period",
    "Test-negative CVD rate, to 10 March 2022",
    "Test-negative CVD rate, after 10 March 2022",
    "First infection, 1 November 2021 to 15 March 2022",
    "Unrecorded share of those infections"
  ),
  simulated = c(
    f0(cv$n),
    q3(unlist(cv[c("age_q1", "age_q2", "age_q3")])),
    paste(pct(v$age_n / cv$n), collapse = " / "),
    f0(cv$py),
    q3(unlist(cv[c("fu_q1", "fu_q2", "fu_q3")])),
    sprintf("%s (%s)", f0(cv$rec_n), pct(cv$rec_n / cv$n)),
    sprintf("%s (%s)", f0(cv$rec_pre_entry), pct(cv$rec_pre_entry / cv$n)),
    f0(cv$cvd_n),
    f0(cv$cvd_pos_n),
    sprintf("%.3f", 1000 * v$tn),
    sprintf(
      "%s of those uninfected on 1 November 2021 (%s of all)",
      pct(cv$erik_n / cv$erik_susc),
      pct(cv$erik_n / cv$n)
    ),
    pct(cv$erik_unrec)
  ),
  published = c(
    f0(N),
    q3(PAPER$age_q),
    paste(pct(PAPER$age_n / N), collapse = " / "),
    f0(PAPER$py),
    q3(PAPER$fu_q),
    sprintf("%s (%s)", f0(CFG$p_infected * N), pct(CFG$p_infected)),
    sprintf("%s (%s)", f0(PAPER$first_pos_n), pct(PAPER$first_pos_n / N)),
    f0(PAPER$cvd_n),
    f0(PAPER$cvd_pos_n),
    sprintf("%.3f", tn_pub),
    "66% (63-70%) of blood donors aged 17-72",
    "one third"
  ),
  source = c(
    rep("Paper, Results", 2),
    "Paper, Table 1",
    rep("Paper, Results", 3),
    "Paper, Table 1 (first test positive)",
    rep("Paper, Results", 2),
    "Supp Table 1",
    "Supp Table 7",
    "Supp Table 1 minus Table 7",
    rep("Erikstrup et al. 2022", 2)
  )
)
age_cmp <- data.table(
  band = c("<40", "40-64", "65+"),
  sim_n = v$age_n,
  pub_n = PAPER$age_n
)
age_cmp[, diff_pp := 100 * (sim_n / cv$n - pub_n / N)]
py_cmp <- data.table(
  window = WIN_LAB,
  sim_t1 = v$py[, "t1"],
  pub_t1 = PUB_PY[, "t1"],
  sim_t7 = v$py[, "t7"],
  pub_t7 = PUB_PY[, "t7"]
)
py_cmp[, `:=`(ratio_t1 = sim_t1 / pub_t1, ratio_t7 = sim_t7 / pub_t7)]
fu_cmp <- data.table(
  quantile = c("Q1", "median", "Q3"),
  sim = unlist(cv[c("fu_q1", "fu_q2", "fu_q3")]),
  pub = PAPER$fu_q
)
fu_cmp[, diff := sim - pub]
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
kab <- function(x, ...) print(knitr::kable(x, format = "pipe", ...))
d12 <- ow[window == 5L]
d12[, y := match(outcome, d12[order(pub), outcome])]
kab(
  d12[
    order(-y),
    .(
      "Outcome, 12 months or more" = pretty_name(outcome),
      "Assumed true HR" = sprintf("%.2f", true),
      "95% prediction interval, one study" = sprintf("%.2f-%.2f", l95, u95),
      "Published HR (95% CI)" = sprintf("%.2f (%.2f-%.2f)", pub, lo, hi)
    )
  ],
  align = "lrrr"
)
cat(sprintf(
  "\nPublished HR inside the 95%% prediction interval: %d of %d outcome-windows with an interval\n",
  sum(ow$pub_in_pi, na.rm = TRUE),
  sum(!is.na(ow$l95))
))
kab(ow[,
  .(n = .N, with_pi = sum(!is.na(l95)), inside = sum(pub_in_pi, na.rm = TRUE)),
  keyby = window
])
cat("\nOutside the interval, or without one:\n")
kab(ow[
  !pub_in_pi | is.na(pub_in_pi),
  .(outcome, window, k, true, l95, u95, pub)
])
cat("\nPooled over the 12 outcomes, per window:\n")
kab(pooled, digits = 4)
cat(sprintf(
  "12 months or more, outcomes with a 95%% CI below 1: mean %.1f of 12 per cohort (seeds: %s), published %d\n",
  mean(sig$n),
  paste(sig$n, collapse = " "),
  n_sig_pub
))
cat(sprintf(
  "Pre-wave / later recorded infected, seeds %d-%d: mean frailty ratio %.4f, modelled baseline CVD rate ratio %.4f (seed 1: %.4f and %.4f)\n",
  min(CFG$seeds),
  max(CFG$seeds),
  sel_ratio["all", "u"],
  sel_ratio["all", "rate"],
  sel_ratio["seed1", "u"],
  sel_ratio["seed1", "rate"]
))

cat("\n## The seed-1 cohort against the published data\n\n")
kab(cmp)
cat("\nAge at the first test (start of follow-up minus 30 days):\n")
kab(age_cmp, digits = 4)
cat(
  "\nFollow-up to the first CVD diagnosis, death or 31 December 2022, months:\n"
)
kab(fu_cmp, digits = 4)
cat("\nPerson-years per window since the recorded positive:\n")
kab(py_cmp, digits = c(0, 0, 0, 0, 0, 4, 4), format.args = list(big.mark = ","))
cat("\nFirst events and rates per 1000 person-years, test-negative time:\n")
kab(rt, digits = 4)
cat("\nRecorded and unrecorded infection by SSI age group at entry:\n")
kab(inf_age, digits = 4)
cat("\nFirst infections per month, from March 2022:\n")
kab(mon[month >= "2022-03"], digits = c(0, 0, 0, 0, 0, 4))
cat(sprintf(
  "Recorded infections in 2020 / 2021 / 2022: simulated %s, SSI %s\n",
  paste(
    pct(mon[, sum(recorded), keyby = substr(month, 1, 4)]$V1 / cv$rec_n),
    collapse = " / "
  ),
  paste(
    pct(mon[, sum(ssi), keyby = substr(month, 1, 4)]$V1 / sum(mon$ssi)),
    collapse = " / "
  )
))
cat(sprintf(
  "SSI first infections, March 2020 to December 2022: %s; scale to the cohort %s / %s = %.4f\n",
  f0(sum(mon$ssi)),
  f0(cv$rec_n),
  f0(sum(mon$ssi)),
  cv$rec_n / sum(mon$ssi)
))
dir.create("results", showWarnings = FALSE)
saveRDS(
  list(est = est, ow = ow, pooled = pooled, sig = sig, sel = sel, desc = v),
  "results/run.rds"
)

## Figure: 12 months or more, with a true hazard ratio column ----
# Rows sorted by the published HR, the highest at the top.
ny <- nrow(d12)
bands <- data.table(y = seq(2L, ny, by = 2L))
q <- ggplot(d12)
q <- q + forest_layers(pretty_name(d12[order(y), outcome]), bands)
q <- q + scale_fill_manual(NULL, values = "#1b9e77")
q <- q + scale_x_log10(breaks = c(0.5, 0.71, 1, 1.41))
q <- q + coord_cartesian(xlim = c(0.38, 1.6), ylim = c(0.5, ny + 0.5))
q <- q + labs(x = "Hazard ratio at 12 months or more (log scale)", y = NULL)
q <- q + theme_forest() + theme(plot.margin = margin(5, 2, 5, 5))
t <- ggplot(
  d12[, .(y, value = sprintf("%.2f", true), x = 0.5)],
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
    breaks = 0.5,
    labels = "True\nhazard ratio",
    limits = c(0, 1),
    expand = c(0, 0)
  )
t <- t + scale_y_continuous(expand = expansion(add = c(0.5, 0.5)))
t <- t + coord_cartesian(ylim = c(0.5, ny + 0.5))
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
qq <- patchwork::wrap_plots(q, t, widths = c(1, 0.19)) +
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
WL <- c(
  "Day 0-1",
  "Day 2 to <1 month",
  "1 to 5 months",
  "6 to 11 months",
  "12 months or more"
)
ord <- ow[window == 5L][order(true, pub), outcome]
ow[, y := match(outcome, ord)]
ow[, window_pretty := factor(WL[window], levels = WL)]
blank <- CJ(window = 3:5, x = c(0.45, 1.6))[,
  window_pretty := factor(WL[window], levels = WL)
]
q <- ggplot(ow)
q <- q + forest_layers(pretty_name(ord), bands)
q <- q + geom_blank(data = blank, aes(x = x))
q <- q +
  scale_fill_manual(
    NULL,
    values = "#1b9e77",
    labels = paste(strwrap(LAB_SIM, 40), collapse = "\n")
  )
# The true HR is written at the right edge of each panel.
q <- q +
  geom_text(
    aes(x = Inf, y = y, label = sprintf("%.2f", true)),
    hjust = 1.15,
    size = 2.6,
    colour = "#d7191c"
  )
q <- q +
  annotate(
    "text",
    x = Inf,
    y = ny + 0.95,
    label = "True",
    hjust = 1.15,
    size = 2.6,
    fontface = "bold",
    colour = "#d7191c"
  )
q <- q +
  scale_x_log10(
    breaks = function(l) {
      if (max(l) > 3) c(0.5, 1, 2, 4, 8, 16, 32) else c(0.5, 0.71, 1, 1.41)
    },
    expand = expansion(mult = c(0.04, 0.3))
  )
q <- q + coord_cartesian(ylim = c(0.5, ny + 1.3), clip = "off")
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
    x = "Hazard ratio (log scale)",
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
