# Compares Norwegian national cardiovascular counts, 2016-2025, with their
# pre-pandemic trend, and with the effect that the simulation's true hazard
# ratios imply for Norway. It reads results/run.rds, so run Run.R first.
#
# Observed: FHI patients per year with at least one contact for the diagnosis
# ("total forekomst"), by age group, in data/fhi/, age-standardised to 2019.
# Expected: a quasi-Poisson trend in the rate per person, per age group, fitted
# on 2016-2019. Implied: (true HR - 1) times the test-negative rate of each
# outcome, per window since infection, times the person-years that Norway
# spends in each window per calendar year. First infections follow the timing
# of the simulated Danish cohort, recorded and unrecorded, to the end of 2022.
# In 2023, 2024 and 2025, a share of persons is reinfected each year, at a
# uniform time in the year: 25%, 50% or 75%. The effects of the infections add
# up: each infection has its own windows, and a reinfection does not end those
# of an earlier one.
#
# A year's effect on these counts lies between that year's implied extra first
# diagnoses and those accumulated since 2020: an extra patient can return in
# later years.
#
# Needs R 4.6 with data.table, ggplot2, knitr, MASS, readxl and csdata.
# Run it from the repository root with: Rscript Norway.R

library(data.table)
library(ggplot2)

# CFG ====
CFG <- list()
CFG$fhi <- "data/fhi/Total forekomst_ Utvalgte diagnoser, antall pasienter.xlsx"
CFG$fhi_md5 <- "6750b08e2fd247555b2f35d061a14475"
CFG$years_fit <- 2016:2019
CFG$years_out <- 2020:2025
CFG$n_draw <- 10000L
CFG$seed <- 1L
# Synthetic persons that carry the infection timing, and the years with one
# reinfection per person.
CFG$n_person <- 200000L
CFG$years_reinf <- 2023:2025
# The reinfection scenarios: the share of persons reinfected in each
# reinfection year.
CFG$p_reinf <- c(0.25, 0.5, 0.75)
# FHI age groups and their ages in csdata.
CFG$age <- list(
  "0-49 år" = 0:49,
  "50-69 år" = 50:69,
  "70-89 år" = 70:89,
  "90 år og eldre" = 90:150
)
# The FHI groups compared, and the simulated outcomes that make each one.
# Stroke is cerebral infarction plus cerebrovascular hemorrhage. Arrhythmias
# is wider than atrial fibrillation and flutter. All CVD is the sum of all 15
# simulated outcomes; a person with two outcomes counts twice.
CFG$map <- list(
  "Hjerte- og karsykdommer (I00-I99)" = c(
    "ischemic_heart_disease",
    "cerebral_infarction",
    "cerebrovascular_hemorrhage",
    "other_cerebrovascular",
    "venous_embolism",
    "arterial_embolism",
    "pulmonary_embolism",
    "aneurysm_dissection",
    "cardiac_arrest",
    "heart_failure",
    "cardiomyopathy",
    "inflammatory_heart_disease",
    "arrhythmias",
    "conduction_disorders",
    "valve_disorders"
  ),
  "Iskemisk hjertesykdom (I20-I25)" = "ischemic_heart_disease",
  "Atrieflimmer og atrieflutter (I48)" = "arrhythmias",
  "Akutt hjerneslag (I61, I63, I64)" = c(
    "cerebral_infarction",
    "cerebrovascular_hemorrhage"
  ),
  "Hjertesvikt (I11.0, I13.0, I13.2, I42.0, I43, I50)" = "heart_failure"
)
CFG$label <- c(
  "All CVD (I00-I99)",
  "Ischemic heart disease (I20-I25)",
  "Atrial fibrillation and flutter (I48)",
  "Stroke (I61, I63, I64)",
  "Heart failure"
)
# The window edges since infection, in days, as in Run.R. The last window is
# open.
EDGE <- c(0, 1, 30, 182.62, 365.25, Inf)

# PART 1 -- DATA CREATION ====
stopifnot(unname(tools::md5sum(CFG$fhi)) == CFG$fhi_md5)
x <- setDT(readxl::read_excel(
  CFG$fhi,
  col_names = FALSE,
  .name_repair = "minimal"
))
setnames(x, c("title", "age", "source", "dx", paste0("y", 2016:2025)))
x <- x[!is.na(dx) & !is.na(y2016)]
x[, age := age[cummax(seq_len(.N) * !is.na(age))]]
d <- melt(
  x[age %in% names(CFG$age) & dx %in% names(CFG$map)],
  id.vars = c("age", "dx"),
  measure.vars = paste0("y", 2016:2025),
  variable.name = "year",
  value.name = "n"
)
d[, `:=`(year = as.integer(sub("y", "", year)), n = as.numeric(n))]
# The age groups sum to the FHI total for every diagnosis and year.
tot <- melt(
  x[age == "Alle aldre" & dx %in% names(CFG$map)],
  id.vars = "dx",
  measure.vars = paste0("y", 2016:2025),
  variable.name = "year",
  value.name = "n"
)
tot[, `:=`(year = as.integer(sub("y", "", year)), n = as.numeric(n))]
stopifnot(all.equal(
  d[, .(n = sum(n)), keyby = .(dx, year)],
  tot[, .(dx, year, n)][order(dx, year)],
  check.attributes = FALSE
))

# Population: the mean of 1 January of the year and of the next year.
pop <- csdata::nor_population_by_age_cats(cats = unname(CFG$age))
pop <- pop[location_code == "nation_nor" & age != "total"]
stopifnot(identical(sort(unique(pop$age)), c("000_049", "050_069", "070_089", "090_150")))
pop[, age := names(CFG$age)[match(age, sort(unique(age)))]]
pop <- pop[order(age, calyear)][,
  .(year = calyear, pop_n = (pop_jan1_n + shift(pop_jan1_n, -1L)) / 2),
  by = age
]
d[pop, on = .(age, year), pop_n := i.pop_n]
stopifnot(!anyNA(d$pop_n))

# Infection timing of the simulated cohort, seed 1: the share of persons with a
# first infection, recorded or unrecorded, per month.
r <- readRDS("results/run.rds")
mon <- r$mon[, .(
  start = as.numeric(as.Date(paste0(month, "-01"))),
  share = (recorded + unrecorded) / r$cfg$n
)]
mon[, end := c(start[-1L], as.numeric(as.Date("2023-01-01")))]

# PART 2 -- ANALYSIS ====
## Implied excess ----
# Infection days of the synthetic persons: a first infection drawn from mon,
# uniform within its month, or none; then one uniform in each reinfection
# year. Each person is reinfected in every reinfection year here; the share p
# below scales that down.
set.seed(CFG$seed)
np <- CFG$n_person
k <- sample(
  c(seq_len(nrow(mon)), NA),
  np,
  replace = TRUE,
  prob = c(mon$share, 1 - sum(mon$share))
)
t0 <- mon$start[k] + runif(np) * (mon$end[k] - mon$start[k])
yd <- function(y) as.numeric(as.Date(paste0(y, "-01-01")))
tr <- sapply(CFG$years_reinf, function(y) yd(y) + runif(np) * (yd(y + 1L) - yd(y)))
tt <- cbind(t0, tr)
# Person-years per person in each window since each infection, per calendar
# year, summed over the infections: py_first after the first infection, and
# py_reinf after the reinfections. The effects of the infections add up: a
# reinfection does not end the windows of an earlier infection.
pyj <- function(j, w, y) {
  s <- pmax(tt[, j] + EDGE[w], yd(y))
  e <- pmin(tt[, j] + EDGE[w + 1L], yd(y + 1L))
  return(sum(pmax(e - s, 0), na.rm = TRUE) / np / 365.25)
}
pyw <- CJ(window = 1:5, year = CFG$years_out)
pyw[, py_first := mapply(function(w, y) pyj(1L, w, y), window, year)]
pyw[, py_reinf := mapply(
  function(w, y) sum(sapply(2:ncol(tt), pyj, w = w, y = y)),
  window,
  year
)]
# The extra rate per 100,000 person-years, per outcome and window: (true HR -
# 1) times the test-negative rate of the outcome in the Danish cohort
# (Supplementary Table 1).
xr <- data.table(
  outcome = rep(rownames(r$cfg$tm), 5L),
  window = rep(1:5, each = nrow(r$cfg$tm)),
  hr = as.vector(r$cfg$tm)
)
xr[, rd := (hr - 1) * r$cfg$rate[outcome] * 100]
# The implied extra patients per 100,000 persons per year, when a share p of
# persons is reinfected in each reinfection year. lo counts each year's extra
# first diagnoses; hi adds those of all earlier years from 2020, since an extra
# patient can return in later years.
implied <- function(p) {
  return(rbindlist(lapply(names(CFG$map), function(g) {
    z <- xr[outcome %in% CFG$map[[g]], .(rd = sum(rd)), keyby = window]
    z <- pyw[z, on = "window"]
    z <- z[, .(x = sum((py_first + p * py_reinf) * rd)), keyby = year]
    z[, xc := cumsum(x)]
    return(z[, .(dx = g, year, lo = pmin(x, xc), hi = pmax(x, xc))])
  })))
}

## Observed against the pre-pandemic trend ----
# Rates are age-standardised to the Norwegian population of 2019 over the 4
# FHI age groups, per 100,000. A quasi-Poisson model per diagnosis fits a
# log-linear trend in the rate per person, per age group, on 2016-2019. Its
# 95% prediction interval draws the coefficients from their estimated
# distribution and adds overdispersed Poisson noise (variance phi * mu).
w <- d[year == 2019 & dx == dx[1L], .(age, w = pop_n / sum(pop_n))]
d[w, on = "age", w := i.w]
set.seed(CFG$seed)
res <- rbindlist(lapply(names(CFG$map), function(g) {
  dg <- d[dx == g][order(year, age)]
  fit <- stats::glm(
    n ~ age * I(year - 2019) + offset(log(pop_n)),
    family = stats::quasipoisson(),
    data = dg[year %in% CFG$years_fit]
  )
  phi <- summary(fit)$dispersion
  mm <- stats::model.matrix(~ age * I(year - 2019), dg)
  b <- MASS::mvrnorm(CFG$n_draw, stats::coef(fit), stats::vcov(fit))
  mu <- exp(mm %*% t(b) + log(dg$pop_n))
  y <- matrix(stats::rnorm(length(mu), mu, sqrt(phi * mu)), nrow(mu))
  sw <- 1e5 * dg$w / dg$pop_n
  ys <- rowsum(y * sw, dg$year)
  ex <- exp(mm %*% stats::coef(fit) + log(dg$pop_n))
  return(data.table(
    dx = g,
    year = sort(unique(dg$year)),
    observed = as.vector(rowsum(dg$n * sw, dg$year)),
    expected = as.vector(rowsum(ex * sw, dg$year)),
    expected_p025 = apply(ys, 1, stats::quantile, 0.025),
    expected_p975 = apply(ys, 1, stats::quantile, 0.975),
    phi = phi
  ))
}))
res[, dx_pretty := factor(CFG$label[match(dx, names(CFG$map))], levels = CFG$label)]
## The reinfection scenarios ----
imp <- rbindlist(lapply(CFG$p_reinf, function(p) implied(p)[, p := p]))
imp[res, on = .(dx, year), `:=`(lo = i.expected + lo, hi = i.expected + hi, dx_pretty = i.dx_pretty)]
imp[, p_pretty := factor(
  sprintf("%d%% reinfected per year", round(100 * p)),
  levels = sprintf("%d%% reinfected per year", round(100 * CFG$p_reinf))
)]

# PART 3 -- OUTPUT ====
tab <- dcast(imp, dx_pretty + year ~ p, value.var = c("lo", "hi"))
tab <- res[, .(dx_pretty, year, observed, expected, expected_p025, expected_p975)][tab, on = .(dx_pretty, year)]
print(knitr::kable(tab[order(dx_pretty, year)], format = "pipe", digits = 0))
dir.create("results", showWarnings = FALSE)
saveRDS(list(res = res, imp = imp, tab = tab, pyw = pyw, d = d, fhi = x, cfg = CFG), "results/norway.rds")

rp <- res[, .(p_pretty = levels(imp$p_pretty)), by = names(res)]
q <- ggplot(rp, aes(x = year))
q <- q + geom_ribbon(aes(ymin = expected_p025, ymax = expected_p975, fill = "Pre-pandemic trend (95% prediction interval)"))
q <- q + geom_line(aes(y = expected), linetype = "dashed")
q <- q + geom_ribbon(data = imp, aes(ymin = lo, ymax = hi, fill = "Trend plus the effect of the true hazard ratios"), alpha = 0.7)
q <- q + geom_point(aes(y = observed, shape = "Observed"), size = 2)
q <- q + scale_fill_manual(NULL, values = c("grey85", "#1b9e77"))
q <- q + scale_shape_manual(NULL, values = 16)
q <- q + scale_x_continuous(breaks = seq(2016, 2025, by = 3))
q <- q + facet_grid(dx_pretty ~ p_pretty, scales = "free_y", axes = "all", axis.labels = "all_y", labeller = labeller(dx_pretty = label_wrap_gen(18)))
q <- q + labs(x = NULL, y = "Patients per 100,000")
q <- q +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(3.5, "pt"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )
ggsave(
  "figures/norway_excess.png",
  q,
  width = 240,
  height = 280,
  units = "mm",
  dpi = 200,
  bg = "white"
)
