# The functions of Run.R: the simulation sim(), the analysis analyse(), the
# cohort description describe() and the forest plot layers. Run.R sources
# this file after it sets CFG and the constants that the functions read.

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

# 16-point Gauss-Legendre nodes t and weights g on (0, 1), by the Golub-Welsch
# method. The weights sum to 1.
GL16 <- local({
  k <- seq_len(15L)
  J <- matrix(0, 16L, 16L)
  J[cbind(k, k + 1L)] <- J[cbind(k + 1L, k)] <- k / sqrt(4 * k^2 - 1)
  e <- eigen(J, symmetric = TRUE)
  list(t = (e$values + 1) / 2, g = e$vectors[1L, ]^2)
})

# One simulated cohort. Follow-up is cut into 7 segments per person: before
# infection, one per exposure window, and one more at the calendar step. The
# hazard of each outcome is Gompertz in age, times frailty, times the calendar
# step, times the outcome's hazard ratio in the segment's window. One Exp(1)
# threshold per person and outcome inverts in closed form to the event time.
sim <- function(seed, cfg = CFG) {
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

  # Recorded infections. Exactly round(n * p_infected) persons are infected.
  # pps_draw() draws each with a probability proportional to the rr of their
  # age group.
  rr <- cfg$inf_age$rr[findInterval(d$age0, cfg$inf_age$edges)]
  hit <- rep(FALSE, n)
  hit[pps_draw(rr, round(n * cfg$p_infected))] <- TRUE
  # The probability of a date before o0 is a logistic tilt on standardised log
  # frailty. uniroot() sets the intercept so that share stays at f_pre. The
  # other dates come from the rest of the curve.
  d[, inf_cal := NA_real_]
  zf <- (log(d$u) + cfg$fsd^2 / 2) / cfg$fsd
  a <- stats::uniroot(
    function(x) mean(stats::plogis(x - cfg$pre_frailty * zf)) - cfg$f_pre,
    interval = c(-20, 20),
    extendInt = "yes"
  )$root
  pre <- runif(n) < stats::plogis(a - cfg$pre_frailty * zf)
  d[hit & pre, inf_cal := runif(.N, 90, cfg$o0)]
  d[
    hit & !pre,
    inf_cal := qcurve(runif(.N, cfg$f_pre, 1), cfg$inf_cm, cfg$inf_x)
  ]
  # bio: infection on the biological clock, in days since entry. ana: the
  # analysis clock, which starts a recorded infection before entry at entry.
  d[, bio := pmax(inf_cal - entry, -30)]
  d[is.na(inf_cal) | bio >= fu, bio := NA_real_]
  d[, ana := pmax(bio, 0)]
  d[, seen := !is.na(ana)]

  # Unrecorded infections. The fitted curve describes recorded infections, so
  # the missed ones follow it weighted by (1 - p) / p, with p the probability
  # of detection. They go to never-recorded persons by pps_draw(), with the
  # same rr by age group.
  cand <- which(!d$seen)
  k <- round(length(cand) * cfg$contam)
  m <- 4L * k
  dt0 <- qcurve(runif(m), cfg$inf_cm, cfg$inf_x)
  pdet <- ifelse(dt0 < 740, cfg$detect_early, cfg$detect_late)
  dt0 <- dt0[sample(seq_len(m), k, replace = TRUE, prob = (1 - pdet) / pdet)]
  pick <- cand[pps_draw(rr[cand], k)]
  # A hidden date at or after day 740 is redrawn from the wastewater shape, on
  # its own stream. The main stream is restored, so later draws do not change.
  rs <- .GlobalEnv$.Random.seed
  set.seed(seed + 6000000L)
  late <- which(dt0 >= 740)
  dt0[late] <- qcurve(runif(length(late)), cfg$ww_cm, cfg$ww_edges)
  .GlobalEnv$.Random.seed <- rs
  set(d, pick, "inf_cal", dt0)
  set(d, pick, "bio", dt0 - d$entry[pick])
  d[!is.na(bio) & bio >= fu, bio := NA_real_]

  # Segment edges in days since entry. A recorded infection splits follow-up
  # on the analysis clock, an unrecorded one on the biological clock. The
  # calendar step adds one more edge. A bubble sort puts the edges in order.
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
  s0 <- cbind(s0, pmin(pmax(cfg$step_day - d$entry, 0), d$fu))
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

  # Death, Gompertz in age and proportional to frailty, then death at
  # infection with a probability that rises with age and frailty. The draws
  # use their own stream, and the main stream is restored after them.
  rs <- .GlobalEnv$.Random.seed
  Ad <- (DAY_YR / bd) * exp(bd * (d$age0 - cfg$age_median))
  kd <- (0.0008 / DAY_YR) * d$u
  Hd <- kd * (Ad * exp(bd * d$fu / DAY_YR) - Ad)
  set.seed(seed + 2000000L)
  ed <- rexp(n)
  d[, tdeath := NA_real_]
  die <- ed < Hd
  d[die, tdeath := (DAY_YR / bd) * log((Ad[die] + ed[die] / kd[die]) / Ad[die])]
  pcd <- cfg$p80 *
    exp(bd * (d$age0 - cfg$age_median)) /
    exp(bd * (80 - cfg$age_median)) *
    d$u
  kill <- !is.na(d$bio) & runif(n) < pmin(pcd, 0.95)
  d[kill, tdeath := pmin(tdeath, pmax(bio, 0), na.rm = TRUE)]
  .GlobalEnv$.Random.seed <- rs

  # sc scales the baseline so the rate over test-negative person-time matches
  # the published 8.766 per 1000 person-years, on a sample of 200,000 persons.
  # Test-negative time runs from entry to the recorded infection, death or fu.
  # The hazard has no infection multiplier, and it is step_mult times higher
  # after the step. The expected first events equal the rate times the expected
  # time at risk, which also stops at the first event. GL16 integrates it.
  sset <- sample(seq_len(n), 200000L)
  tcal <- pmin(ifelse(is.na(d$ana), d$fu, d$ana), d$tdeath, na.rm = TRUE)
  ltot <- sum(lam)
  k0 <- ltot *
    d$u[sset] *
    (DAY_YR / b) *
    exp(b * (d$age0[sset] - cfg$age_median))
  rt <- cfg$rate_true / 1000 / DAY_YR
  beta <- b / DAY_YR
  s_end <- tcal[sset]
  s_day <- pmin(pmax(cfg$step_day - d$entry[sset], 0), s_end)
  # The probability W of a first event by span, and the expected time at risk
  # to span, from a start with hazard constant K.
  tar <- function(K, span) {
    W <- -expm1(-K * expm1(beta * span))
    tm <- W *
      as.vector((1 / (K - log1p(-outer(W, GL16$t)))) %*% GL16$g) /
      beta
    return(list(W = W, tm = tm))
  }
  calib <- function(x) {
    K <- x * k0
    p1 <- tar(K, s_day)
    s1 <- exp(-K * expm1(beta * s_day))
    p2 <- tar(cfg$step_mult * K * exp(beta * s_day), s_end - s_day)
    return(
      sum(p1$W) + sum(s1 * p2$W) - rt * (sum(p1$tm) + sum(s1 * p2$tm))
    )
  }
  sc <- stats::uniroot(calib, interval = c(0.1, 100), extendInt = "yes")$root
  cmult <- matrix(
    ifelse(as.vector(d$entry + mid) >= cfg$step_day, cfg$step_mult, 1),
    nrow(mid),
    ncol(mid)
  )

  # Event time of each outcome: where the cumulative hazard reaches E. Row k
  # of cfg$tm gives outcome k its hazard ratio per window.
  set.seed(seed + 1000000L)
  E <- matrix(rexp(n * length(lam)), nrow = n)
  for (k in seq_along(lam)) {
    kk <- d$u *
      lam[k] *
      sc *
      matrix(c(1, cfg$tm[k, ])[wb + 1L], nrow(mid)) *
      cmult
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
  return(list(d = d, s0 = s0, win = win, sc = sc))
}

# The paper's analysis. Follow-up stops at the first cardiovascular event of
# any type. A Poisson model per outcome adjusts the exposure window for 5-year
# age band and calendar year.
analyse <- function(z, seed, cfg = CFG) {
  d <- z$d
  d[, tfirst := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = names(cfg$rate)]
  sadm <- pmin(d$fu, d$tdeath, na.rm = TRUE)
  stop_t <- pmin(d$tfirst, sadm, na.rm = TRUE)

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
  p[, py := (q1 - q0) / DAY_YR]
  p[, last := !duplicated(row, fromLast = TRUE)]

  retval <- list()
  for (o in cfg$outcomes) {
    tk <- d[[o]]
    isev <- !is.na(tk) & !is.na(d$tfirst) & tk == d$tfirst & tk <= sadm
    p[, ev := as.integer(isev[row] & abs(q1 - stop_t[row]) < 1e-8 & last)]
    a <- p[, .(py = sum(py), ev = sum(ev)), keyby = .(win, ab, cp)]
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

# The cohort statistics of one seed. Every seed gives the sums behind the
# selection ratio. Recorded infections split into pre-wave, a date before o0,
# and later. A recorded date clamped at 30 days before entry is pre-wave when
# that bound is at or before o0, and otherwise unknown. The modelled baseline
# rate of a person is the expected count of CVD events over potential
# follow-up without infection, per person-year, up to a constant. The first
# seed also returns the comparisons with the paper, the supplement and SSI.
describe <- function(z, seed, cfg = CFG) {
  d <- z$d
  n <- nrow(d)
  b <- log(2) / 7.5
  clip <- d$seen & d$bio == -30
  dt <- d$entry + d$bio
  grp <- fcase(
    !d$seen                             ,
    NA_character_                       ,
    dt < cfg$o0 | (clip & dt <= cfg$o0) ,
    "pre"                               ,
    clip                                ,
    "unknown"                           ,
    default = "later"
  )
  h0 <- z$sc *
    d$u *
    exp(b * (d$age0 - cfg$age_median)) *
    (exp(b * d$fu / DAY_YR) - 1)
  sel <- data.table(grp = grp, u = d$u, h0 = h0, py = d$fu / DAY_YR)[
    !is.na(grp),
    .(n = .N, su = sum(u), sh = sum(h0), spy = sum(py)),
    keyby = grp
  ][, seed := seed][]
  if (seed != cfg$seeds[1]) {
    return(list(sel = sel))
  }

  # Stop times and first events, as analyse() sets them. wev is the window of
  # the first event on the analysis clock, 0 for test-negative time.
  tfirst <- do.call(pmin, c(d[, names(cfg$rate), with = FALSE], na.rm = TRUE))
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

  # Infections by calendar month, recorded or not.
  i <- which(!is.na(d$inf_cal))
  x <- d$inf_cal[i]
  s <- d$seen[i]
  m <- findInterval(x, MON_DAY)
  erik <- x >= ERIK[["lo"]] & x < ERIK[["hi"]]
  return(list(
    sel = sel,
    coh = c(
      n = n,
      age_q = stats::quantile(d$age0, c(0.25, 0.5, 0.75), names = FALSE),
      py = sum(stop_t) / DAY_YR,
      fu_q = stats::quantile(
        stop_t / (DAY_YR / 12),
        c(0.25, 0.5, 0.75),
        names = FALSE
      ),
      rec_n = sum(d$seen),
      rec_pre_entry = sum(d$seen & d$inf_cal < d$entry),
      cvd_n = sum(ev),
      cvd_pos_n = sum(wev >= 1L),
      erik_n = sum(erik),
      erik_susc = n - sum(x < ERIK[["lo"]]),
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
      alpha = 0.55,
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
