# Bias and the long-term cardiovascular risk after SARS-CoV-2 infection

This repository holds the simulation behind a letter to EHJ Open on Boyd et al. 2026 (doi:10.1093/ehjopen/oeag121). The paper is "SARS-CoV-2 infection and long-term risk of cardiovascular and renal morbidity". The paper reports hazard ratios below 1 for most cardiovascular outcomes from 12 months after infection. The simulation shows that an approximation of the study design can produce those estimates when the true effect of infection is harmful.

`Run.R` simulates 10 cohorts of 4,508,489 synthetic persons, one per member of the paper's cardiovascular cohort. It analyses each cohort as the paper does, and compares the estimates with the published ones. It then compares one simulated cohort with the paper, its supplement and national data. Every number below is printed by `Run.R`, is a ratio of printed numbers, or is a value that `Run.R` holds with its source.

## Result

![Assumed true, biased and published hazard ratios](figures/forest_outcome_truths.png)

Each row is one of the 12 cardiovascular outcomes with the highest baseline rate. Each panel is a window since the first positive test.

- **Red dot:** the assumed true hazard ratio. It is at least 1.01 in every outcome and window.
- **Green band:** the 95% prediction interval for the estimate of one study of this size. It comes from the 10 simulated cohorts.
- **Black square and bar:** the published estimate and its 95% CI, from Supplementary Table 1.

The simulated estimates reproduce the published ones:

- 56 of the 58 published estimates that have a band lie inside it (97%). The 2 outside are ischaemic heart disease at 1 to 5 months and heart failure at 6 to 11 months.
- 2 day-1 estimates have no band. Aneurysm or dissection and inflammatory heart disease have an estimate in only 4 and 3 of the 10 cohorts, and a band needs 6.
- At 12 months or more, the geometric mean of the 12 simulated estimates is 0.82 (0.816), against 0.82 (0.825) published. A t-test over the 10 cohorts gives p = 0.53, from Monte Carlo error only.
- At 12 months or more, a mean of 6.1 of the 12 outcomes per cohort have a 95% CI below 1, against 6 published.

At 12 months or more, the true hazard ratio is 1.03 to 1.26 for 7 of the 12 outcomes, and 1.01 for the other 5. The published estimates are below 1 for 11 of them. The difference is bias from four sources in the simulated study:

- unrecorded differences in cardiovascular risk between persons of the same age (shared frailty);
- infections that testing did not record, in the comparison group;
- censoring at the first diagnosis in any of the 15 cardiovascular groups;
- selection into early infection: persons infected before the Omicron wave have a lower baseline cardiovascular risk.

The letter figure is [figures/forest_12m.png](figures/forest_12m.png). It shows the 12-month window only, with rows sorted by the published hazard ratio.

## The model

Each cohort follows its persons from 1 March 2020 to 31 December 2022. The table gives every input, with its source or with the word "assumed".

| Component | How `Run.R` draws it | Source |
|---|---|---|
| Age | Age at the start of follow-up, from a monotone spline. Its knots are the median 34.2, the quartiles 18.5 and 52.3, and the Table 1 band edges at 40 and 65 years. Above the last band, the knots are 72 years (95th percentile), 97 (99.9th) and 105. | Paper, Results and Table 1. The knots above 65 years are assumed. |
| Follow-up | Start of follow-up, from a monotone spline. Its quartiles are day 1036 minus the paper's follow-up quartiles. Its 95th percentile is day 759.7, and its maximum is day 1035. Follow-up stops at the first CVD diagnosis, death or 31 December 2022. | Paper, Results, for the quartiles. The 95th percentile is fitted to the person-years of Supplementary Tables 1 and 7. |
| Frailty | One log-normal multiplier on all 15 CVD hazards, SD 1.65 on the log scale. | Assumed. |
| Outcomes | 15 CVD groups. Each hazard is Gompertz in age, doubling every 7.5 years, times frailty. The level of each outcome is its events per test-negative person-year. One scale factor sets the test-negative first-event rate to 8.766 per 1000 person-years. | Supplementary Table 1: 54,247 events in 6,188,401 test-negative person-years. The age slope is assumed. |
| Calendar step | From 11 March 2022, every outcome hazard is 1.175 times higher. The multiplier and the scale factor are solved against the published test-negative rates of 8.119 per 1000 person-years before that day and 10.870 after. They are solved without an infection effect or unrecorded infections. | Supplementary Table 7, and Table 1 minus Table 7. The step is assumed. |
| Recorded infection | 2,698,261 of 4,508,489 persons (59.8%) have one. The probability is proportional to the SSI confirmed cases per resident of the person's age group. The date is piecewise uniform: 16.4% to 21 December 2021, the Omicron wave to 10 March 2022, and 3.9% after, split by month. | Paper, Results. Statens Serum Institut (SSI) cases by age group and first infections by month, in `data/ssi/`. The dates to 10 March 2022 are fitted with the start of follow-up. |
| Unrecorded infection | 74.5% of persons without a recorded infection get an unrecorded one, so one third of all infections are unrecorded. They follow the same age pattern, with each probability capped at 1. Their dates follow the recorded curve, weighted by (1 - p)/p. p is 0.40 to 10 March 2022 and 0.05 after. Dates after 10 March 2022 follow the weekly SSI wastewater index. | Erikstrup et al. 2022 for the one third. SSI wastewater index in `data/ssi/`. p is assumed. |
| Selection | A person 1 SD of log frailty less frail has exp(0.15) times the odds of a recorded infection before 21 December 2021 (`pre_frailty` 0.15). | Assumed. |
| Effect of infection | Each outcome has five true hazard ratios, one per window: day 0-1, day 2 to 1 month, months 1-5, 6-11 and 12 or more. Recorded and unrecorded infections carry them alike. | Fitted so the simulated estimates reproduce Supplementary Table 1. |
| Death | Gompertz in age: 0.8 per 1000 person-years at the median age, doubling every 8 years, times frailty. Infection kills with probability 0.03 at age 80 and average frailty, more at older ages and higher frailty. | Assumed. |
| Analysis | Person-time split by window since the recorded test, 5-year age band and calendar year, censored at the first CVD diagnosis. One Poisson model per outcome. | Approximates the paper's Cox model. Sex and comorbidity are not simulated. |

## How the simulated cohort matches the known data

`Run.R` compares the seed-1 cohort with the published data. In the Role column, "input" marks a quantity that the model takes from the source or is fitted to. "Check" marks a quantity that no input sets.

| Quantity | Simulated | Published | Source | Role |
|---|---|---|---|---|
| Persons | 4,508,489 | 4,508,489 | Paper, Results | input |
| Age, median (IQR), years | 34.2 (18.5-52.3) | 34.2 (18.5-52.3) | Paper, Results | input |
| Age at first test <40 / 40-64 / 65+ | 57.7% / 32.2% / 10.1% | 57.7% / 32.2% / 10.1% | Paper, Table 1 | input |
| Person-years | 8,813,170 | 8,909,627 | Paper, Results | input |
| Follow-up, median (IQR), months | 25.1 (21.4-27.4) | 25.2 (21.7-27.5) | Paper, Results | input |
| Recorded positives | 2,698,261 (59.8%) | 2,698,261 (59.8%) | Paper, Results | input |
| Recorded share by SSI age group | within 1.7 binomial SE of the target in 9 of 9 groups | the age pattern of SSI cases | SSI | input |
| Test-negative person-years | 6,164,052 | 6,188,401 | Supplementary Table 1 | input |
| Test-negative CVD rate per 1000 person-years, whole period | 8.939 | 8.766 | Supplementary Table 1 | input |
| Same, to 10 March 2022 | 8.233 | 8.119 | Supplementary Table 7 | input |
| Same, after 10 March 2022 | 11.305 | 10.870 | Supplementary Tables 1 and 7 | input |
| Baseline rate per outcome, simulated / published | 0.93 to 1.13 | 1 | Supplementary Table 1 | input |
| Persons with a first CVD diagnosis | 70,614 | 70,885 | Paper, Results | check |
| ... after a recorded positive | 15,516 | 16,475 | Paper, Results | check |
| Person-years, 12 months or more, to 10 March 2022 | 43,481 | 50,715 | Supplementary Table 7 | check |
| Recorded infection before the start of follow-up | 306,571 (6.8%) | 216,898 (4.8%), first test positive | Paper, Table 1 | check |
| Recorded infections in 2020 / 2021 / 2022 | 6.2% / 20.6% / 73.2% | 5.3% / 20.2% / 74.5% | SSI first infections | check |
| First infection, 1 November 2021 to 15 March 2022 | 75.2% of those uninfected on 1 November 2021 | 66% (63-70%) of blood donors aged 17-72 | Erikstrup et al. 2022 | check |
| Unrecorded share of those infections | 26.3% | one third | Erikstrup et al. 2022 | check |

![Recorded and unrecorded infections per month, against SSI](figures/cohort_infection_timing.png)

The SSI line is the national count of first infections per month, times 2,698,261 / 3,151,231. With that factor, the SSI total over the study period equals the paper's recorded positives. So the figure compares the timing, and not the level. The grey vertical line marks 10 March 2022, when widespread testing ended.

![Person-years per exposure window, simulated and published](figures/cohort_person_time.png)

![Baseline CVD rate per outcome in test-negative time, simulated and published](figures/cohort_baseline_rates.png)

The cohort matches the paper in age, person-time, recorded infection and the baseline rate of each outcome. First CVD diagnoses are 0.996 times the published count.

The cohort misses the known data in three places:

- **Follow-up.** The first quartile is 21.4 months, against 21.7. The start of follow-up is fitted to potential follow-up, and death and the first CVD diagnosis shorten the actual follow-up.
- **Recorded infections after March 2022.** From April to December 2022, each month has 0.33 times the scaled SSI count, and March 2022 has 1.38 times. The model puts 3.9% of recorded infections after 10 March 2022. That share is fitted to the person-time of the supplementary tables, and is not taken from SSI.
- **The Erikstrup window.** 75.2% of the persons still uninfected on 1 November 2021 get a first infection before 16 March 2022, above Erikstrup's 63-70%. Of all persons, the share is 66.2%.

Smaller differences remain:

- **CVD diagnoses after a positive.** The cohort has 15,516, which is 0.94 times the published 16,475.
- **Test-negative rate after 10 March 2022.** The rate is 11.305 per 1000 person-years, 1.04 times the published 10.870. The calendar step was solved without an infection effect.
- **Supplementary Table 7.** At 12 months or more before 10 March 2022, the cohort has 0.86 times the published person-years.
- **Infection before the Omicron wave.** The recorded dates before 21 December 2021 are uniform, so the 2020 and 2021 waves are absent.
- **Recorded infection before the start of follow-up.** The cohort has 6.8% of persons, against 4.8% with a positive first test in Table 1.
- **Unrecorded infection by age.** At ages 6-19, every person without a recorded infection gets an unrecorded one, because the cap at 1 binds.

The Erikstrup comparison is weak in three ways. Erikstrup counts healthy blood donors aged 17-72, and the cohort has all ages. The model has first infections only, and a donor can seroconvert on a reinfection. The one-third share applies to the whole study period, but the detection weights move unrecorded infections late. So in the Erikstrup window only 26.3% of infections are unrecorded.

## Assumptions and limits

The result depends on these assumptions:

- **The true hazard ratios are fitted to the published estimates.** The simulation shows that bias can produce the published estimates from harmful true effects. It does not show that the fitted hazard ratios are the true effects of infection.
- **The selection into early infection is assumed.** At `pre_frailty` 0.15, persons with a recorded infection before the Omicron wave have a 22% lower modelled baseline CVD rate than persons infected later. Over the 10 cohorts, the rate ratio is 0.7826 and the ratio of mean frailty is 0.7823. No source measures this in the cohort.
- **One third of infections are unrecorded.** Erikstrup et al. 2022 measured this in blood donors aged 17-72 during one Omicron wave. The model applies it to all ages and the whole study period.
- **The frailty SD of 1.65 is assumed.** So are the detection probabilities 0.40 and 0.05, the death hazard and the death at infection.
- **The calendar step is assumed.** It models a rise in the test-negative rate that the simulation does not explain.
- **The national SSI series are applied to the cohort.** The SSI cases by age group and the wastewater index both include reinfections.
- **The analysis approximates the paper's Cox model.** It uses Poisson models in 5-year age bands and calendar year, without sex or comorbidity.

## How to run it

Run `Rscript Run.R` from the repository root. It needs R 4.6 with data.table, ggplot2, patchwork and knitr.

`Run.R` sets the model in its CFG section, and sources the functions in `R/functions.R`. It checks the md5 of the 3 SSI files in `data/ssi/`. It stops if a value in CFG differs from the SSI file that it comes from. It prints the tables above, draws the 5 figures into `figures/`, and saves the estimates to `results/run.rds`.

With 2 workers the run took 299 seconds on a 20-core Linux machine. The largest R process peaked at 8.66 GiB of memory, and all R processes together at 17.23 GiB. Memory was sampled every 2 seconds, so these are lower bounds. Run it on a machine with at least 20 GB of free memory. On Windows the run uses 1 worker.

## Layout

| Path | Content |
|---|---|
| `Run.R` | The model settings, the published values with their sources, the analysis and the output. |
| `R/functions.R` | The simulation, the analysis and the cohort description. |
| `data/ssi/` | The SSI source files, byte for byte, with their sources and md5s. |
| `figures/` | The 5 figures that `Run.R` draws. |

## Licence

MIT. See `LICENSE`.
