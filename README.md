

# Bias and the long-term cardiovascular risk after SARS-CoV-2 infection

This repository holds the simulation behind a letter to EHJ Open on Boyd et al. 2026 (doi:10.1093/ehjopen/oeag121). The paper is “SARS-CoV-2 infection and long-term risk of cardiovascular and renal morbidity”. The paper reports hazard ratios below 1 for most cardiovascular outcomes from 12 months after infection. The simulation shows that an approximation of the study design can produce those estimates when the true effect of infection is harmful.

`Run.R` simulates 10 cohorts of 4,508,489 synthetic persons, one per member of the paper’s cardiovascular cohort. It analyses each cohort as the paper does, and compares the estimates with the published ones. It then compares one simulated cohort with the paper, its supplement and national data. `README.qmd` reads every simulation result in this README from `results/run.rds`.

## Result

At 12 months or more, the assumed true hazard ratio is at least 1.01 for all 12 outcomes. The published estimates are below 1 for 11 of them. The simulated study, analysed as the paper does, reproduces the published estimates: 55 of the 58 published estimates with a prediction interval lie inside it.

![Assumed true, biased and published hazard ratios](figures/forest_outcome_truths.png)

Each row is one of the 12 cardiovascular outcomes with the highest baseline rate. Each panel is a window since the first positive test.

- **Red dot:** the assumed true hazard ratio. It is the effect of infection on the hazard of each person, given age and frailty. It is not a population (marginal) hazard ratio.
- **Green band:** the 95% prediction interval (PI) for the estimate of one study of this size. It comes from the 10 simulated cohorts. A cohort has an estimate only if the window has at least 5 events, and a band needs estimates from at least 6 cohorts.
- **Black square and bar:** the published estimate and its 95% CI, from Supplementary Table 1.

### Values per outcome

The tables give the values in the figure, one table per outcome. The last column gives the extra diagnoses that infection causes, counted to the end of the window. So it adds up the effect of the true hazard ratios in that window and in all earlier windows.

The extra diagnoses compare the same persons with and without infection. They are per 100,000 persons with a recorded infection. The count starts at the recorded infection, or at entry for an infection before entry. It ends 1 day, 30 days, 6 months, 12 months and 24 months later, one end per window. In this comparison infection does not kill: each person has the same death time with and without infection. A death before the diagnosis counts as no diagnosis. A diagnosis of another outcome does not stop the count. The base is the persons whose study period reaches the end of the window. The value is the mean over the 10 cohorts.

#### Pulmonary embolism

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 22.75 | 11.03-21.70 | 16.10 (11.00-23.50) | 1.5 |
| Day 2 to \<1 month | 6.48 | 3.83-5.02 | 4.14 (3.57-4.79) | 13.0 |
| 1 to 5 months | 1.57 | 1.06-1.35 | 1.21 (1.08-1.36) | 19.8 |
| 6 to 11 months | 1.26 | 0.81-1.08 | 0.89 (0.78-1.02) | 20.2 |
| 12 months or more | 1.26 | 0.49-1.11 | 0.77 (0.61-0.98) | 20.4 |

#### Venous embolism

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 4.27 | 1.59-5.44 | 2.67 (1.48-4.84) | 0.4 |
| Day 2 to \<1 month | 1.69 | 1.20-1.69 | 1.42 (1.22-1.66) | 3.3 |
| 1 to 5 months | 1.23 | 1.04-1.30 | 1.12 (1.03-1.21) | 9.0 |
| 6 to 11 months | 1.23 | 1.03-1.20 | 1.09 (1.00-1.18) | 14.3 |
| 12 months or more | 1.20 | 0.66-1.12 | 0.90 (0.77-1.05) | 18.9 |

#### Inflammatory heart disease

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 5.52 | none | 4.90 (2.04-11.80) | 0.1 |
| Day 2 to \<1 month | 1.37 | 0.72-1.93 | 1.08 (0.76-1.54) | 0.4 |
| 1 to 5 months | 1.16 | 0.84-1.32 | 0.99 (0.84-1.18) | 1.3 |
| 6 to 11 months | 1.16 | 0.91-1.34 | 1.05 (0.89-1.23) | 2.4 |
| 12 months or more | 1.16 | 0.58-1.17 | 0.91 (0.67-1.23) | 2.7 |

#### Conduction disorders

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 11.13 | 5.40-17.86 | 8.85 (5.01-15.60) | 0.5 |
| Day 2 to \<1 month | 1.16 | 0.69-1.65 | 0.76 (0.53-1.09) | 0.8 |
| 1 to 5 months | 1.16 | 0.89-1.41 | 1.08 (0.94-1.24) | 2.2 |
| 6 to 11 months | 1.16 | 0.86-1.30 | 1.13 (0.98-1.30) | 3.7 |
| 12 months or more | 1.16 | 0.60-1.37 | 1.13 (0.88-1.44) | 5.2 |

#### Arrhythmias

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 9.17 | 5.56-9.95 | 7.67 (6.17-9.53) | 3.1 |
| Day 2 to \<1 month | 1.38 | 1.10-1.39 | 1.19 (1.07-1.32) | 7.5 |
| 1 to 5 months | 1.13 | 0.97-1.18 | 1.07 (1.02-1.12) | 16.0 |
| 6 to 11 months | 1.13 | 0.97-1.09 | 1.05 (1.00-1.11) | 23.0 |
| 12 months or more | 1.12 | 0.79-0.92 | 0.90 (0.82-0.99) | 28.8 |

#### Ischemic heart disease

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 7.26 | 4.54-7.77 | 5.84 (4.21-8.11) | 1.5 |
| Day 2 to \<1 month | 1.15 | 0.94-1.26 | 1.01 (0.87-1.17) | 2.4 |
| 1 to 5 months | 1.08 | 0.99-1.11 | 0.93 (0.87-1.00) | 6.0 |
| 6 to 11 months | 1.08 | 0.90-1.10 | 1.03 (0.93-1.14) | 9.4 |
| 12 months or more | 1.08 | 0.67-1.07 | 0.89 (0.78-1.00) | 14.1 |

#### Valve disorders

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 3.22 | 1.12-6.04 | 2.51 (1.20-5.28) | 0.3 |
| Day 2 to \<1 month | 1.04 | 0.85-1.23 | 1.03 (0.83-1.28) | 0.4 |
| 1 to 5 months | 1.04 | 0.87-1.17 | 0.98 (0.89-1.08) | 1.2 |
| 6 to 11 months | 1.04 | 0.89-1.12 | 1.03 (0.93-1.14) | 1.8 |
| 12 months or more | 1.03 | 0.65-1.09 | 0.86 (0.71-1.04) | 2.4 |

#### Cerebrovascular hemorrhage

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 7.50 | 3.59-10.41 | 6.93 (3.45-13.90) | 0.3 |
| Day 2 to \<1 month | 1.28 | 0.85-1.78 | 1.19 (0.87-1.63) | 0.6 |
| 1 to 5 months | 1.01 | 0.82-1.26 | 0.86 (0.73-1.01) | 0.7 |
| 6 to 11 months | 1.01 | 0.83-1.13 | 0.96 (0.82-1.13) | 0.6 |
| 12 months or more | 1.01 | 0.51-1.20 | 0.85 (0.63-1.14) | 0.7 |

#### Other cerebrovascular

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 11.38 | 5.41-12.32 | 8.75 (5.95-12.90) | 1.1 |
| Day 2 to \<1 month | 1.27 | 0.82-1.51 | 1.15 (0.94-1.40) | 2.0 |
| 1 to 5 months | 1.08 | 0.98-1.10 | 1.01 (0.92-1.11) | 3.9 |
| 6 to 11 months | 1.08 | 0.92-1.11 | 1.02 (0.93-1.12) | 4.7 |
| 12 months or more | 1.01 | 0.63-0.99 | 0.77 (0.64-0.94) | 4.1 |

#### Cerebral infarction

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 10.35 | 6.04-12.13 | 8.29 (5.97-11.50) | 1.6 |
| Day 2 to \<1 month | 1.27 | 1.04-1.45 | 1.16 (0.99-1.37) | 3.2 |
| 1 to 5 months | 1.01 | 0.89-1.08 | 0.88 (0.81-0.96) | 3.2 |
| 6 to 11 months | 1.01 | 0.90-1.00 | 0.91 (0.83-0.99) | 3.1 |
| 12 months or more | 1.01 | 0.69-0.91 | 0.74 (0.62-0.88) | 2.6 |

#### Aneurysm dissection

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 2.98 | none | 3.85 (1.83-8.09) | 0.2 |
| Day 2 to \<1 month | 1.01 | 0.72-1.15 | 0.90 (0.68-1.21) | 0.2 |
| 1 to 5 months | 1.01 | 0.88-1.18 | 1.02 (0.90-1.15) | 0.3 |
| 6 to 11 months | 1.01 | 0.79-1.15 | 0.98 (0.86-1.11) | 0.4 |
| 12 months or more | 1.01 | 0.54-1.18 | 0.73 (0.56-0.95) | 0.9 |

#### Heart failure

| Window | True HR | Simulated, 95% PI | Published HR (95% CI) | Extra per 100,000 |
|:---|---:|---:|---:|---:|
| Day 0-1 | 16.21 | 7.85-20.89 | 10.80 (6.51-18.00) | 0.8 |
| Day 2 to \<1 month | 1.11 | 0.61-1.36 | 1.00 (0.73-1.37) | 0.9 |
| 1 to 5 months | 1.01 | 0.84-1.10 | 0.87 (0.73-1.01) | 1.1 |
| 6 to 11 months | 1.01 | 0.76-1.19 | 0.73 (0.62-0.87) | 0.9 |
| 12 months or more | 1.01 | 0.54-1.31 | 0.57 (0.41-0.80) | 0.7 |

### Agreement with the published estimates

The true hazard ratios are fitted to the published estimates. So the agreement below shows that the fit succeeded. It is not an independent validation.

- 55 of the 58 published estimates that have a band lie inside it. The 3 outside are cerebral infarction at 1 to 5 months, ischemic heart disease at 1 to 5 months, and heart failure at 6 to 11 months.
- 2 day-1 estimates have no band, because fewer than 6 of the 10 cohorts had 5 or more events: aneurysm dissection and inflammatory heart disease.
- At 12 months or more, the geometric mean of the 12 simulated estimates is 0.821, against 0.825 published.
- At 12 months or more, a mean of 5.7 of the 12 outcomes per cohort have a 95% CI below 1, against 6 published.

### Why the estimates fall below 1

The difference between the true and the estimated hazard ratios is bias from five sources in the simulated study:

- unrecorded differences in cardiovascular risk between persons of the same age (shared frailty);
- infections that testing did not record, in the comparison group;
- censoring at the first diagnosis in any of the 15 cardiovascular groups;
- selection into early infection: persons infected before the Omicron wave have a lower baseline cardiovascular risk;
- death at infection: it is more probable at older ages and higher frailty, so it removes high-risk infected persons from follow-up.

### Absolute excess

`Run.R` simulates each cohort again without infection. The persons, their natural death times and their event thresholds stay the same. With infection, a mean of 2,381 more persons per cohort have a CVD diagnosis. That is 52.8 per 100,000 persons, or 3.36% of the 70,885 persons with a CVD diagnosis in the paper. Of them, 1,488 have a recorded infection and 893 an unrecorded one.

In the 12 months after a recorded infection, the risk of a first CVD diagnosis is 0.536% with infection and 0.486% without. The base is the 665,687 persons per cohort with a recorded infection whose study period runs at least 12 months past it.

The letter figure is [figures/forest_12m.png](figures/forest_12m.png). It shows the window of 12 months or more only, with rows sorted by the published hazard ratio. Its right column repeats the extra diagnoses at 12 months or more from the tables above.

## The model

Each cohort follows its persons from 1 March 2020 to 31 December 2022. The table gives every input, with its source or with the word “assumed”.

| Component | How `Run.R` draws it | Source |
|----|----|----|
| Age | Age at the start of follow-up, from a monotone spline. Its knots are the median 34.2, the quartiles 18.5 and 52.3, and the Table 1 band edges at 40 and 65 years. Above the last band, the knots are 72 years (95th percentile), 97 (99.9th) and 105. | Paper, Results and Table 1. The knots above 65 years are assumed. |
| Follow-up | Start of follow-up, from a monotone spline. Its quartiles are day 1036 minus the paper’s follow-up quartiles. Its 95th percentile is day 759.7, and its maximum is day 1035. Follow-up stops at the first CVD diagnosis, death or 31 December 2022. | Paper, Results, for the quartiles. The 95th percentile is fitted to the person-years of Supplementary Tables 1 and 7. |
| Frailty | One log-normal multiplier, SD 1.65 on the log scale. It multiplies all 15 CVD hazards, the death hazard and the probability of death at infection. | Assumed. |
| Outcomes | 15 CVD groups. Each hazard is Gompertz in age, doubling every 7.5 years, times frailty. The level of each outcome is its events per test-negative person-year. One scale factor sets the test-negative first-event rate to 8.766 per 1000 person-years. | Supplementary Table 1: 54,247 events in 6,188,401 test-negative person-years. The age slope is assumed. |
| Calendar step | From 11 March 2022, every outcome hazard is 1.175 times higher. The multiplier and the scale factor are solved against the published test-negative rates of 8.119 per 1000 person-years before that day and 10.870 after. They are solved without an infection effect or unrecorded infections. | Supplementary Table 7, and Table 1 minus Table 7. The step is assumed. |
| Recorded infection | The model draws one for 2,698,261 of 4,508,489 persons (59.8%). The probability is proportional to the SSI confirmed cases per resident of the person’s age group. The date is piecewise uniform: 16.4% to 21 December 2021, the Omicron wave to 10 March 2022, and 3.9% after, split by month. | Paper, Results. Statens Serum Institut (SSI) cases by age group and first infections by month, in `data/ssi/`. The dates to 10 March 2022 are fitted with the start of follow-up. |
| Unrecorded infection | 74.5% of persons without a recorded infection get an unrecorded one, so one third of all infections are unrecorded. They follow the same age pattern, with each probability capped at 1. Their dates follow the recorded curve, weighted by (1 - p)/p. The timing weight p is 0.40 to 10 March 2022 and 0.05 after. It sets only the dates, and the one-third share sets the number. Dates after 10 March 2022 follow the weekly SSI wastewater index. | Erikstrup et al. 2022 for the one third. SSI wastewater index in `data/ssi/`. The timing weights are assumed. |
| Selection | A person 1 SD of log frailty less frail has exp(0.15) times the odds of a recorded infection before 21 December 2021 (`pre_frailty` 0.15). | Assumed. |
| Effect of infection | Each outcome has five true hazard ratios, one per window: day 0-1, day 2 to 1 month, months 1-5, 6-11 and 12 or more. Recorded and unrecorded infections carry them alike. Arterial embolism, cardiac arrest and cardiomyopathy have the 3 lowest rates, and are not among the 12 fitted outcomes. They share one row: 7.29, 1.27, 1.11, 1.11, 1.02. | Fitted so the simulated estimates of the 12 outcomes reproduce Supplementary Table 1. The shared row of the other 3 is assumed. |
| Death | Gompertz in age: 0.8 per 1000 person-years at the median age, doubling every 8 years, times frailty. Only a living person is infected. An infection dated after natural death moves to a living person of the same age group, or of any age group if none is free. An infection during follow-up kills with probability 0.03 at age 80 and average frailty, more at older ages and higher frailty. Every person is alive at the start of follow-up, so an earlier infection does not kill. | Assumed. |
| Analysis | Person-time split by window since the recorded test, 5-year age band and calendar year, censored at the first CVD diagnosis. One Poisson model per outcome. | Approximates the paper’s Cox model. Sex and comorbidity are not simulated. |

## How the simulated cohort matches the known data

`Run.R` compares the seed-1 cohort with the published data. In the Role column, “input” marks a quantity that the model takes from the source or is fitted to. “Check” marks a quantity that no input sets.

| Quantity | Simulated | Published | Source | Role |
|:---|:---|:---|:---|:---|
| Persons | 4,508,489 | 4,508,489 | Paper, Results | input |
| Age, median (IQR), years | 34.2 (18.5-52.3) | 34.2 (18.5-52.3) | Paper, Results | input |
| Age at first test \<40 / 40-64 / 65+ | 57.7% / 32.2% / 10.1% | 57.7% / 32.2% / 10.1% | Paper, Table 1 | input |
| Person-years | 8,813,676 | 8,909,627 | Paper, Results | input |
| Follow-up, median (IQR), months | 25.1 (21.4-27.4) | 25.2 (21.7-27.5) | Paper, Results | input |
| Recorded positives | 2,698,261 (59.8%) | 2,698,261 (59.8%) | Paper, Results | input |
| Recorded share by SSI age group | within 1.7 binomial SE of the target in 9 of 9 groups | the age pattern of SSI cases | SSI | input |
| Test-negative person-years | 6,157,664 | 6,188,401 | Supplementary Table 1 | input |
| Test-negative CVD rate per 1000 person-years, whole period | 8.947 | 8.766 | Supplementary Table 1 | input |
| Same, to 10 March 2022 | 8.253 | 8.119 | Supplementary Table 7 | input |
| Same, after 10 March 2022 | 11.281 | 10.870 | Supplementary Tables 1 and 7 | input |
| Baseline rate per outcome, simulated / published | 0.94 to 1.13 | 1 | Supplementary Table 1 | input |
| Person-years, 12 months or more, to 10 March 2022 | 43,515 | 50,715 | Supplementary Table 7 | input |
| Persons with a first CVD diagnosis | 71,025 | 70,885 | Paper, Results | check |
| … after a recorded positive | 15,934 | 16,475 | Paper, Results | check |
| Recorded infections in 2020 / 2021 / 2022 | 6.2% / 20.6% / 73.2% | 5.3% / 20.2% / 74.5% | SSI first infections | check |
| Infected, 1 November 2021 to 15 March 2022 | 66.2% of all persons | 66% (63-70%) of all blood donors aged 17-72 | Erikstrup et al. 2022 | check |
| Unrecorded share of those infections | 26.3% | one third | Erikstrup et al. 2022 | check |

![Recorded and unrecorded infections per month, against SSI](figures/cohort_infection_timing.png)

The SSI line is the national count of first infections per month, scaled so that its total over the study period equals the paper’s 2,698,261 recorded positives. So the figure compares the timing, and not the level. The grey vertical line marks 10 March 2022, when widespread testing ended.

![Person-years per exposure window, simulated and published](figures/cohort_person_time.png)

![Baseline CVD rate per outcome in test-negative time, simulated and published](figures/cohort_baseline_rates.png)

The cohort agrees closely with the published data in age, follow-up, person-time, recorded infection, baseline rates and CVD diagnoses. Three differences are larger:

- **Recorded infections after widespread testing ended.** Widespread testing ended on 10 March 2022, but SSI continued to record first infections. The model puts 3.9% of recorded infections after 10 March 2022. That share is fitted to the person-years of the paper’s supplementary tables, not to SSI. The national SSI series has more infections in this period. From April to December 2022, the model has 0.33 times the scaled SSI count. In March 2022 it has 1.38 times.
- **Recorded infections before the Omicron wave.** The model spreads them evenly up to 21 December 2021, so the 2020 and 2021 waves are absent from the figure above. The share per calendar year still agrees with SSI.
- **Person-years at 12 months or more to 10 March 2022.** The cohort has 0.86 times the published person-years in Supplementary Table 7. The infection dates and the start of follow-up are fitted to these person-years, but the fit does not reach them.

The Erikstrup comparison is weak. Erikstrup counts healthy blood donors aged 17-72, and the cohort has all ages. The model has first infections only, and a donor can seroconvert on a reinfection. The one-third unrecorded share applies to the whole study period, but the timing weights move unrecorded infections later. So in the Erikstrup window only 26.3% of infections are unrecorded.

## Assumptions and limits

The result depends on these assumptions:

- **The true hazard ratios are fitted to the published estimates.** The simulation shows that bias can produce the published estimates from harmful true effects. It does not show that the fitted hazard ratios are the true effects of infection.
- **The selection into early infection is assumed.** At `pre_frailty` 0.15, persons with a recorded infection before the Omicron wave are less frail than persons infected later. Over the 10 cohorts, the ratio of mean frailty is 0.80. Frailty is the unrecorded cardiovascular risk of persons of the same age. The ratio of modelled baseline CVD rates, which includes age, is 0.91. No source measures this in the cohort.
- **One third of infections are unrecorded.** Erikstrup et al. 2022 measured this in blood donors aged 17-72 during one Omicron wave. The model applies it to all ages and the whole study period.
- **The frailty SD of 1.65 is assumed.** At the same age, the 95th percentile of the hazard multiplier is 227.7 times the 5th. The same frailty drives every outcome, death and death at infection, so these risks rise together.
- **The timing weights 0.40 and 0.05, the death hazard and the death at infection are assumed.**
- **The calendar step is assumed.** It models a rise in the test-negative rate that the simulation does not explain.
- **The analysis adjusts for calendar year only.** This approximates the paper’s adjustment for calendar period. The calendar step falls inside 2022, so this adjustment cannot remove it.
- **The national SSI series are applied to the cohort.** The SSI cases by age group and the wastewater index both include reinfections.
- **The analysis approximates the paper’s Cox model.** It uses Poisson models in 5-year age bands and calendar year, without sex or comorbidity.

## How to run it

1.  Run `Rscript Run.R` from the repository root. It needs R 4.6 with data.table, ggplot2, patchwork and knitr.
2.  Run `quarto render README.qmd` to rebuild this README from `results/run.rds`.

`Run.R` sets the model in its CFG section, and sources the functions in `R/functions.R`. It checks the md5 of the 3 SSI files in `data/ssi/`. It stops if a value in CFG differs from the SSI file that it comes from. It prints its results, draws the 5 figures into `figures/`, and saves the results to `results/run.rds`.

In one measured run with 2 workers, it took 657 seconds on a 20-core Linux machine. The largest R process peaked at 11.05 GiB of memory, and all R processes together at 21.28 GiB. Memory was sampled every second, so these are lower bounds. Run it on a machine with at least 25 GB of free memory. On Windows the run uses 1 worker.

## Layout

| Path | Content |
|----|----|
| `Run.R` | The model settings, the published values with their sources, the analysis and the output. |
| `R/functions.R` | The simulation, the analysis and the cohort description. |
| `README.qmd` | The source of this README. |
| `data/ssi/` | The SSI source files, byte for byte, with their sources and md5s. |
| `figures/` | The 5 figures that `Run.R` draws. |
| `results/run.rds` | The results that `Run.R` saves and `README.qmd` reads. |

## Licence

The code and figures are under the MIT licence. See `LICENSE`. The MIT licence does not cover the SSI files in `data/ssi/`: they are © Copyright Statens Serum Institut.
