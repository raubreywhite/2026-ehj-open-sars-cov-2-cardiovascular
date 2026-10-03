

# Long-term cardiovascular risk after SARS-CoV-2 infection: bias in Boyd et al. 2026

## Summary

Boyd et al. 2026, “SARS-CoV-2 infection and long-term risk of cardiovascular and renal morbidity” (EHJ Open, doi:10.1093/ehjopen/oeag121), followed 4.5 million Danish persons with a SARS-CoV-2 test. At 12 months or more after a positive test, 11 of the 12 most common cardiovascular outcomes have a hazard ratio below 1, and 6 have a 95% CI below 1. The authors conclude that infection carries little long-term cardiovascular risk.

**The published estimates do not support that conclusion.** We simulated a cohort in which infection increases the risk of every cardiovascular outcome at every time since infection. At 12 months or more, the true hazard ratios average 1.19. Analysed as the paper does, the simulated cohort gives 0.83, against 0.82 published. 57 of the 59 published estimates lie inside the simulated 95% prediction intervals, and the estimates by virus variant also match. A design that reports harm as protection cannot be used to rule out harm.

Three biases explain the gap. The largest is in the comparison group. By the end of 2022, 90% of the cohort had been infected. Those still uninfected were at higher underlying cardiovascular risk, because higher-risk persons were slightly more likely to avoid infection.

## Terms

| Term | Meaning |
|----|----|
| Time window | Time since the first positive test: day 0-1, day 2 to \<1 month, 1 to 5 months, 6 to 11 months, or 12 months or more. |
| Recorded infection | An infection with a positive test. |
| Unrecorded infection | An infection without a positive test. The person stays in the comparison group. |
| Comparison group | Persons whose tests have all been negative so far. The paper calls it test-negative person-time. |
| Underlying cardiovascular risk | The part of a person’s cardiovascular risk that age does not explain, for example from smoking, body weight or medication. The paper cannot adjust for it. |
| True hazard ratio | In the simulation, the factor by which infection multiplies a person’s hazard of a cardiovascular diagnosis. |
| Average | The geometric mean over the 12 most common outcomes. At day 0-1 it is over the 8 outcomes with an estimate in every simulated cohort. |
| Prediction interval | The range in which 95% of estimates from one study of this size would fall. |
| Omicron wave | 21 December 2021 to 10 March 2022. Widespread PCR testing in Denmark ended on 10 March 2022. |
| Supplementary Tables 1, 5 and 7 | Boyd et al.’s estimates by time window (Table 1), by virus variant (Table 5), and with follow-up ended on 10 March 2022 (Table 7). |

## What the published estimates show

Follow-up ended on 31 December 2022, so the window of 12 months or more holds only persons infected in 2020 and 2021.

- **By time window (Supplementary Table 1):** the published average is 7.77 at day 0-1, 1.20 at day 2 to \<1 month, 1.00 at 1 to 5 months, 0.98 at 6 to 11 months, and 0.82 at 12 months or more.
- **By variant (Supplementary Table 5):** in the first year after infection, the published average is 1.00 for the original strain, 1.02 for alpha, 1.05 for delta, and 0.97 for omicron.

These patterns point to bias, not protection. Persons infected in 2020 and 2021 show no lower risk in their first year. Their estimates fall below 1 only later, in 2022. A protective effect that starts after a year is implausible. If persons infected early were healthier from the start, their first-year estimates would be lower too. A comparison group that becomes riskier in 2022 explains both patterns.

## How the simulation works

1.  **Cohort.** `Run.R` simulates 10 cohorts of 4,508,489 persons, with the paper’s age distribution and follow-up from 1 March 2020 to 31 December 2022.
2.  **Infections.** 2,698,261 persons get a recorded infection, as in the paper, timed by the Danish case counts. One third of all infections are unrecorded.
3.  **True effects.** Each infection multiplies the hazard of each of 15 cardiovascular outcomes by a true hazard ratio for each time window. Every true hazard ratio is at least 1.01.
4.  **Biases.** Unrecorded infections, depletion of susceptibles and avoidance of infection are built in. The model also lets infection kill some persons. This changes the estimates by 3% or less, so it is not treated as a bias.
5.  **Analysis.** Each cohort is analysed approximately as the paper does: one Poisson model per outcome, with person-time split by time window, 5-year age band and calendar year.

The true hazard ratios are fitted, so that the simulated estimates match the published ones. They are at least 1.01 and do not increase with time since infection. [The model](#the-model) gives every input and its source.

## Results

### By time window (Supplementary Table 1)

At 12 months or more, every true hazard ratio is at least 1.01, and the average is 1.19. Yet 11 of the 12 published estimates are below 1. The simulation reproduces them: 57 of the 59 published estimates (97%) lie inside the 95% prediction interval.

![True, simulated and published hazard ratios](figures/forest_outcome_truths.png)

Red dot: true hazard ratio. Green band: 95% prediction interval for one study. Black square and bar: published estimate and 95% CI.

| Time window        | Simulated average | Published average |
|--------------------|------------------:|------------------:|
| Day 0-1            |              7.91 |              7.77 |
| Day 2 to \<1 month |              1.23 |              1.20 |
| 1 to 5 months      |              1.01 |              1.00 |
| 6 to 11 months     |              0.94 |              0.98 |
| 12 months or more  |              0.83 |              0.82 |

- **Outside the prediction interval:** aneurysm dissection at day 2 to \<1 month and ischemic heart disease at 6 to 11 months.
- **No interval:** inflammatory heart disease at day 0-1, because too few cohorts had 5 or more events. At day 0-1, the averages are over the 8 outcomes with an estimate in every simulated cohort.
- **95% CI below 1 at 12 months or more:** 4.4 of the 12 outcomes per simulated study, against 6 published.

For several outcomes the published estimate at 6 to 11 months is higher than at 1 to 5 months. The true hazard ratios cannot increase with time, so the fit gives both windows the same value. The appendix gives the [values per outcome](#appendix-values-per-outcome).

### By variant (Supplementary Table 5)

`Bias.R` analyses the simulated cohorts as Supplementary Table 5 does: the first year after infection, by the variant period of the positive test.

| Variant  | Simulated average | Published average | Outcomes |
|:---------|------------------:|------------------:|---------:|
| Original |              1.04 |              1.00 |       12 |
| Alpha    |              1.01 |              1.02 |       12 |
| Delta    |              0.98 |              1.03 |       11 |
| Omicron  |              0.95 |              0.97 |       12 |

Each average is over the outcomes with an estimate in every simulated cohort, so delta has 11. The simulation reproduces the pattern: about 1 before the Omicron wave, and lower for omicron. Omicron against the earlier variants gives 0.94 in the simulation and 0.96 in the paper. The largest difference is for delta, 0.98 against 1.03.

### Extra diagnoses

`Run.R` also simulates each cohort without infection, with everything else the same. Infection adds 6,181 persons with a cardiovascular diagnosis per cohort: 137.1 per 100,000 persons, or 8.7% of the 70,885 persons diagnosed in the paper. Of them, 3,738 have a recorded infection and 2,443 an unrecorded one.

In the 12 months after a recorded infection, the risk of a first cardiovascular diagnosis is 0.586% with infection and 0.439% without. The base is the 665,704 persons per cohort followed for at least 12 months after a recorded infection.

[figures/forest_12m.png](figures/forest_12m.png) is the figure for the letter. It shows 12 months or more, and the extra diagnoses per 100,000 persons over 24 months after a recorded infection, among the 101,299 persons per cohort followed that long.

## Why the estimates fall below 1

Each bias below has an example with 100 persons per group. In each example, infection multiplies every person’s risk by 1.20, so the true ratio is 1.20. A high-risk person has a 50% risk and a low-risk person a 5% risk. After each example comes what the bias does in the simulation, at 12 months or more. The biases are added one at a time, in this order.

### Unrecorded infections

About one third of infections in the Omicron wave were not recorded, so some infected persons stay in the comparison group and raise its risk. In the example, 50 of 150 infections are unrecorded, and the study sees 1.09.

**In the simulation:** 1.18 to 1.04.

![Unrecorded infections](figures/illustration_unrecorded.png)

### Depletion of susceptibles

Infection brings diagnoses forward, most of all in high-risk persons. A person diagnosed is no longer followed. So after a year, the infected group has fewer high-risk persons left than the comparison group. In the example, year 1 shows the true ratio, and year 2 shows 1.10. The paper adds to this: it stops follow-up for every outcome at the first diagnosis of any of the 15 outcomes.

**In the simulation:** 1.04 to 0.98. Depletion needs time to act, so it lowers the estimates only from 6 months on. The table gives the simulated average without and with depletion, with unrecorded infections in both:

| Time window        | Without depletion | With depletion |
|:-------------------|------------------:|---------------:|
| Day 0-1            |              8.87 |           9.15 |
| Day 2 to \<1 month |              1.43 |           1.46 |
| 1 to 5 months      |              1.18 |           1.18 |
| 6 to 11 months     |              1.16 |           1.10 |
| 12 months or more  |              1.04 |           0.98 |

![Depletion of susceptibles](figures/illustration_depletion.png)

### Avoidance of infection

If persons at high underlying risk are more careful, fewer of them are infected. By 2022, when most of the population had been infected, the persons still uninfected include a larger share of high-risk persons. The comparison group then looks riskier than the infected group. In the example, the study sees 0.81.

**In the simulation:** 0.98 to 0.83. This is the largest bias.

The difference per person is small: a person with twice the underlying risk of another of the same age has a 3.3% lower chance of infection. It matters because almost everyone was infected, so the comparison group became the persons who avoided infection. In the simulation, the average underlying risk of the comparison group was 15% higher than that of infected persons in 2020, and 28% higher in late 2022.

In the paper, 83% of the person-time at 12 months or more lies after 10 March 2022, so these estimates are made against the comparison group of 2022.

![Avoidance of infection](figures/illustration_avoidance.png)

### All three together

`Bias.R` switches the biases on one at a time in the simulated cohorts, with the true hazard ratios unchanged.

![From the true hazard ratios to the published estimates](figures/bias_steps.png)

Each point is the simulated average, with its 95% prediction interval for one study. At day 0-1, 8 outcomes have an estimate in every cohort. The order of the steps changes the size of each step, but not the first or last point.

| Step | Day 0-1 | Day 2 to \<1 month | 1 to 5 months | 6 to 11 months | 12 months or more |
|:---|---:|---:|---:|---:|---:|
| True | 11.61 | 1.62 | 1.34 | 1.32 | 1.19 |
| All biases off | 10.03 | 1.60 | 1.33 | 1.31 | 1.18 |
| \+ unrecorded infections | 8.87 | 1.43 | 1.18 | 1.16 | 1.04 |
| \+ depletion of susceptibles | 9.15 | 1.46 | 1.18 | 1.10 | 0.98 |
| \+ avoidance of infection | 7.91 | 1.23 | 1.01 | 0.94 | 0.83 |
| Published | 7.77 | 1.20 | 1.00 | 0.98 | 0.82 |

In the first month, the estimate with all biases off is already below the true value. The paper starts follow-up 30 days after the first test, so for a person whose first test was positive, these windows lie about a month after infection.

## Assumptions

The argument needs one thing: that a cohort in which infection is harmful can give the published estimates. The true hazard ratios and the size of each bias do not need to be exact.

- **The true hazard ratios are fitted** to the published estimates. Their average is 11.61 at day 0-1, 1.62 at day 2 to \<1 month, 1.34 at 1 to 5 months, 1.32 at 6 to 11 months, and 1.19 at 12 months or more.
- **Avoidance of infection is assumed.** No source measures it in Denmark. Its strength was chosen to reproduce the omicron against pre-Omicron contrast of Supplementary Table 5, in a screen of model variants on seeds 1 to 5 of the 10. Any other mechanism that makes the comparison group riskier in 2022 acts in the same way.
- **How much underlying risk varies is assumed.** At the same age, the 95th percentile is 227.7 times the 5th. It stands for the risk left after the paper’s adjustment for age, sex and comorbidity.
- **One third of infections are unrecorded.** Erikstrup et al. 2022 measured this in blood donors aged 17-72 during the Omicron wave. The model applies it to all ages and the whole period.
- **The timing of unrecorded infections and the death hazards are assumed.**
- **The analysis approximates the paper’s Cox model.** The paper uses age as the time scale and adjusts for sex and comorbidity. The simulation uses 5-year age bands, within which the hazard still rises 1.59-fold, and has no sex or comorbidity.
- **The national case counts are applied to the cohort.** They include reinfections.

Not in the simulation:

- **Washout:** the paper starts follow-up 30 days after the first test, so for a positive first test the acute phase is left out.
- **Testing at hospital contact:** a person admitted with a cardiovascular event is tested on arrival. The paper names this as the likely cause of its day 0-1 results.
- **Vaccination and variant severity:** infections in 2020 and early 2021 came before most vaccination, and were on average more severe.
- **Reinfections:** the simulation has first infections only.
- **Delayed or caught-up diagnoses in 2020 and 2021.**

## The model

Every input, with its source. The name in brackets is the setting in `Run.R`.

| Component | How `Run.R` draws it | Source |
|----|----|----|
| Age | Age at the start of follow-up, from a monotone spline. Its knots are the median 34.2, the quartiles 18.5 and 52.3, and the Table 1 band edges at 40 and 65 years. Above the last band, the knots are 72 years (95th percentile), 97 (99.9th) and 105. | Paper, Results and Table 1. The knots above 65 years are assumed. |
| Follow-up | Start of follow-up, from a monotone spline. Day 1036 is 31 December 2022, the end of follow-up, so the quartiles of the start day are day 1036 minus the paper’s follow-up quartiles. Its 95th percentile is day 759.7, and its maximum is day 1035. Follow-up stops at the first cardiovascular diagnosis, death or 31 December 2022. | Paper, Results, for the quartiles. The 95th percentile is fitted to the person-years of Supplementary Tables 1 and 7. |
| Underlying cardiovascular risk (`fsd`) | One log-normal multiplier per person, SD 1.65 on the log scale, mean 1. It multiplies all 15 cardiovascular hazards, the death hazard and the probability of death caused by infection. | Assumed. |
| Outcomes | 15 cardiovascular groups. Each hazard is Gompertz in age, doubling every 7.5 years, times the underlying risk, times one common scale (`sc`). The level of each outcome is its rate in the comparison group in Supplementary Table 1. The scale is fitted with the true hazard ratios, so that the first-event rate in the comparison group is 8.766 per 1000 person-years. | Supplementary Table 1: 54,247 events in 6,188,401 person-years in the comparison group. The age slope is assumed. |
| Recorded infection (`p_infected`) | 2,698,261 of 4,508,489 persons (59.8%). Persons are drawn with a weight: the SSI confirmed cases per resident of the person’s age group, times the avoidance factor. The probability of being drawn is proportional to the weight, and at most 1. The date is piecewise uniform: 16.4% to 21 December 2021, the Omicron wave to 10 March 2022, and 3.9% after, split by month. | Paper, Results. SSI cases by age group and first infections by month, in `data/ssi/`. The dates to 10 March 2022 are fitted with the start of follow-up. |
| Unrecorded infection (`contam`) | 74.5% of persons without a recorded infection get an unrecorded one, so one third of all infections are unrecorded. They follow the same age pattern and avoidance factor. Their dates follow the recorded curve, weighted by (1 - p)/p. The timing weight p is 0.40 to 10 March 2022 and 0.05 after. Dates after 10 March 2022 follow the weekly SSI wastewater index. | Erikstrup et al. 2022 for the one third. SSI wastewater index in `data/ssi/`. The timing weights are assumed. |
| Avoidance of infection (`avoid`) | The weight of a person in the draws of recorded and unrecorded infections is multiplied by exp(-0.08 z), where z is the person’s log underlying risk in SD units. Twice the underlying risk gives a 3.3% lower weight. An infection that moves to another person does not use the factor. A recorded infection moves when its first person died before it, and an unrecorded one also when its first person already has a recorded infection: 0.25% of recorded and 0.80% of unrecorded infections. | Assumed. Chosen in a screen of model variants to reproduce the Omicron against pre-Omicron contrast of Supplementary Table 5. |
| Effect of infection (`tm`) | Five true hazard ratios per outcome, one per time window. Recorded and unrecorded infections carry them alike. Arterial embolism, cardiac arrest and cardiomyopathy have the 3 lowest rates and are not among the 12 fitted outcomes. Their true hazard ratios are the geometric mean of the 12 fitted ones, per window. | Fitted so that the simulated estimates of the 12 outcomes reproduce Supplementary Table 1. |
| Death | Gompertz in age: 0.8 per 1000 person-years at the median age, doubling every 8 years, times the underlying risk. Only a living person is infected. An infection during follow-up kills with probability 0.03 at age 80 and average underlying risk (`p80`), more at older ages and higher underlying risk. Every person is alive at the start of follow-up, so an earlier infection does not kill. | Assumed. |
| Analysis | Person-time split by time window since the recorded test, 5-year age band and calendar year, and stopped at the first diagnosis of any of the 15 outcomes. One Poisson model per outcome. | Approximates the paper’s Cox model. |

A person’s cardiovascular hazard changes only with age, underlying risk and infection. There is no calendar-time effect.

## How the simulated cohort matches the paper

The first simulated cohort against the paper and national data. “Input” marks a quantity the model is set or fitted to. “Check” marks one that no input sets.

| Quantity | Simulated | Published | Source | Role |
|:---|:---|:---|:---|:---|
| Persons | 4,508,489 | 4,508,489 | Paper, Results | input |
| Age, median (IQR), years | 34.2 (18.5-52.3) | 34.2 (18.5-52.3) | Paper, Results | input |
| Age at first test \<40 / 40-64 / 65+ | 57.7% / 32.2% / 10.1% | 57.7% / 32.2% / 10.1% | Paper, Table 1 | input |
| Person-years | 8,815,402 | 8,909,627 | Paper, Results | input |
| Follow-up, median (IQR), months | 25.1 (21.4-27.4) | 25.2 (21.7-27.5) | Paper, Results | input |
| Recorded positives | 2,698,261 (59.8%) | 2,698,261 (59.8%) | Paper, Results | input |
| Recorded share by SSI age group | within 1.7 binomial SE of the target in 9 of 9 groups | the age pattern of SSI cases | SSI | input |
| Person-years in the comparison group | 6,156,784 | 6,188,401 | Supplementary Table 1 | input |
| Comparison group, first cardiovascular diagnoses per 1000 person-years, whole period | 8.767 | 8.766 | Supplementary Table 1 | input |
| Same, to 10 March 2022 | 8.073 | 8.119 | Supplementary Table 7 | check |
| Same, after 10 March 2022 | 11.102 | 10.870 | Supplementary Tables 1 and 7 | check |
| Rate per outcome in the comparison group, simulated / published | 0.94 to 1.13 | 1 | Supplementary Table 1 | input |
| Person-years, 12 months or more, to 10 March 2022 | 43,536 | 50,715 | Supplementary Table 7 | input |
| Persons with a first cardiovascular diagnosis | 69,053 | 70,885 | Paper, Results | check |
| … after a recorded positive | 15,079 | 16,475 | Paper, Results | check |
| Recorded infections in 2020 / 2021 / 2022 | 6.2% / 20.6% / 73.2% | 5.3% / 20.2% / 74.5% | SSI first infections | check |
| Infected, 1 November 2021 to 15 March 2022 | 66.2% of all persons | 66% (63-70%) of all blood donors aged 17-72 | Erikstrup et al. 2022 | check |
| Unrecorded share of those infections | 26.3% | one third | Erikstrup et al. 2022 | check |

![Recorded and unrecorded infections per month, against SSI](figures/cohort_infection_timing.png)

The SSI line is the national count of first infections, scaled to the paper’s 2,698,261 recorded positives, so it compares timing only. The grey line marks 10 March 2022, when widespread testing ended.

![Person-years per time window, simulated and published](figures/cohort_person_time.png)

![Rate per outcome in the comparison group, simulated and published](figures/cohort_baseline_rates.png)

The figures call the comparison group “test-negative”, as the paper does. The cohort agrees closely with the paper. The larger differences:

- **Recorded infections after 10 March 2022.** The model puts 3.9% of recorded infections after that day, fitted to the paper’s person-years rather than to SSI. From April to December 2022 it has 0.33 times the scaled SSI count, and in March 2022 1.38 times.

- **Recorded infections before the Omicron wave** are spread evenly, so the 2020 and 2021 waves are absent. The share per calendar year still agrees with SSI.

- **Person-years at 12 months or more to 10 March 2022:** 0.86 times the published value in Supplementary Table 7.

- **Persons with a first cardiovascular diagnosis after a recorded positive:** 15,079, against 16,475.

- **Unrecorded share in the Erikstrup window:** 26.3%, because the model moves unrecorded infections later. Erikstrup et al. studied blood donors aged 17-72, and a donor can seroconvert on a reinfection.

## Appendix: values per outcome

The values in the first figure, in its order. The last column is the extra rate of first diagnoses caused by infection, per 100,000 person-years, among persons with a recorded infection, each simulated with and without infection. It is the mean over the 10 cohorts. Where the true hazard ratio is 1.01, it can fall just below 0, from Monte Carlo error and because earlier diagnoses leave fewer persons at risk later.

### Pulmonary embolism

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 27.79 | 11.20-20.53 | 16.10 (11.00-23.50) | 568.7 |
| Day 2 to \<1 month | 7.46 | 3.63-4.89 | 4.14 (3.57-4.79) | 138.8 |
| 1 to 5 months | 2.03 | 1.08-1.29 | 1.21 (1.08-1.36) | 23.3 |
| 6 to 11 months | 1.63 | 0.76-1.04 | 0.89 (0.78-1.02) | 14.1 |
| 12 months or more | 1.47 | 0.60-1.15 | 0.77 (0.61-0.98) | 11.2 |

### Conduction disorders

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 12.01 | 4.54-18.42 | 8.85 (5.01-15.60) | 193.6 |
| Day 2 to \<1 month | 1.39 | 0.68-1.60 | 0.76 (0.53-1.09) | 5.6 |
| 1 to 5 months | 1.39 | 0.92-1.17 | 1.08 (0.94-1.24) | 6.5 |
| 6 to 11 months | 1.39 | 0.83-1.18 | 1.13 (0.98-1.30) | 6.8 |
| 12 months or more | 1.39 | 0.77-1.30 | 1.13 (0.88-1.44) | 5.3 |

### Venous embolism

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 4.64 | 1.28-6.54 | 2.67 (1.48-4.84) | 151.7 |
| Day 2 to \<1 month | 2.00 | 1.14-1.75 | 1.42 (1.22-1.66) | 43.7 |
| 1 to 5 months | 1.60 | 1.03-1.25 | 1.12 (1.03-1.21) | 27.5 |
| 6 to 11 months | 1.60 | 0.99-1.15 | 1.09 (1.00-1.18) | 27.7 |
| 12 months or more | 1.34 | 0.73-1.07 | 0.90 (0.77-1.05) | 12.8 |

### Arrhythmias

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 11.37 | 6.08-9.90 | 7.67 (6.17-9.53) | 1244.2 |
| Day 2 to \<1 month | 1.55 | 1.07-1.30 | 1.19 (1.07-1.32) | 69.3 |
| 1 to 5 months | 1.50 | 1.03-1.15 | 1.07 (1.02-1.12) | 62.3 |
| 6 to 11 months | 1.50 | 0.98-1.09 | 1.05 (1.00-1.11) | 58.6 |
| 12 months or more | 1.31 | 0.81-1.02 | 0.90 (0.82-0.99) | 36.5 |

### Valve disorders

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 3.09 | 1.26-4.74 | 2.51 (1.20-5.28) | 86.7 |
| Day 2 to \<1 month | 1.30 | 0.78-1.32 | 1.03 (0.83-1.28) | 11.3 |
| 1 to 5 months | 1.30 | 0.92-1.13 | 0.98 (0.89-1.08) | 11.7 |
| 6 to 11 months | 1.30 | 0.80-1.13 | 1.03 (0.93-1.14) | 11.2 |
| 12 months or more | 1.22 | 0.70-1.13 | 0.86 (0.71-1.04) | 7.8 |

### Ischemic heart disease

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 8.00 | 4.47-9.07 | 5.84 (4.21-8.11) | 521.6 |
| Day 2 to \<1 month | 1.26 | 0.88-1.21 | 1.01 (0.87-1.17) | 20.6 |
| 1 to 5 months | 1.26 | 0.87-1.09 | 0.93 (0.87-1.00) | 19.2 |
| 6 to 11 months | 1.26 | 0.84-1.01 | 1.03 (0.93-1.14) | 19.1 |
| 12 months or more | 1.22 | 0.75-1.02 | 0.89 (0.78-1.00) | 13.7 |

### Inflammatory heart disease

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 5.52 | none | 4.90 (2.04-11.80) | 40.6 |
| Day 2 to \<1 month | 1.25 | 0.55-1.67 | 1.08 (0.76-1.54) | 2.0 |
| 1 to 5 months | 1.25 | 0.78-1.16 | 0.99 (0.84-1.18) | 2.6 |
| 6 to 11 months | 1.25 | 0.78-1.18 | 1.05 (0.89-1.23) | 2.3 |
| 12 months or more | 1.20 | 0.62-1.17 | 0.91 (0.67-1.23) | 2.1 |

### Other cerebrovascular

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 12.10 | 5.44-12.65 | 8.75 (5.95-12.90) | 402.2 |
| Day 2 to \<1 month | 1.56 | 0.92-1.52 | 1.15 (0.94-1.40) | 21.8 |
| 1 to 5 months | 1.41 | 0.96-1.21 | 1.01 (0.92-1.11) | 16.9 |
| 6 to 11 months | 1.41 | 0.86-1.14 | 1.02 (0.93-1.12) | 15.4 |
| 12 months or more | 1.12 | 0.58-1.10 | 0.77 (0.64-0.94) | 5.0 |

### Cerebrovascular hemorrhage

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 7.82 | 2.63-13.30 | 6.93 (3.45-13.90) | 96.1 |
| Day 2 to \<1 month | 1.47 | 0.70-1.72 | 1.19 (0.87-1.63) | 5.6 |
| 1 to 5 months | 1.17 | 0.79-1.05 | 0.86 (0.73-1.01) | 2.3 |
| 6 to 11 months | 1.17 | 0.69-1.13 | 0.96 (0.82-1.13) | 2.9 |
| 12 months or more | 1.07 | 0.45-1.30 | 0.85 (0.63-1.14) | 0.8 |

### Cerebral infarction

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 12.83 | 5.37-15.67 | 8.29 (5.97-11.50) | 609.5 |
| Day 2 to \<1 month | 1.38 | 0.90-1.44 | 1.16 (0.99-1.37) | 21.7 |
| 1 to 5 months | 1.15 | 0.81-1.02 | 0.88 (0.81-0.96) | 8.0 |
| 6 to 11 months | 1.15 | 0.77-0.92 | 0.91 (0.83-0.99) | 7.7 |
| 12 months or more | 1.01 | 0.59-0.88 | 0.74 (0.62-0.88) | 1.2 |

### Aneurysm dissection

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 3.97 | 2.00-5.69 | 3.85 (1.83-8.09) | 65.0 |
| Day 2 to \<1 month | 1.28 | 0.91-1.38 | 0.90 (0.68-1.21) | 6.1 |
| 1 to 5 months | 1.28 | 0.87-1.18 | 1.02 (0.90-1.15) | 6.7 |
| 6 to 11 months | 1.28 | 0.80-1.08 | 0.98 (0.86-1.11) | 5.5 |
| 12 months or more | 1.01 | 0.46-1.12 | 0.73 (0.56-0.95) | 0.0 |

### Heart failure

| Window | True<br>HR | Simulated<br>95% PI | Published HR<br>(95% CI) | Extra per 100,000<br>person-years |
|:---|---:|---:|---:|---:|
| Day 0-1 | 15.14 | 5.91-17.65 | 10.80 (6.51-18.00) | 195.0 |
| Day 2 to \<1 month | 1.21 | 0.77-1.32 | 1.00 (0.73-1.37) | 4.5 |
| 1 to 5 months | 1.03 | 0.73-1.06 | 0.87 (0.73-1.01) | 0.5 |
| 6 to 11 months | 1.01 | 0.64-1.03 | 0.73 (0.62-0.87) | 0.1 |
| 12 months or more | 1.01 | 0.45-1.19 | 0.57 (0.41-0.80) | 0.7 |

## How to run it

From the repository root, with R 4.6, data.table, ggplot2, patchwork and knitr:

1.  `Rscript Run.R`: the main simulation. Each cohort needs about 11 GiB of memory, and it runs 2 at a time, so it needs at least 25 GB free.
2.  `Rscript Bias.R`: the biases and the analysis by variant. It simulates 130 more cohorts.
3.  `Rscript Illustrations.R`: the three examples.
4.  `quarto render README.qmd`: this README, from `results/`.

`Run.R` checks the md5 of the SSI files in `data/ssi/`, and stops if a value in its CFG section differs from the file it comes from.

## Layout

| Path | Content |
|----|----|
| `Run.R` | The model settings, the published values with their sources, the analysis and the output. |
| `R/functions.R` | The simulation, the analysis and the cohort description. |
| `Bias.R` | What each bias contributes, and the analysis by variant. |
| `Illustrations.R` | The three examples with 100 persons per group. |
| `README.qmd` | The source of this README. |
| `data/ssi/` | The SSI source files, byte for byte, with their sources and md5s. |
| `figures/` | The figures that `Run.R`, `Bias.R` and `Illustrations.R` draw. |
| `results/` | The results that the three scripts save and `README.qmd` reads. |

## Licence

The code and figures are under the MIT licence. See `LICENSE`. The MIT licence does not cover the SSI files in `data/ssi/`: they are © Copyright Statens Serum Institut.
