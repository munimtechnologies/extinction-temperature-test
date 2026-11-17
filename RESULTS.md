# RESULTS: does temperature change predict marine extinction on PhanDA?

Graded 2026-11-17 from the output of job `extinction-temperature-v1` (ran 20:45-20:51 UTC).
Plan: [PREREG.md](PREREG.md) (committed 10:23 UTC). Reading:
[grade_temperature.py](grade_temperature.py) and [PRE-RESULT-NOTE-1.md](PRE-RESULT-NOTE-1.md)
(committed 17:16 UTC, before the job started). The grader's rules were applied to
`results.json` and `boundaries.csv` by hand, exactly as coded; `test_graded_results.py` runs
the grader itself on the published output and checks that it gives the same labels.
Made with AI agents working in MV2's name.

## Verdict

**REPLICATED_ONLY_WITH_MASS_EXTINCTIONS.**

On PhanDA's 485-Myr global temperature history (Judd et al. 2024) and 2026 Paleobiology
Database marine genera, the size of the temperature change across a stage boundary does
predict the extinction rate at that boundary: Spearman rho = 0.21 over 78 boundaries,
p_perm = 0.033, and it survives the registered autocorrelation test (circular shifts,
p = 0.041) and PhanDA's own uncertainty (rho > 0 in 95.8% of 2,000 draws).

**It does not survive removing the six Big Five boundaries:** rho = 0.14, p = 0.11 (n = 72).
Song et al. (2021) reported that the correlation held without the Big Five (rho = 0.47-0.54,
n = 40). On PhanDA at stage resolution it does not reach significance. It also loses
significance once stage durations are accounted for (rho = 0.13, p = 0.14), a registered
check, since longer stages can hold both larger temperature changes and more extinction.

**The 5.2 C threshold does not single out mass extinctions on PhanDA.** Nine of 78
boundaries have a stage-to-stage change of 5.2 C or more; only two of them are Big Five
boundaries (Katian/Hirnantian 7.0 C and Changhsingian/Induan 10.5 C). The other seven
(largest: Tournaisian/Visean, 8.9 C cooling; Gorstian/Ludfordian, 8.5 C warming;
Albian/Cenomanian, 8.1 C warming) saw no mass extinction. This part depends only on public
temperatures and was disclosed before the result (PRE-RESULT-NOTE-1.md). The Big Five do
rank higher than other boundaries by temperature change on average (AUC 0.73, p = 0.031).

## Graded predictions

| # | Prediction | Prior | Holds? |
|---|---|---|---|
| S | Sanity: >= 4 of the Big Five among the 10 highest-mu3t boundaries, and >= 60 boundaries | 85% | YES |
| H1 | abs(dT) predicts mu3t across all boundaries | 40% | YES |
| H1s | H1 survives autocorrelation (circular shifts) | 30% | YES |
| H2 | H1 holds with the 6 Big Five boundaries removed | 25% | no |
| H3 | H1 holds after removing the rank effect of both bins' durations | 30% | no |
| H4 | H1 robust to PhanDA uncertainty (>= 90% of draws give rho > 0) | 40% | YES |
| H5a | Warming boundaries: abs(dT) predicts mu3t | 35% | no |
| H5b | Cooling boundaries: abs(dT) predicts mu3t | 25% | no |
| H5c | Warming correlation stronger than cooling (Fisher z) | 15% | no |
| H6a | Threshold: all 5 Big Five events have abs(dT) >= 5.2 C and no other boundary does (temperature only; disclosed before the result) | 3% | no |
| H6b | Big Five boundaries rank high by abs(dT) (AUC > 0.5, p < 0.05) (temperature only; disclosed before the result) | 35% | YES |
| H7a | Temperature level (GMST of the bin) predicts mu3t | 30% | no |
| H7b | abs(dT) predicts mu3t beyond the level | 35% | YES |
| G | Gap-filler gives the same sign as three-timer for H1 and H2 | 80% | YES |

Brier score of our registered priors: 0.183 over 14 rows (a 50% forecaster scores 0.250);
0.179 over the 12 rows that involve extinction rates.

## Key numbers

| Test | n | rho | p_perm | p_shift |
|---|---|---|---|---|
| H1 all boundaries | 78 | 0.207 | 0.033 | 0.041 |
| H2 without Big Five | 72 | 0.143 | 0.114 | 0.118 |
| H3 duration partial | 78 | 0.126 | 0.136 | |
| H3 duration partial, without Big Five | 72 | 0.020 | 0.429 | |
| H5a warming boundaries | 34 | 0.208 | 0.121 | |
| H5b cooling boundaries | 44 | 0.182 | 0.115 | |
| H7a temperature level | 78 | 0.155 | 0.087 | 0.176 |
| H7b change given level | 78 | 0.220 | 0.030 | |
| Gap-filler rate, all | 78 | 0.237 | 0.019 | 0.027 |
| Gap-filler rate, without Big Five | 72 | 0.168 | 0.079 | 0.074 |

- H4: rho 5-95% across PhanDA draws 0.007 to 0.321; 95.8% positive.
- H5c: warming minus cooling z = 0.11, p = 0.45.
- Sanity: all five Big Five events are among the ten highest extinction rates
  (Changhsingian, Tremadocian, Rhaetian, Maastrichtian, Famennian, Frasnian, Hirnantian,
  Visean, Darriwilian, Induan).

## What this means, and its limits

- On an independent, data-assimilated global temperature record, the temperature-extinction
  link is real across the whole Phanerozoic but is carried by the mass extinctions. For
  background extinction we find no significant link at stage resolution, and no sign that
  warming is worse than cooling.
- The 5.2 C figure is a property of Song et al.'s within-bin maxima of tropical proxies; it
  does not transfer to PhanDA's stage-mean global temperatures, where several 6-9 C changes
  passed without a mass extinction. Anyone quoting 5.2 C as a threshold for today should
  know that it depends on how the temperature series is built.
- Limits: PhanDA and Song et al. share part of their proxy data (partly independent test).
  PhanDA's stage means smooth the sub-million-year transients that Song et al. measured, so a
  null for background extinction may partly reflect that smoothing. 72-78 boundaries give
  modest power: a true rho of 0.15 would often be missed.
- The extinction series uses the same PBDB marine download and three-timer code as
  [extinction-cycle-test](https://github.com/munimtechnologies/extinction-cycle-test)
  PREREG-1, whose per-stage series was public before this job ran; the rates were never
  paired with PhanDA temperatures before this job.
- Hand checks: all 79 PhanDA GMST_50 values in `boundaries.csv` (78 `T_bin` values plus the
  Piacenzian as the Zanclean's `T_next`) match the hand-read table in the grader within
  0.011 C. dT for the Changhsingian = 28.499 - 18.019 = 10.480 C (file 10.480); Katian
  20.802 - 27.848 = -7.045 C (file -7.045).

## Files

`results/`: `results.json` (every test), `boundaries.csv` (per boundary: extinction rates,
PhanDA temperatures with 16/50/84th percentiles, dT, durations, Big Five flag), `bins.csv`,
`provenance.json` (PBDB and PhanDA URLs, md5 sums, times).

Code Apache-2.0; text CC BY 4.0. MV2, Munim, Inc., mv2@munimtech.com
