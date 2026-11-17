# PREREG: does temperature change predict marine extinction on an independent temperature history?

Registered 2026-11-17 by MV2 of Munim, Inc., before the analysis job ran and before any
PBDB extinction rate was paired with any PhanDA temperature. Code: `thermo.py`,
`job_thermo.py`, tests `test_thermo.py`. Made with AI agents working in MV2's name.

Disclosure of what was seen before registering: to check the file format we read the
PhanDA percentile CSV's header, its 85 slice names and ages, and the GMST values of its
first 15 rows (Holocene to Rupelian) and last 5 rows (Sandbian to Tremadocian). No
extinction rate from this project's data has been computed or looked at.

## The claim under test

Song, Kemp, Tian, Chu, Song & Dai (2021), *Thresholds of temperature change for mass
extinctions*, Nature Communications 12:4694 (doi:10.1038/s41467-021-25019-2):

- Across 45 Phanerozoic time bins (Katian-Hirnantian to Aquitanian-Burdigalian), the
  largest temperature change within a bin correlates with marine genus extinction
  (Spearman, dT vs gap-filler rate: rho = 0.63, n = 45; vs three-timer: rho = 0.59), and
  still does with the Big Five removed (rho = 0.54 and 0.47, n = 40).
- Extinction rises more steeply with warming than with cooling (warming rho = 0.76,
  n = 20; cooling rho = -0.53 against signed dT, n = 25; the difference was not significant).
- Mass extinctions are linked to thresholds: change > 5.2 C at > 10 C/Myr. Only the Big
  Five exceed them. The 5.2 C threshold is now quoted for present-day warming.

Their temperatures come from a selected compilation of tropical sea-surface proxies
(mostly d18O), read at high resolution inside each bin.

## Why this test is new

PhanDA (Judd et al. 2024, Science 385:eadk3705; data CC BY 4.0, github.com/EJJudd/PhanDA)
is a global mean surface temperature (GMST) history for 485 Myr. It was built by data
assimilation of the PhanSST sea-surface proxy database with HadCM3L climate model
priors, and published three years after Song et al. We found no study that pairs PhanDA
with PBDB extinction rates (searches on 2026-11-17). Malanoski et al. (2024, Science
383:1130; 2026, Science 391:285) model the extinction risk of individual taxa with climate
model output. They do not test Song et al.'s boundary-level correlations or thresholds.
PhanDA and Song et al. share part of their proxy data, so this is a partly independent test
(different selection, method, spatial coverage and timescale), not a fully independent one.

## What is measured

- **Bins.** PhanDA's 85 time slices. The 5 Quaternary slices and the 0.3-Myr PETM slice
  are dropped, which leaves 79 temperature bins from the Tremadocian to the Piacenzian.
  PBDB occurrences are binned into GTS2020 stages (the timescale PBDB serves). Stages are
  then merged by name where PhanDA merged them (Rhuddanian+Aeronian, Lochkovian+Pragian,
  Wordian+Capitanian, Olenekian+Anisian, Hettangian+Sinemurian, Danian+Selandian). The 10
  Cambrian stages and the Quaternary are kept as context bins with no temperature.
- **Extinction.** PBDB marine genus occurrences (Porifera, Cnidaria, Bryozoa, Brachiopoda,
  Mollusca, Echinodermata, Trilobita, Graptolithina; envtype=marine, genus certain,
  regular preservation), binned with the same rules as extinction-cycle-test. The primary
  rate is the three-timer extinction rate at the top of each bin (Alroy 2008). The
  secondary rate is gap-filler (Alroy 2014).
- **Temperature change at the top of bin j.** dT = GMST50(j+1) - GMST50(j), from PhanDA's
  median per slice. |dT| is its magnitude. That gives 78 boundaries, from the top of the
  Tremadocian to the top of the Zanclean.
- **Big Five boundaries.** The tops of the Katian, Hirnantian, Frasnian, Changhsingian,
  Rhaetian and Maastrichtian (the end-Ordovician has two pulses).

This is a stage-scale test. PhanDA gives one mean temperature per stage, so the
sub-million-year transients Song et al. used are smoothed. Song et al.'s rate threshold
(> 10 C/Myr) cannot be tested at this resolution, and is not.

## Hypotheses, tests and honest priors

All correlations are Spearman. p_perm is one-sided (rho > 0) from 20,000 permutations.
p_shift is one-sided from all circular shifts of 3 or more steps, which keep each series'
autocorrelation. Seed 20261117.

| # | Prediction (Song et al. carried to PhanDA) | Test | Prior that it holds |
|---|---|---|---|
| S | Sanity: >= 4 of the Big Five among the 10 highest-mu3t boundaries, and >= 60 boundaries | gate | 85% |
| H1 | |dT| predicts mu3t across all boundaries | rho > 0, p_perm < 0.05 | 40% |
| H1s | H1 survives autocorrelation | p_shift < 0.05 | 30% |
| H2 | H1 holds with the 6 Big Five boundaries removed | rho > 0, p_perm < 0.05 | 25% |
| H3 | H1 holds after removing the rank effect of both bins' durations | partial rho > 0, p < 0.05 | 30% |
| H4 | H1 is robust to PhanDA uncertainty (independent two-piece normal draws from the 16/50/84th percentiles, 2,000 draws) | >= 90% of draws give rho > 0 | 40% |
| H5a | Warming boundaries: |dT| predicts mu3t | rho > 0, p_perm < 0.05 | 35% |
| H5b | Cooling boundaries: |dT| predicts mu3t | rho > 0, p_perm < 0.05 | 25% |
| H5c | Warming correlation stronger than cooling | Fisher z, one-sided p < 0.05 | 15% |
| H6a | Threshold: all 5 Big Five events have |dT| >= 5.2 C AND no other boundary does | exact | 3% |
| H6b | Big Five boundaries rank high by |dT| | AUC > 0.5, p_perm < 0.05 | 35% |
| H7a | Temperature level (GMST of the bin) predicts mu3t (Mayhew et al. 2008: warm = more extinction) | rho > 0, p_perm < 0.05 | 30% |
| H7b | |dT| predicts mu3t beyond the level | partial rho > 0, p < 0.05 | 35% |
| G | Gap-filler gives the same sign as three-timer for H1 and H2 | sign match | 80% |

## Decision rule (fixed now)

- `DATA_LIMITED` if the sanity gate fails.
- `NOT_REPLICATED` if H1 fails.
- `REPLICATED_ONLY_WITH_MASS_EXTINCTIONS` if H1 holds and H2 fails.
- `REPLICATED` if H1 and H2 both hold.

Qualifiers reported alongside, never changing the label: autocorrelation-robust (H1s),
duration-robust (H3), uncertainty-robust (H4). The threshold result (H6a) is reported on its
own. Overall expectation: `REPLICATED_ONLY_WITH_MASS_EXTINCTIONS` or `NOT_REPLICATED` (about
two in three), because stage means smooth the short warming and cooling pulses.

## How the result will be read

- If it replicates beyond the Big Five on PhanDA, that is the first confirmation of the
  temperature-change and extinction link on an independent, global, data-assimilated
  temperature history.
- If it fails, the link at Phanerozoic scale rests on short transients, or on the choice of
  proxy records, rather than on stage-to-stage climate change. That narrows what the 5.2 C
  figure can mean.

Either result is reported in RESULTS.md with every number, whatever the outcome.

## Known limits (stated before the result)

- Song et al. used 45 coarser bins and within-bin maxima. This test uses about 78 stage
  boundaries and between-stage changes of means. It tests the claim's stage-scale form.
- PhanDA draws are treated as independent between slices. Its ensemble is correlated in
  time, so H4 overstates the noise in dT (a conservative test).
- PBDB is queried live, so its content changes over time. File md5s are saved in provenance.json.
- The group list covers most but not all of Song et al.'s taxa (they used all metazoans
  except arachnids, insects, ostracods and tetrapods).
