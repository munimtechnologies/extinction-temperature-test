# Pre-result note 1 (2026-11-17, about 17:30 UTC): what the temperature record alone already says

Written by MV2 of Munim, Inc. (made with AI agents working in his name) **before the analysis
job `extinction-temperature-v1` had produced any output**. No PBDB extinction rate from this
project has been computed or seen. Nothing in `PREREG.md` is changed by this note.

## Why this note exists

Two of the registered predictions, **H6a** (the 5.2 C threshold) and **H6b** (Big Five
boundaries rank high by |dT|), depend only on PhanDA's temperatures and on which boundaries
are Big Five. They do not depend on extinction rates, except that a boundary is dropped if its
three-timer rate is undefined. PhanDA is public, so both could be worked out before the result.
We did that today, while writing the grading script, and disclose it here so that nobody reads
those two rows as a test that was passed or failed blind. `grade_temperature.py` marks them
"temperature only" and gives the Brier score of our priors with and without them.

## What we read

GMST_50 (median global mean surface temperature) for all 85 slices of
`PhanDA_GMSTandCO2_percentiles.csv` (github.com/EJJudd/PhanDA, Judd et al. 2024, CC BY 4.0),
read on 2026-11-17 through a web reader and rounded to 2 decimals. The 79 values used are in
`grade_temperature.py` (`HAND_GMST50`); the grader compares them with the job's
`boundaries.csv`, so a reading error on our side would show up there. Values near the
threshold were read again with all digits: Kungurian 18.6979, Roadian 13.5035,
Capitanian/Wordian 18.6812, Ypresian 32.2821, Lutetian 26.9858, Aptian 32.0025,
Albian 26.0115, Thanetian 25.6034.

## H6a: already false (registered prior 3%)

dT is the change in median GMST from one PhanDA stage to the next (PETM slice dropped, as
registered). Of the 78 boundaries, **9 reach 5.2 C or more**, and only **2 of them are Big Five
events**:

| Boundary (top of) | dT (C) | Big Five? |
|---|---|---|
| Changhsingian (end-Permian) | +10.48 | yes |
| Tournaisian | -8.92 | no |
| Gorstian | +8.45 | no |
| Albian | +8.13 | no |
| Katian (end-Ordovician, first pulse) | -7.05 | yes |
| Thanetian (PETM slice dropped) | +6.68 | no |
| Anisian/Olenekian | -6.64 | no |
| Aptian | -5.99 | no |
| Ypresian | -5.30 | no |

The other Big Five boundaries fall short: Hirnantian -1.32, Frasnian -4.15, Rhaetian +2.74,
Maastrichtian +4.34. Kungurian (-5.19) and Roadian (+5.18) are just under.

So at stage resolution on PhanDA, a 5.2 C change neither singles out the Big Five nor is it
needed for one. This does **not** refute Song et al. (2021): their threshold is for the largest
change *within* a time bin, from high-resolution tropical sea-surface records, and PREREG.md
already says stage means smooth such transients. It does say the 5.2 C figure cannot be carried
over to stage-to-stage changes in global mean temperature.

## H6b: probably true (registered prior 35%)

Ranking all 78 boundaries by |dT|, the six Big Five boundaries beat 72, 69, 59, 57, 38 and 21
of the 72 others (counted by hand), so AUC = 316/432 = **0.73**. The normal approximation to
the Mann-Whitney null (mean 216, SD 53.3) gives one-sided p of about **0.03**. The job's
permutation p (20,000 permutations, only boundaries with a defined rate) is the registered
number. If a few boundaries drop out, p can move across 0.05, so H6b stays open until the job
reports.

## What is still blind

Every prediction that involves an extinction rate (S, H1, H1s, H2, H3, H4, H5a, H5b, H5c, H7a,
H7b, G) and the decision label. Those are what the test is for.

Files: `grade_temperature.py`, `test_grade_temperature.py` (published with this note).
