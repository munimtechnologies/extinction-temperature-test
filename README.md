# extinction-temperature-test

**Does global temperature change predict mass extinction? A pre-registered test on an
independent 485-million-year temperature history.**

In 2021 Song et al. (Nature Communications 12:4694) reported that marine extinction rates
across the Phanerozoic track the size of temperature change, even without the Big Five
mass extinctions. They also found that only the Big Five crossed a change of 5.2 C. That
figure is now cited as a threshold for present-day warming.

In 2024 Judd et al. (Science 385:eadk3705) published PhanDA, a new global mean surface
temperature history for the last 485 million years, built by data assimilation of the
PhanSST proxy database with climate-model priors. This project asks whether Song et al.'s
result holds on PhanDA, using genus extinction rates from the Paleobiology Database on
PhanDA's own time slices.

- `PREREG.md`: hypotheses, tests, decision rule and honest priors, fixed before the job ran.
- `thermo.py`: binning, three-timer and gap-filler extinction rates, and the statistics.
- `job_thermo.py`: downloads PBDB and PhanDA and writes `results.json`, `boundaries.csv`,
  `bins.csv` and `provenance.json`.
- `test_thermo.py`: offline tests (`python -m unittest -v test_thermo`).
- `grade_temperature.py` (added 2026-11-17, before any result): turns `results.json` into the
  graded prior table and the registered label, mechanically, and checks the job's PhanDA
  values against a hand reading. Tests: `python -m unittest -v test_grade_temperature`.
- `PRE-RESULT-NOTE-1.md` (2026-11-17, before any result): the two predictions that depend only
  on temperature, worked out from public PhanDA. On PhanDA's stage-to-stage changes, 9 of 78
  boundaries reach 5.2 C and only 2 of them are Big Five events, so the 5.2 C threshold
  (prediction H6a) already fails at this resolution. Every prediction that involves
  extinction rates is still blind.

Status: **registered 2026-11-17; job queued.** Results will be posted in `RESULTS.md`
whatever they show.

Sister projects: [extinction-cycle-test](https://github.com/munimtechnologies/extinction-cycle-test)
(is extinction periodic?) and
[diversity-temperature-test](https://github.com/munimtechnologies/diversity-temperature-test)
(is a warmer Earth richer in life?).

## Data

- Paleobiology Database occurrences (paleobiodb.org, CC BY 4.0), queried live; md5s are
  recorded in the provenance file.
- PhanDA GMST percentiles, Judd et al. 2024 (github.com/EJJudd/PhanDA, CC BY 4.0).

## Licence and credit

Code: Apache-2.0. Text: CC BY 4.0. By MV2 of Munim, Inc. (mv2@munimtech.com); made with AI
agents working in his name. Corrections are welcome as GitHub issues.
