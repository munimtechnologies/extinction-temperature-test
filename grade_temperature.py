# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""grade_temperature.py: the mechanical reading of extinction-temperature-test.

Written and published 2026-11-17 (about 17:30 UTC), BEFORE the analysis job
(extinction-temperature-v1) had produced any result and before any PBDB
extinction rate was paired with any PhanDA temperature.

It turns the job's results.json into the graded prior table of PREREG.md
("Hypotheses, tests and honest priors") and the registered decision label,
row by row, exactly as registered. No new statistic is computed here; it only
reads numbers the job wrote. Two reading choices are fixed now, in writing:

  * H1s ("H1 survives autocorrelation") is graded TRUE only if rho > 0 AND
    p_shift < 0.05, as for every other row of the table. The job's own
    qualifier `autocorrelation_robust_H1` uses p_shift alone; both are shown.
  * G ("gap-filler gives the same sign as three-timer for H1 and H2") is TRUE
    only if both rhos are finite and non-zero and their signs match for H1
    and for H2.

H6a and H6b depend only on PhanDA temperatures and the Big Five labels, not on
extinction rates (except that rows with an undefined three-timer rate are
dropped). They were computed by hand from public PhanDA before the result and
disclosed in PRE-RESULT-NOTE-1.md. This script marks them as such.

Optional: --boundaries boundaries.csv checks the job's PhanDA medians against
the GMST_50 values read by hand on 2026-11-17 (PRE-RESULT-NOTE-1.md, 2
decimals, tolerance 0.011 C). A mismatch is reported, never silently fixed;
it does not change the label but stops the result from being graded until a
person has looked (status CHECK_DATA).

    python grade_temperature.py --results results.json \
        [--boundaries boundaries.csv] --out-dir .

Outputs: graded.json and GRADED.md. Made with AI agents working in MV2's name.
"""
import argparse
import csv
import json
import math
import os

ALPHA = 0.05

# (id, prior that it holds, short text) -- copied from PREREG.md.
PRIORS = [
    ("S", 0.85, "Sanity: >= 4 of the Big Five among the 10 highest-mu3t boundaries, and >= 60 boundaries"),
    ("H1", 0.40, "|dT| predicts mu3t across all boundaries"),
    ("H1s", 0.30, "H1 survives autocorrelation (circular shifts)"),
    ("H2", 0.25, "H1 holds with the 6 Big Five boundaries removed"),
    ("H3", 0.30, "H1 holds after removing the rank effect of both bins' durations"),
    ("H4", 0.40, "H1 robust to PhanDA uncertainty (>= 90% of draws give rho > 0)"),
    ("H5a", 0.35, "Warming boundaries: |dT| predicts mu3t"),
    ("H5b", 0.25, "Cooling boundaries: |dT| predicts mu3t"),
    ("H5c", 0.15, "Warming correlation stronger than cooling (Fisher z)"),
    ("H6a", 0.03, "Threshold: all 5 Big Five events have |dT| >= 5.2 C and no other boundary does"),
    ("H6b", 0.35, "Big Five boundaries rank high by |dT| (AUC > 0.5, p < 0.05)"),
    ("H7a", 0.30, "Temperature level (GMST of the bin) predicts mu3t"),
    ("H7b", 0.35, "|dT| predicts mu3t beyond the level"),
    ("G", 0.80, "Gap-filler gives the same sign as three-timer for H1 and H2"),
]
TEMPERATURE_ONLY = ("H6a", "H6b")
LABELS = ("DATA_LIMITED", "NOT_REPLICATED", "REPLICATED_ONLY_WITH_MASS_EXTINCTIONS",
          "REPLICATED")

# PhanDA GMST_50 (deg C), read by hand from PhanDA_GMSTandCO2_percentiles.csv on
# 2026-11-17, rounded to 2 decimals; the 79 temperature bins, oldest first.
HAND_GMST50 = {
    "Tremadocian": 32.10, "Floian": 31.53, "Dapingian": 30.90, "Darriwilian": 31.00,
    "Sandbian": 27.31, "Katian": 27.85, "Hirnantian": 20.80, "Aeronian/Rhuddanian": 19.48,
    "Telychian": 24.36, "Sheinwoodian": 19.90, "Homerian": 20.97, "Gorstian": 25.60,
    "Ludfordian": 34.05, "Pridoli": 33.25, "Pragian/Lochkovian": 28.81, "Emsian": 25.29,
    "Eifelian": 22.49, "Givetian": 25.67, "Frasnian": 27.91, "Famennian": 23.76,
    "Tournaisian": 21.75, "Visean": 12.83, "Serpukhovian": 12.02, "Bashkirian": 14.08,
    "Moscovian": 18.09, "Kasimovian": 14.46, "Gzhelian": 16.55, "Asselian": 20.25,
    "Sakmarian": 20.80, "Artinskian": 20.76, "Kungurian": 18.70, "Roadian": 13.50,
    "Capitanian/Wordian": 18.68, "Wuchiapingian": 22.19, "Changhsingian": 18.02,
    "Induan": 28.50, "Anisian/Olenekian": 25.82, "Ladinian": 19.18, "Carnian": 23.43,
    "Norian": 25.02, "Rhaetian": 22.32, "Sinemurian/Hettangian": 25.06,
    "Pliensbachian": 24.10, "Toarcian": 26.37, "Aalenian": 23.20, "Bajocian": 21.39,
    "Bathonian": 17.77, "Callovian": 20.11, "Oxfordian": 23.88, "Kimmeridgian": 24.47,
    "Tithonian": 27.27, "Berriasian": 27.00, "Valanginian": 28.25, "Hauterivian": 27.27,
    "Barremian": 30.52, "Aptian": 32.00, "Albian": 26.01, "Cenomanian": 34.14,
    "Turonian": 35.56, "Coniacian": 32.26, "Santonian": 31.42, "Campanian": 28.13,
    "Maastrichtian": 24.51, "Selandian/Danian": 28.85, "Thanetian": 25.60,
    "Ypresian": 32.28, "Lutetian": 26.99, "Bartonian": 28.03, "Priabonian": 26.08,
    "Rupelian": 23.40, "Chattian": 23.44, "Aquitanian": 22.07, "Burdigalian": 24.06,
    "Langhian": 23.78, "Serravallian": 23.60, "Tortonian": 20.62, "Messinian": 16.94,
    "Zanclean": 16.27, "Piacenzian": 15.63,
}
HAND_TOL = 0.011


def _f(x):
    """float or NaN."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def _finite(x):
    return not math.isnan(_f(x)) and not math.isinf(_f(x))


def _pos_sig(d, p_key="p_perm"):
    """rho > 0 and p < ALPHA; False if missing or NaN."""
    if not isinstance(d, dict):
        return False
    rho, p = _f(d.get("rho")), _f(d.get(p_key))
    return _finite(rho) and _finite(p) and rho > 0 and p < ALPHA


def _sign(x):
    x = _f(x)
    if not _finite(x) or x == 0:
        return 0
    return 1 if x > 0 else -1


def grade_rows(a):
    """{id: True/False/None} from results.json's "analysis" block. None means
    the quantity was not computed (e.g. too few boundaries)."""
    g = {k: None for k, _, _ in PRIORS}
    if not isinstance(a, dict):
        return g
    san = a.get("sanity")
    if isinstance(san, dict) and "passed" in san:
        g["S"] = bool(san["passed"])
    if "H1_all" not in a:          # analyze() stopped at MIN_RUN
        return g
    h1 = a["H1_all"]
    g["H1"] = _pos_sig(h1)
    g["H1s"] = bool(_finite(h1.get("rho")) and _f(h1.get("rho")) > 0
                    and _finite(h1.get("p_shift")) and _f(h1.get("p_shift")) < ALPHA)
    g["H2"] = _pos_sig(a.get("H2_without_big_five"))
    g["H3"] = _pos_sig(a.get("H3_duration_partial"))
    h4 = a.get("H4_uncertainty") or {}
    g["H4"] = bool(_finite(h4.get("frac_rho_pos")) and _f(h4.get("frac_rho_pos")) >= 0.9)
    g["H5a"] = _pos_sig(a.get("H5_warming"))
    g["H5b"] = _pos_sig(a.get("H5_cooling"))
    fd = a.get("H5_warming_minus_cooling") or {}
    g["H5c"] = bool(_finite(fd.get("p_one_sided")) and _f(fd.get("p_one_sided")) < ALPHA
                    and _finite(fd.get("z")) and _f(fd.get("z")) > 0)
    th = a.get("H6_threshold") or {}
    if "big_five_events_above" in th:
        g["H6a"] = bool(int(th["big_five_events_above"]) == 5
                        and len(th.get("other_boundaries_above", [None])) == 0)
    au = a.get("H6_big_five_auc") or {}
    g["H6b"] = bool(_finite(au.get("auc")) and _f(au.get("auc")) > 0.5
                    and _finite(au.get("p_perm")) and _f(au.get("p_perm")) < ALPHA)
    g["H7a"] = _pos_sig(a.get("H7_level"))
    g["H7b"] = _pos_sig(a.get("H7_change_given_level"))
    gf = a.get("gap_filler") or {}
    s1 = (_sign(h1.get("rho")), _sign((gf.get("H1_all") or {}).get("rho")))
    s2 = (_sign((a.get("H2_without_big_five") or {}).get("rho")),
          _sign((gf.get("H2_without_big_five") or {}).get("rho")))
    g["G"] = bool(0 not in s1 and 0 not in s2 and s1[0] == s1[1] and s2[0] == s2[1])
    return g


def label_from(g):
    """PREREG decision rule from the graded rows."""
    if not g.get("S"):
        return "DATA_LIMITED"
    if not g.get("H1"):
        return "NOT_REPLICATED"
    if g.get("H2"):
        return "REPLICATED"
    return "REPLICATED_ONLY_WITH_MASS_EXTINCTIONS"


def brier(g):
    """Mean squared error of the registered priors over graded rows, and the
    same for a know-nothing 0.5 forecaster."""
    sq, sq05, n = 0.0, 0.0, 0
    for k, p, _ in PRIORS:
        if g.get(k) is None:
            continue
        o = 1.0 if g[k] else 0.0
        sq += (p - o) ** 2
        sq05 += 0.25
        n += 1
    if n == 0:
        return {"n": 0, "brier": float("nan"), "brier_coinflip": float("nan")}
    return {"n": n, "brier": sq / n, "brier_coinflip": sq05 / n}


def check_boundaries(rows):
    """Compare boundaries.csv (list of dicts) with HAND_GMST50. Returns a
    dict with the number of bins checked and every mismatch."""
    seen, bad = {}, []
    for r in rows:
        name = (r.get("bin") or "").strip()
        if name:
            seen[name] = _f(r.get("T_bin"))
    # The last temperature bin (Piacenzian) only appears as T_next of Zanclean.
    for r in rows:
        if (r.get("bin") or "").strip() == "Zanclean":
            seen.setdefault("Piacenzian", _f(r.get("T_next")))
    for name, v in sorted(seen.items()):
        if name not in HAND_GMST50:
            bad.append({"bin": name, "job": v, "hand": None})
        elif not (_finite(v) and abs(v - HAND_GMST50[name]) <= HAND_TOL):
            bad.append({"bin": name, "job": v, "hand": HAND_GMST50[name]})
    missing = sorted(set(HAND_GMST50) - set(seen))
    return {"n_checked": len(seen), "mismatches": bad, "missing_bins": missing,
            "passed": not bad}


def grade(results, boundary_rows=None):
    """Everything graded.json holds."""
    out = {"source": "grade_temperature.py (published before results)",
           "job_status": results.get("status"), "label": None, "rows": [],
           "temperature_only_rows": list(TEMPERATURE_ONLY)}
    if results.get("status") != "ok" or "analysis" not in results:
        out["status"] = "JOB_FAILED"
        out["label"] = None
        return out
    a = results["analysis"]
    g = grade_rows(a)
    out["label"] = label_from(g)
    out["label_from_job"] = results.get("decision", a.get("decision"))
    out["label_agrees_with_job"] = out["label"] == out["label_from_job"]
    q = a.get("qualifiers") or {}
    out["qualifiers"] = {
        "autocorrelation_robust_H1 (job: p_shift only)": q.get("autocorrelation_robust_H1"),
        "duration_robust_H1": q.get("duration_robust_H1"),
        "uncertainty_robust_H1": q.get("uncertainty_robust_H1")}
    for k, p, text in PRIORS:
        out["rows"].append({"id": k, "prior": p, "prediction": text, "holds": g[k],
                            "temperature_only": k in TEMPERATURE_ONLY})
    out["brier"] = brier(g)
    out["brier_excluding_temperature_only"] = brier(
        {k: v for k, v in g.items() if k not in TEMPERATURE_ONLY})
    status = "GRADED"
    if not g.get("S"):
        status = "GRADED_SANITY_FAILED"
    if not out["label_agrees_with_job"]:
        status = "CHECK_CODE"
    if boundary_rows is not None:
        out["phanda_check"] = check_boundaries(boundary_rows)
        if not out["phanda_check"]["passed"]:
            status = "CHECK_DATA"
    out["status"] = status
    out["n_boundaries"] = a.get("n")
    return out


def _fmt(x, nd=3):
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return "yes" if x else "no"
    v = _f(x)
    if not _finite(v):
        return "n/a"
    return f"{v:.{nd}f}"


def render_md(gr, results):
    a = results.get("analysis", {}) if isinstance(results, dict) else {}
    L = ["# Graded result: extinction-temperature-test", "",
         "Produced mechanically by `grade_temperature.py`, which was published before the",
         "result existed. Made with AI agents working in MV2's name.", "",
         f"- Status: **{gr.get('status')}**",
         f"- Decision label (PREREG rule): **{gr.get('label')}**"]
    if gr.get("status") == "JOB_FAILED":
        L += ["", f"The job did not finish its analysis (status `{gr.get('job_status')}`).",
              "Nothing is graded."]
        return "\n".join(L) + "\n"
    L += [f"- Label written by the job: {gr.get('label_from_job')} "
          f"(agrees: {_fmt(gr.get('label_agrees_with_job'))})",
          f"- Boundaries with a three-timer rate: {gr.get('n_boundaries')}", ""]
    if gr.get("status") == "GRADED_SANITY_FAILED":
        L += ["The sanity gate failed, so the label is DATA_LIMITED and the rows below are",
              "reported but not interpreted.", ""]
    L += ["| # | Prediction | Prior | Holds? |", "|---|---|---|---|"]
    for r in gr["rows"]:
        tag = " (temperature only; computed before the result, PRE-RESULT-NOTE-1.md)" \
            if r["temperature_only"] else ""
        L.append(f"| {r['id']} | {r['prediction']}{tag} | {int(round(100 * r['prior']))}% | "
                 f"{'n/a' if r['holds'] is None else ('YES' if r['holds'] else 'no')} |")
    b, bx = gr["brier"], gr["brier_excluding_temperature_only"]
    L += ["", f"Brier score of the registered priors: {_fmt(b['brier'])} over {b['n']} rows "
          f"(a 50% forecaster scores {_fmt(b['brier_coinflip'])}); without the two "
          f"temperature-only rows: {_fmt(bx['brier'])} over {bx['n']}.", "",
          "## Key numbers", ""]
    for key, name in (("H1_all", "H1 all boundaries"),
                      ("H2_without_big_five", "H2 without Big Five"),
                      ("H3_duration_partial", "H3 duration partial"),
                      ("H5_warming", "H5a warming"), ("H5_cooling", "H5b cooling"),
                      ("H7_level", "H7a level"), ("H7_change_given_level", "H7b change | level")):
        d = a.get(key) or {}
        L.append(f"- {name}: rho {_fmt(d.get('rho'))}, n {d.get('n', 'n/a')}, "
                 f"p_perm {_fmt(d.get('p_perm'), 4)}"
                 + (f", p_shift {_fmt(d.get('p_shift'), 4)}" if "p_shift" in d else ""))
    h4 = a.get("H4_uncertainty") or {}
    L.append(f"- H4: fraction of PhanDA draws with rho > 0: {_fmt(h4.get('frac_rho_pos'))} "
             f"(5-95%: {_fmt(h4.get('rho_p05'))} to {_fmt(h4.get('rho_p95'))})")
    au = a.get("H6_big_five_auc") or {}
    th = a.get("H6_threshold") or {}
    L.append(f"- H6b: AUC {_fmt(au.get('auc'))}, p_perm {_fmt(au.get('p_perm'), 4)}")
    L.append(f"- H6a: Big Five events at or above 5.2 C: {th.get('big_five_events_above', 'n/a')} "
             f"of 5; other boundaries at or above: {', '.join(th.get('other_boundaries_above', [])) or 'none'}")
    gf = a.get("gap_filler") or {}
    for key in ("H1_all", "H2_without_big_five"):
        d = gf.get(key) or {}
        L.append(f"- gap-filler {key}: rho {_fmt(d.get('rho'))}, p_perm {_fmt(d.get('p_perm'), 4)}")
    L += ["", "## Qualifiers (never change the label)", ""]
    for k, v in (gr.get("qualifiers") or {}).items():
        L.append(f"- {k}: {_fmt(v)}")
    if "phanda_check" in gr:
        pc = gr["phanda_check"]
        L += ["", "## PhanDA hand check", "",
              f"{pc['n_checked']} bins compared with the hand-read GMST_50 values; "
              f"mismatches: {len(pc['mismatches'])}; missing: {len(pc['missing_bins'])}."]
        for m in pc["mismatches"]:
            L.append(f"- {m['bin']}: job {_fmt(m['job'], 2)}, hand {_fmt(m['hand'], 2)}")
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--boundaries", default=None)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args(argv)
    with open(args.results, encoding="utf-8") as f:
        results = json.load(f)
    rows = None
    if args.boundaries and os.path.exists(args.boundaries):
        with open(args.boundaries, newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
    gr = grade(results, rows)
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "graded.json"), "w", encoding="utf-8") as f:
        json.dump(gr, f, indent=1, default=str)
    with open(os.path.join(args.out_dir, "GRADED.md"), "w", encoding="utf-8") as f:
        f.write(render_md(gr, results))
    print(json.dumps({"status": gr["status"], "label": gr["label"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
