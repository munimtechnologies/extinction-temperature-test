# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Offline tests for grade_temperature.py (no data, no network).
Run: python -m unittest -v test_grade_temperature"""
import copy
import json
import os
import tempfile
import unittest

import grade_temperature as gt

NAN = float("nan")


def sig(rho=0.5, p=0.001, **kw):
    d = {"n": 78, "rho": rho, "p_perm": p}
    d.update(kw)
    return d


def make_results():
    """A results.json in the shape job_thermo.py writes, where every
    registered prediction holds."""
    a = {
        "n": 78,
        "sanity": {"top10_mu3t": [], "big_five_events_in_top10": 5, "n_ok": True,
                   "passed": True},
        "H1_all": sig(0.6, 0.0001, p_shift=0.01, n_shifts=73),
        "H2_without_big_five": sig(0.4, 0.002, p_shift=0.03, n_shifts=67),
        "H3_duration_partial": sig(0.5, 0.001),
        "H3_duration_partial_without_big_five": sig(0.3, 0.02),
        "H4_uncertainty": {"n_draws": 2000, "rho_p05": 0.2, "rho_median": 0.5,
                           "rho_p95": 0.7, "frac_rho_pos": 0.99},
        "H5_warming": sig(0.7, 0.001),
        "H5_cooling": sig(0.3, 0.04),
        "H5_warming_minus_cooling": {"z": 2.5, "p_one_sided": 0.006},
        "H6_threshold": {"boundaries_at_or_above": [], "big_five_events_above": 5,
                         "other_boundaries_above": [], "big_five_abs_dT": {},
                         "max_abs_dT_other": 4.0},
        "H6_big_five_auc": {"auc": 0.9, "p_perm": 0.001},
        "H7_level": sig(0.3, 0.01, p_shift=0.04),
        "H7_change_given_level": sig(0.5, 0.001),
        "gap_filler": {"H1_all": sig(0.5, 0.001, p_shift=0.02),
                       "H2_without_big_five": sig(0.3, 0.01, p_shift=0.04)},
        "decision": "REPLICATED",
        "qualifiers": {"autocorrelation_robust_H1": True, "duration_robust_H1": True,
                       "uncertainty_robust_H1": True},
    }
    return {"prereg": "PREREG.md", "status": "ok", "analysis": a, "decision": "REPLICATED"}


def holds(gr):
    return {r["id"]: r["holds"] for r in gr["rows"]}


def hand_rows(shift=None):
    names = list(gt.HAND_GMST50)
    rows = []
    for j in range(len(names) - 1):
        t0, t1 = gt.HAND_GMST50[names[j]], gt.HAND_GMST50[names[j + 1]]
        if shift and names[j] == shift:
            t0 += 0.5
        rows.append({"bin": names[j], "T_bin": repr(t0 + 0.004), "T_next": repr(t1 - 0.004)})
    return rows


class Table(unittest.TestCase):
    def test_priors_match_prereg(self):
        self.assertEqual([k for k, _, _ in gt.PRIORS],
                         ["S", "H1", "H1s", "H2", "H3", "H4", "H5a", "H5b", "H5c", "H6a",
                          "H6b", "H7a", "H7b", "G"])
        self.assertEqual([int(round(100 * p)) for _, p, _ in gt.PRIORS],
                         [85, 40, 30, 25, 30, 40, 35, 25, 15, 3, 35, 30, 35, 80])

    def test_hand_table(self):
        self.assertEqual(len(gt.HAND_GMST50), 79)
        names = list(gt.HAND_GMST50)
        self.assertEqual(names[0], "Tremadocian")
        self.assertEqual(names[-1], "Piacenzian")
        self.assertNotIn("PETM", gt.HAND_GMST50)


class Grading(unittest.TestCase):
    def test_all_hold(self):
        gr = gt.grade(make_results())
        self.assertEqual(gr["status"], "GRADED")
        self.assertEqual(gr["label"], "REPLICATED")
        self.assertTrue(gr["label_agrees_with_job"])
        self.assertTrue(all(v is True for v in holds(gr).values()))
        # Brier of the registered priors when every row holds: sum((1-p)^2)/14,
        # worked by hand = 6.3084 / 14.
        self.assertAlmostEqual(gr["brier"]["brier"], 6.3084 / 14, places=6)
        self.assertEqual(gr["brier"]["n"], 14)
        self.assertEqual(gr["brier_excluding_temperature_only"]["n"], 12)
        self.assertAlmostEqual(gr["brier"]["brier_coinflip"], 0.25, places=12)

    def test_only_with_mass_extinctions(self):
        r = make_results()
        r["analysis"]["H2_without_big_five"] = sig(0.1, 0.2, p_shift=0.3)
        r["decision"] = r["analysis"]["decision"] = "REPLICATED_ONLY_WITH_MASS_EXTINCTIONS"
        gr = gt.grade(r)
        self.assertEqual(gr["label"], "REPLICATED_ONLY_WITH_MASS_EXTINCTIONS")
        self.assertEqual(gr["status"], "GRADED")
        self.assertFalse(holds(gr)["H2"])

    def test_not_replicated_negative_rho(self):
        r = make_results()
        r["analysis"]["H1_all"] = sig(-0.4, 0.99, p_shift=0.01)
        r["decision"] = "NOT_REPLICATED"
        gr = gt.grade(r)
        self.assertEqual(gr["label"], "NOT_REPLICATED")
        h = holds(gr)
        self.assertFalse(h["H1"])
        self.assertFalse(h["H1s"])        # rho must be > 0 as well
        self.assertFalse(h["G"])          # gap-filler H1 sign differs

    def test_sanity_failed(self):
        r = make_results()
        r["analysis"]["sanity"]["passed"] = False
        r["decision"] = "DATA_LIMITED"
        gr = gt.grade(r)
        self.assertEqual(gr["label"], "DATA_LIMITED")
        self.assertEqual(gr["status"], "GRADED_SANITY_FAILED")
        self.assertFalse(holds(gr)["S"])

    def test_too_few_boundaries(self):
        r = {"status": "ok", "decision": "DATA_LIMITED",
             "analysis": {"n": 5, "boundaries": [], "sanity": {"passed": False, "n": 5},
                          "decision": "DATA_LIMITED"}}
        gr = gt.grade(r)
        self.assertEqual(gr["label"], "DATA_LIMITED")
        h = holds(gr)
        self.assertFalse(h["S"])
        self.assertTrue(all(h[k] is None for k in h if k != "S"))
        self.assertEqual(gr["brier"]["n"], 1)
        self.assertAlmostEqual(gr["brier"]["brier"], 0.85 ** 2, places=12)

    def test_job_failed(self):
        gr = gt.grade({"status": "pbdb_incomplete: missing Trilobita"})
        self.assertEqual(gr["status"], "JOB_FAILED")
        self.assertIsNone(gr["label"])
        self.assertIn("Nothing is graded", gt.render_md(gr, {}))

    def test_label_disagreement_flags_code(self):
        r = make_results()
        r["decision"] = "NOT_REPLICATED"
        self.assertEqual(gt.grade(r)["status"], "CHECK_CODE")

    def test_nan_and_threshold_rows(self):
        r = make_results()
        a = r["analysis"]
        a["H5_cooling"] = {"n": 3, "rho": NAN, "p_perm": 1.0}
        a["H5_warming_minus_cooling"] = {"z": NAN, "p_one_sided": 1.0}
        a["H6_threshold"]["big_five_events_above"] = 2
        a["H6_threshold"]["other_boundaries_above"] = ["Gorstian", "Tournaisian"]
        a["H4_uncertainty"]["frac_rho_pos"] = 0.85
        a["gap_filler"]["H2_without_big_five"]["rho"] = NAN
        h = holds(gt.grade(r))
        for k in ("H5b", "H5c", "H6a", "H4", "G"):
            self.assertFalse(h[k], k)
        self.assertTrue(h["H6b"])

    def test_h6a_needs_no_other_boundary(self):
        r = make_results()
        r["analysis"]["H6_threshold"]["other_boundaries_above"] = ["Albian"]
        self.assertFalse(holds(gt.grade(r))["H6a"])

    def test_does_not_mutate_input(self):
        r = make_results()
        before = json.dumps(r, sort_keys=True)
        gt.grade(copy.deepcopy(r))
        gt.grade(r)
        self.assertEqual(json.dumps(r, sort_keys=True), before)


class PhanDACheck(unittest.TestCase):
    def test_hand_values_pass(self):
        pc = gt.check_boundaries(hand_rows())
        self.assertTrue(pc["passed"])
        self.assertEqual(pc["n_checked"], 79)
        self.assertEqual(pc["missing_bins"], [])

    def test_mismatch_flags_data(self):
        rows = hand_rows(shift="Induan")
        gr = gt.grade(make_results(), rows)
        self.assertEqual(gr["status"], "CHECK_DATA")
        self.assertEqual([m["bin"] for m in gr["phanda_check"]["mismatches"]], ["Induan"])
        self.assertEqual(gr["label"], "REPLICATED")   # the label itself is unchanged


class Output(unittest.TestCase):
    def test_main_writes_files(self):
        with tempfile.TemporaryDirectory() as d:
            rp = os.path.join(d, "results.json")
            with open(rp, "w", encoding="utf-8") as f:
                json.dump(make_results(), f)
            bp = os.path.join(d, "boundaries.csv")
            with open(bp, "w", encoding="utf-8", newline="") as f:
                f.write("bin,T_bin,T_next\n")
                for r in hand_rows():
                    f.write(f"{r['bin']},{r['T_bin']},{r['T_next']}\n")
            self.assertEqual(gt.main(["--results", rp, "--boundaries", bp, "--out-dir", d]), 0)
            with open(os.path.join(d, "graded.json"), encoding="utf-8") as f:
                g = json.load(f)
            self.assertEqual(g["status"], "GRADED")
            with open(os.path.join(d, "GRADED.md"), encoding="utf-8") as f:
                md = f.read()
            self.assertIn("REPLICATED", md)
            self.assertIn("| H6a |", md)
            self.assertIn("PhanDA hand check", md)


if __name__ == "__main__":
    unittest.main()
