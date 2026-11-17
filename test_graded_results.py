# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""test_graded_results.py: runs grade_temperature.py (published before the result) on the
published job output in results/ and checks that it gives exactly the label and row grades
reported in RESULTS.md (which were graded by hand from the same rules), including the
PhanDA hand check of all 79 temperature bins.

    python -m unittest -v test_graded_results

Standard library only. Made with AI agents working in MV2's name.
"""
import csv
import json
import os
import unittest

import grade_temperature

HERE = os.path.dirname(os.path.abspath(__file__))


class TestResult(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(os.path.join(HERE, "results", "results.json"), encoding="utf-8") as f:
            cls.res = json.load(f)
        with open(os.path.join(HERE, "results", "boundaries.csv"), newline="",
                  encoding="utf-8") as f:
            cls.rows = list(csv.DictReader(f))
        cls.gr = grade_temperature.grade(cls.res, cls.rows)

    def test_label(self):
        self.assertEqual(self.gr["label"], "REPLICATED_ONLY_WITH_MASS_EXTINCTIONS")
        self.assertTrue(self.gr["label_agrees_with_job"])
        self.assertEqual(self.gr["status"], "GRADED")

    def test_phanda_hand_check(self):
        pc = self.gr["phanda_check"]
        self.assertEqual(pc["n_checked"], 79)
        self.assertEqual(pc["mismatches"], [])
        self.assertEqual(pc["missing_bins"], [])

    def test_rows(self):
        got = {r["id"]: r["holds"] for r in self.gr["rows"]}
        want = {"S": True, "H1": True, "H1s": True, "H2": False, "H3": False, "H4": True,
                "H5a": False, "H5b": False, "H5c": False, "H6a": False, "H6b": True,
                "H7a": False, "H7b": True, "G": True}
        self.assertEqual(got, want)

    def test_brier(self):
        self.assertEqual(self.gr["brier"]["n"], 14)
        self.assertAlmostEqual(self.gr["brier"]["brier"], 2.5684 / 14, places=6)
        bx = self.gr["brier_excluding_temperature_only"]
        self.assertEqual(bx["n"], 12)
        self.assertAlmostEqual(bx["brier"], 2.145 / 12, places=6)

    def test_markdown_renders(self):
        md = grade_temperature.render_md(self.gr, self.res)
        self.assertIn("REPLICATED_ONLY_WITH_MASS_EXTINCTIONS", md)


if __name__ == "__main__":
    unittest.main()
