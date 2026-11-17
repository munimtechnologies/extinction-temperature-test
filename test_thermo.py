# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Offline tests for thermo.py. Run: python -m unittest -v test_thermo"""
import math
import unittest

import numpy as np

import thermo


class Bins(unittest.TestCase):
    def test_structure(self):
        bins = thermo.build_bins()
        self.assertEqual(len(bins), 90)
        self.assertEqual(sum(b["phanda"] is not None for b in bins), 79)
        self.assertEqual(bins[10]["name"], "Tremadocian")
        self.assertEqual(bins[-2]["name"], "Piacenzian")
        self.assertEqual(bins[-1]["name"], "Quaternary")
        used = [s for b in bins for s in b["stages"]]
        self.assertEqual(sorted(used), sorted(s[0] for s in thermo.STAGES))
        for a, b in zip(bins[:-1], bins[1:]):
            self.assertAlmostEqual(a["top"], b["base"], places=9)
        merged = {b["name"]: b for b in bins}["Selandian/Danian"]
        self.assertEqual((merged["base"], merged["top"]), (66.0, 59.24))

    def test_slices(self):
        self.assertEqual(len(thermo.PHANDA_SLICES), 85)
        self.assertEqual(len(set(thermo.PHANDA_SLICES)), 85)
        for k in list(thermo.MERGED) + list(thermo.DROPPED_SLICES):
            self.assertIn(k, thermo.PHANDA_SLICES)

    def test_boundary_rows(self):
        bins = thermo.build_bins()
        temps = {b["phanda"]: (k - 1.0, float(k), k + 2.0)
                 for k, b in enumerate(bins) if b["phanda"]}
        z = [0.0] * len(bins)
        rows = thermo.boundary_rows(bins, z, z, temps)
        self.assertEqual(len(rows), 78)
        self.assertEqual(rows[0]["bin"], "Tremadocian")
        self.assertEqual(rows[-1]["bin"], "Zanclean")
        self.assertTrue(all(r["dT"] == 1.0 and r["abs_dT"] == 1.0 for r in rows))
        self.assertEqual(sum(r["big_five"] for r in rows), 6)
        self.assertEqual({r["bin"] for r in rows if r["big_five"]}, set(thermo.BIG_FIVE_TOPS))


class Rates(unittest.TestCase):
    def test_three_timer(self):
        S = [{"a", "b", "c"}, {"a", "b", "c"}, {"a", "b"}, {"a", "b"}, {"a", "b"}]
        mu = thermo.three_timer(S)
        self.assertAlmostEqual(mu[1], math.log(1.5), places=12)
        self.assertTrue(math.isnan(mu[0]) and math.isnan(mu[4]))

    def test_gap_filler(self):
        S = [{"a", "b", "c"}, {"a", "b", "c"}, {"a", "b"}, {"a", "b"}, {"a", "b"}]
        self.assertAlmostEqual(thermo.gap_filler(S)[1], math.log(1.5), places=12)
        # c is missing from bin 2 but seen again in bin 3: a gap, not an extinction.
        S2 = [{"a", "b", "c"}, {"a", "b", "c"}, {"a", "b"}, {"a", "b", "c"}, {"a", "b", "c"}]
        self.assertAlmostEqual(thermo.gap_filler(S2)[1], 0.0, places=12)
        self.assertTrue(math.isnan(thermo.gap_filler(S2)[3]))

    def test_merge_sets(self):
        bins = thermo.build_bins()
        stage_sets = [set() for _ in thermo.STAGES]
        names = [s[0] for s in thermo.STAGES]
        stage_sets[names.index("Danian")].add("x")
        stage_sets[names.index("Selandian")].add("y")
        counts = [1] * len(thermo.STAGES)
        sets, n = thermo.merge_sets(stage_sets, bins, counts)
        j = [b["name"] for b in bins].index("Selandian/Danian")
        self.assertEqual(sets[j], {"x", "y"})
        self.assertEqual(n[j], 2)


class Temperature(unittest.TestCase):
    def test_parse(self):
        text = ("Period,Epoch,Stage,UpperAge,LowerAge,AverageAge,GMST_05,GMST_16,"
                "GMST_50,GMST_84,GMST_95,CO2_05\n"
                "Neogene,Pliocene,Zanclean,3.6,5.33,4.465,15.0,15.3,16.3,17.6,18.0,243\n"
                "Paleogene,Paleocene,Selandian/Danian,59.24,66,62.6,1,2,3,4,5,6\n")
        t = thermo.parse_phanda(text)
        self.assertEqual(t["Zanclean"], (15.3, 16.3, 17.6))
        self.assertEqual(t["Selandian/Danian"], (2.0, 3.0, 4.0))

    def test_split_normal(self):
        v = thermo.split_normal([9.0, 9.0, 9.0], [10.0] * 3, [12.0] * 3, [-1.0, 0.0, 1.0])
        self.assertEqual(list(v), [9.0, 10.0, 12.0])


class Stats(unittest.TestCase):
    def test_rankdata(self):
        self.assertEqual(list(thermo.rankdata([10, 20, 20, 30])), [1.0, 2.5, 2.5, 4.0])
        self.assertEqual(list(thermo.rankdata([3, 1, 2])), [3.0, 1.0, 2.0])

    def test_spearman(self):
        self.assertAlmostEqual(thermo.spearman([1, 2, 3, 4], [10, 20, 30, 45]), 1.0)
        self.assertAlmostEqual(thermo.spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)

    def test_permutation_and_shift(self):
        rng = np.random.default_rng(1)
        x = np.arange(60.0)
        y = x + 0.5 * rng.standard_normal(60)
        r = thermo.spearman_test(x, y, 999, rng)
        self.assertGreater(r["rho"], 0.95)
        self.assertAlmostEqual(r["p_perm"], 0.001)
        self.assertEqual(r["n_shifts"], 55)
        self.assertTrue(0.0 < r["p_shift"] <= 1.0)
        q = thermo.spearman_test(rng.standard_normal(60), rng.standard_normal(60), 999, rng)
        self.assertTrue(0.001 <= q["p_perm"] <= 1.0)

    def test_partial(self):
        rng = np.random.default_rng(2)
        z = rng.standard_normal(200)
        x = z + 0.3 * rng.standard_normal(200)
        y = z + 0.3 * rng.standard_normal(200)
        self.assertGreater(thermo.spearman(x, y), 0.7)
        p = thermo.partial_spearman_test(x, y, z, 199, rng)
        self.assertLess(abs(p["rho"]), 0.3)

    def test_auc(self):
        self.assertEqual(thermo.auc([3, 2, 1, 0], [1, 1, 0, 0]), 1.0)
        self.assertEqual(thermo.auc([0, 1, 2, 3], [1, 1, 0, 0]), 0.0)
        self.assertEqual(thermo.auc([1, 1], [1, 0]), 0.5)

    def test_fisher(self):
        f = thermo.fisher_diff(0.5, 40, 0.5, 40)
        self.assertAlmostEqual(f["z"], 0.0)
        self.assertAlmostEqual(f["p_one_sided"], 0.5)


class Analyze(unittest.TestCase):
    def rows(self, seed, planted=True, big_mu=3.0):
        rng = np.random.default_rng(seed)
        bins = thermo.build_bins()
        temps, k = {}, 0.0
        for b in bins:
            if b["phanda"]:
                k += rng.choice([-1.0, 1.0]) * rng.uniform(0.0, 6.0)
                temps[b["phanda"]] = (k - 1.0, k, k + 1.0)
        z = [0.0] * len(bins)
        rows = thermo.boundary_rows(bins, z, z, temps)
        for r in rows:
            base = 0.3 * r["abs_dT"] if planted else rng.uniform(0.0, 2.0)
            r["mu3t"] = base + 0.05 * rng.standard_normal()
            r["mugf"] = r["mu3t"]
            if r["big_five"]:
                r["mu3t"] = r["mugf"] = big_mu
        return rows

    def test_planted(self):
        res = thermo.analyze(self.rows(3), 500, np.random.default_rng(4), n_draws=50)
        self.assertEqual(res["n"], 78)
        self.assertTrue(res["sanity"]["passed"])
        self.assertEqual(res["sanity"]["big_five_events_in_top10"], 5)
        self.assertGreater(res["H2_without_big_five"]["rho"], 0.9)
        self.assertEqual(res["H2_without_big_five"]["n"], 72)
        self.assertEqual(res["decision"], "REPLICATED")
        self.assertGreater(res["gap_filler"]["H2_without_big_five"]["rho"], 0.9)
        self.assertTrue(res["qualifiers"]["duration_robust_H1"])

    def test_null_runs(self):
        res = thermo.analyze(self.rows(5, planted=False), 300, np.random.default_rng(6),
                             n_draws=50)
        for key in ("H1_all", "H2_without_big_five", "H7_level"):
            self.assertTrue(-1.0 <= res[key]["rho"] <= 1.0)
            self.assertTrue(0.0 < res[key]["p_perm"] <= 1.0)
        self.assertTrue(0.0 <= res["H4_uncertainty"]["frac_rho_pos"] <= 1.0)
        self.assertIn(res["decision"], ("REPLICATED", "REPLICATED_ONLY_WITH_MASS_EXTINCTIONS",
                                        "NOT_REPLICATED", "DATA_LIMITED"))

    def test_sanity_gate(self):
        res = thermo.analyze(self.rows(3, big_mu=-5.0), 200, np.random.default_rng(1),
                             n_draws=20)
        self.assertFalse(res["sanity"]["passed"])
        self.assertEqual(res["decision"], "DATA_LIMITED")

    def test_too_few(self):
        rows = self.rows(3)
        for r in rows:
            r["mu3t"] = float("nan")
        res = thermo.analyze(rows, 100, np.random.default_rng(1))
        self.assertEqual(res["decision"], "DATA_LIMITED")
        self.assertEqual(res["n"], 0)


if __name__ == "__main__":
    unittest.main()
