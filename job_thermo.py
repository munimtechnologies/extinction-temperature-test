# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Run the pre-registered test in PREREG.md end to end.

    python job_thermo.py --data-dir /data --out-dir /out

Downloads PBDB marine invertebrate genus occurrences (same queries and cache
file names as extinction-cycle-test/job.py) and the PhanDA GMST percentiles
(Judd et al. 2024, CC BY 4.0), computes boundary extinction rates on PhanDA's
time slices and runs the registered tests.
Made with AI agents working in MV2's name.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.request

import numpy as np

import thermo

PBDB = "https://paleobiodb.org/data1.2/occs/list.csv"
GROUPS = ["Porifera", "Cnidaria", "Bryozoa", "Brachiopoda", "Mollusca",
          "Echinodermata", "Trilobita", "Graptolithina"]
REQUIRED = {"Brachiopoda", "Mollusca", "Echinodermata", "Trilobita"}
COMMON = "envtype=marine&idreso=genus&idqual=genus_certain&pres=regular&show=class"
PHANDA_URL = ("https://raw.githubusercontent.com/EJJudd/PhanDA/main/5_Outputs/"
              "PhanDA_GMSTandCO2_percentiles.csv")
UA = "MV2-extinction-temperature/0.1 (Munim, Inc.; mv2@munimtech.com)"


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, dest, timeout=1800, tries=3):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return "cached"
    tmp = f"{dest}.part{os.getpid()}"
    last = ""
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            os.replace(tmp, dest)
            return "downloaded"
        except Exception as e:  # noqa: BLE001 - report and retry
            last = repr(e)
            print(f"fetch failed ({k + 1}/{tries}) {url}: {last}", flush=True)
            time.sleep(30 * (k + 1))
    raise RuntimeError(last)


def describe(path, url, how):
    return {"url": url, "status": how, "bytes": os.path.getsize(path), "md5": md5(path),
            "file_mtime_utc": dt.datetime.fromtimestamp(
                os.path.getmtime(path), dt.timezone.utc).isoformat()}


def get_pbdb(data_dir, prov):
    d = os.path.join(data_dir, "pbdb")
    os.makedirs(d, exist_ok=True)
    paths = []
    for g in GROUPS:
        url = f"{PBDB}?base_name={g}&{COMMON}"
        dest = os.path.join(d, f"{g}.csv")
        try:
            how = fetch(url, dest)
            with open(dest, encoding="utf-8", errors="replace") as f:
                head = f.readline()
            if "occurrence_no" not in head:
                os.remove(dest)
                raise RuntimeError("unexpected header: " + head[:300])
            prov[g] = describe(dest, url, how)
            paths.append((g, dest))
        except Exception as e:  # noqa: BLE001
            prov[g] = {"url": url, "error": repr(e)[:500]}
        print(g, prov[g].get("status", prov[g].get("error")), flush=True)
    return paths


def get_phanda(data_dir, prov):
    d = os.path.join(data_dir, "phanda")
    os.makedirs(d, exist_ok=True)
    dest = os.path.join(d, "PhanDA_GMSTandCO2_percentiles.csv")
    how = fetch(PHANDA_URL, dest)
    prov.update(describe(dest, PHANDA_URL, how))
    with open(dest, encoding="utf-8", errors="replace") as f:
        return thermo.parse_phanda(f.read())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20261117)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    prov = {"started_utc": now(), "pbdb": {}, "phanda": {}, "groups": GROUPS,
            "query_common": COMMON, "n_perm": args.n_perm, "seed": args.seed}
    res = {"prereg": "PREREG.md (extinction-temperature-test v0.1.0)", "status": "ok"}

    try:
        temps = get_phanda(args.data_dir, prov["phanda"])
    except Exception as e:  # noqa: BLE001
        prov["phanda"]["error"] = repr(e)[:500]
        temps = {}
    got_slices = set(temps)
    prov["phanda"]["n_slices"] = len(got_slices)
    prov["phanda"]["missing_expected"] = sorted(set(thermo.PHANDA_SLICES) - got_slices)
    prov["phanda"]["unexpected"] = sorted(got_slices - set(thermo.PHANDA_SLICES))
    bins = thermo.build_bins()
    needed = {b["phanda"] for b in bins if b["phanda"]}
    if not needed <= got_slices:
        res["status"] = "phanda_incomplete: missing " + ",".join(sorted(needed - got_slices))
        return _write(args.out_dir, res, prov, 2)

    paths = get_pbdb(args.data_dir, prov["pbdb"])
    got = {g for g, _ in paths}
    if not REQUIRED <= got:
        res["status"] = "pbdb_incomplete: missing " + ",".join(sorted(REQUIRED - got))
        return _write(args.out_dir, res, prov, 2)

    assigner = thermo.BinAssigner()
    acc = thermo.Accumulator(assigner.n)
    for _, p in paths:
        with open(p, newline="", encoding="utf-8", errors="replace") as f:
            acc.add_rows(csv.DictReader(f), assigner)
    prov["parse"] = acc.stats

    sets, n_occ = thermo.merge_sets(acc.sets, bins, acc.n_occ)
    mu3t = thermo.three_timer(sets)
    mugf = thermo.gap_filler(sets)
    rows = thermo.boundary_rows(bins, mu3t, mugf, temps)

    with open(os.path.join(args.out_dir, "bins.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["bin", "stages", "phanda_slice", "base_ma", "top_ma", "n_occ",
                    "n_genera", "mu3t_top", "mugf_top"])
        for j, b in enumerate(bins):
            w.writerow([b["name"], "+".join(b["stages"]), b["phanda"] or "", b["base"],
                        b["top"], n_occ[j], len(sets[j]), mu3t[j], mugf[j]])
    keys = ["bin", "top_ma", "mu3t", "mugf", "T16_bin", "T_bin", "T84_bin", "T16_next",
            "T_next", "T84_next", "dT", "abs_dT", "dur_bin", "dur_next", "big_five"]
    with open(os.path.join(args.out_dir, "boundaries.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in rows:
            w.writerow([r[k] for k in keys])

    res["analysis"] = thermo.analyze(rows, args.n_perm, rng)
    res["decision"] = res["analysis"]["decision"]
    prov["finished_utc"] = now()
    return _write(args.out_dir, res, prov, 0)


def _write(out_dir, res, prov, code):
    def conv(o):
        return o.item() if hasattr(o, "item") else str(o)
    with open(os.path.join(out_dir, "results.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, default=conv)
    with open(os.path.join(out_dir, "provenance.json"), "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=1, default=str)
    a = res.get("analysis", {})
    print(json.dumps({"status": res.get("status"), "decision": res.get("decision"),
                      "sanity": a.get("sanity"), "H1": a.get("H1_all"),
                      "H2": a.get("H2_without_big_five")}, default=conv), flush=True)
    return code


if __name__ == "__main__":
    sys.exit(main())
