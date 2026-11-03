#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download the LR-Fuzzy-Quantile code from GitHub and run all experiments.

Works in a plain Python session and in Google Colab (paste the file into a cell, or run
    !wget -q https://raw.githubusercontent.com/yamenetoo/Statistical-Repository/main/LR-Fuzzy-Quantile/get_and_run.py
    !python get_and_run.py --quick
).

Usage
-----
    python get_and_run.py                 # download + install packages + run everything (full run, ~10 min)
    python get_and_run.py --quick         # same with small replication numbers (about a minute)
    python get_and_run.py --only horizontal inference
    python get_and_run.py --workdir my_run --no-install

Steps: (1) download the repository archive (no git needed), keep the folder LR-Fuzzy-Quantile;
       (2) pip install numpy scipy pandas openpyxl matplotlib (skip with --no-install);
       (3) run: simulation (Excel), inference (JSON), horizontal (JSON), base_paper (figures);
       (4) print the list of output files.
"""
import argparse
import io
import os
import subprocess
import sys
import time
import urllib.request
import zipfile

OWNER, REPO, BRANCH, FOLDER = "yamenetoo", "Statistical-Repository", "main", "LR-Fuzzy-Quantile"
ZIP_URL = "https://github.com/%s/%s/archive/refs/heads/%s.zip" % (OWNER, REPO, BRANCH)
PACKAGES = ["numpy", "scipy", "pandas", "openpyxl", "matplotlib"]


def download(workdir):
    """Download the repository zip and extract only the LR-Fuzzy-Quantile folder."""
    target = os.path.join(workdir, FOLDER)
    if os.path.isdir(target):
        print("folder already present:", target, "(delete it to download again)")
        return target
    print("downloading", ZIP_URL)
    req = urllib.request.Request(ZIP_URL, headers={"User-Agent": "lr-fuzzy-quantile-runner"})
    data = urllib.request.urlopen(req, timeout=120).read()
    prefix = "%s-%s/%s/" % (REPO, BRANCH, FOLDER)
    n = 0
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for name in z.namelist():
            if name.startswith(prefix) and not name.endswith("/"):
                out = os.path.join(target, name[len(prefix):])
                os.makedirs(os.path.dirname(out), exist_ok=True)
                with z.open(name) as src, open(out, "wb") as dst:
                    dst.write(src.read())
                n += 1
    if n == 0:
        sys.exit("folder %s not found in the archive (is the repository public and the branch name right?)" % FOLDER)
    print("extracted %d files to %s" % (n, target))
    return target


def pip_install():
    print("installing packages:", " ".join(PACKAGES))
    cmd = [sys.executable, "-m", "pip", "install", "-q"] + PACKAGES
    r = subprocess.run(cmd)
    if r.returncode != 0:     # externally-managed environments
        subprocess.run(cmd + ["--break-system-packages"], check=True)


def run(cmd, cwd):
    print("\n$ " + " ".join(cmd))
    t0 = time.time()
    r = subprocess.run(cmd, cwd=cwd)
    print("-> exit code %d in %.0f s" % (r.returncode, time.time() - t0))
    return r.returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", default="lr_fuzzy_run", help="where to put the code and the results")
    ap.add_argument("--quick", action="store_true", help="small replication numbers (fast test)")
    ap.add_argument("--only", nargs="+", choices=["simulation", "inference", "horizontal", "base_paper"],
                    help="run only these parts (default: all)")
    ap.add_argument("--no-install", action="store_true", help="do not pip install the packages")
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()

    os.makedirs(a.workdir, exist_ok=True)
    code = download(os.path.abspath(a.workdir))
    if not a.no_install:
        pip_install()
    res = os.path.join(os.path.abspath(a.workdir), "results")
    os.makedirs(res, exist_ok=True)
    py = sys.executable
    q = ["--quick"] if a.quick else []
    todo = a.only or ["simulation", "inference", "horizontal", "base_paper"]
    status = {}
    if "simulation" in todo:
        status["simulation"] = run([py, "lr_fuzzy_simulation.py", "--seed", str(a.seed), "--out",
                                    os.path.join(res, "LR_fuzzy_simulation_results.xlsx")] + q, code)
    if "inference" in todo:
        status["inference"] = run([py, "lr_fuzzy_inference.py", "--seed", str(a.seed), "--out",
                                   os.path.join(res, "inference_results.json")] + q, code)
    if "horizontal" in todo:
        status["horizontal"] = run([py, "lr_fuzzy_horizontal.py", "--seed", str(a.seed), "--out",
                                    os.path.join(res, "horizontal_results.json"),
                                    "--datadir", os.path.join(res, "figdata")] + q, code)
    if "base_paper" in todo:
        status["base_paper"] = run([py, "parchami2024_base_paper.py", "--example", "all", "--outdir",
                                    os.path.join(res, "base_paper_figures")],
                                   os.path.join(code, "base_paper"))

    print("\n================ summary ================")
    for k, v in status.items():
        print("%-12s %s" % (k, "OK" if v == 0 else "FAILED (exit %d)" % v))
    print("\nresult files in", res)
    for root, _, files in os.walk(res):
        for f in sorted(files):
            p = os.path.join(root, f)
            print("  %-70s %8.1f KB" % (os.path.relpath(p, res), os.path.getsize(p) / 1024))
    if any(v != 0 for v in status.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
