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
       (3b) make_tables.py writes the LaTeX tables of the article (results/latex_tables);
       (4) merge everything into ONE Excel workbook, LR_all_results.xlsx (one sheet per table);
       (5) save the workbook, the individual files and the figures to Google Drive.

Google Drive
------------
In Google Colab the script mounts Drive and writes to  My Drive/LR_Fuzzy_Quantile_results/  (override with
--drive-folder NAME). Outside Colab use --drive-dir PATH (e.g. a synced Google Drive folder on your computer);
without it the results stay in <workdir>/results. --no-drive disables saving to Drive.
"""
import argparse
import io
import json
import os
import shutil
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


# ---------------------------------------------------------------- one Excel workbook with all results
def _frames_from_json(path, prefix):
    """Turn a result JSON file into {sheet_name: DataFrame} (lists of dicts -> tables, dicts -> Parameter/Value)."""
    import pandas as pd
    data = json.load(open(path))
    frames = {}
    for key, val in data.items():
        name = ("%s_%s" % (prefix, key))[:31]
        if isinstance(val, list) and val and isinstance(val[0], dict):
            frames[name] = pd.DataFrame(val)
        elif isinstance(val, dict):
            frames[name] = pd.DataFrame({"Parameter": list(val.keys()), "Value": list(val.values())})
        else:
            frames[name] = pd.DataFrame({"Parameter": [key], "Value": [val]})
    return frames


def build_workbook(res_dir, out_path):
    """Merge simulation workbook + inference JSON + horizontal JSON into one .xlsx (styled header, frozen row)."""
    import pandas as pd
    from openpyxl import load_workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    frames = {}
    sim = os.path.join(res_dir, "LR_fuzzy_simulation_results.xlsx")
    if os.path.exists(sim):
        for name, df in pd.read_excel(sim, sheet_name=None).items():
            frames["sim_" + name[:27]] = df
    for fname, prefix in (("inference_results.json", "inf"), ("horizontal_results.json", "hor")):
        path = os.path.join(res_dir, fname)
        if os.path.exists(path):
            frames.update(_frames_from_json(path, prefix))
    if not frames:
        print("no result tables found; workbook not written")
        return None
    with pd.ExcelWriter(out_path, engine="openpyxl") as xw:
        for name, df in frames.items():
            df.to_excel(xw, sheet_name=name[:31], index=False)
    wb = load_workbook(out_path)
    head = PatternFill("solid", fgColor="1F4E78")
    for ws in wb.worksheets:
        for c in ws[1]:
            c.font = Font(bold=True, color="FFFFFF"); c.fill = head
            c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        for col in ws.columns:
            w = max(len(str(c.value)) if c.value is not None else 0 for c in col[:60])
            ws.column_dimensions[col[0].column_letter].width = min(max(12, w + 2), 45)
        ws.freeze_panes = "A2"
    wb.save(out_path)
    print("workbook with %d sheets -> %s" % (len(frames), out_path))
    return out_path


# ---------------------------------------------------------------- Google Drive
def in_colab():
    try:
        import google.colab  # noqa: F401
        return True
    except ImportError:
        return False


def drive_target(args):
    """Return the destination folder on Drive (created), or None if saving to Drive is not possible/wanted."""
    if args.no_drive:
        return None
    if args.drive_dir:
        path = os.path.abspath(os.path.expanduser(args.drive_dir))
    elif in_colab():
        root = "/content/drive/MyDrive"
        if not os.path.isdir(root):
            try:
                from google.colab import drive
                drive.mount("/content/drive")
            except Exception:
                pass
        if not os.path.isdir(root):
            print("\nGoogle Drive is not mounted (mounting only works from a notebook cell, not from '!python').\n"
                  "Run  from google.colab import drive; drive.mount('/content/drive')  in a cell first, then run again.\n"
                  "Results are kept in the work folder meanwhile.")
            return None
        path = os.path.join(root, args.drive_folder)
    else:
        print("\nNot running in Colab and --drive-dir not given: results stay in the work folder "
              "(use --drive-dir PATH for a synced Google Drive folder).")
        return None
    os.makedirs(path, exist_ok=True)
    return path


def copy_to_drive(res_dir, workbook, dest):
    """Copy the combined workbook (top level) and the whole results folder (sub-folder 'details')."""
    if workbook:
        shutil.copy2(workbook, os.path.join(dest, os.path.basename(workbook)))
    det = os.path.join(dest, "details")
    if os.path.isdir(det):
        shutil.rmtree(det)
    shutil.copytree(res_dir, det)
    print("\nsaved to Drive folder:", dest)
    for root, _, files in os.walk(dest):
        for f in sorted(files):
            print("  ", os.path.relpath(os.path.join(root, f), dest))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", default="lr_fuzzy_run", help="where to put the code and the results")
    ap.add_argument("--quick", action="store_true", help="small replication numbers (fast test)")
    ap.add_argument("--only", nargs="+", choices=["simulation", "inference", "horizontal", "base_paper"],
                    help="run only these parts (default: all)")
    ap.add_argument("--no-install", action="store_true", help="do not pip install the packages")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--drive-folder", default="LR_Fuzzy_Quantile_results",
                    help="Colab: folder name inside My Drive (default LR_Fuzzy_Quantile_results)")
    ap.add_argument("--drive-dir", default=None,
                    help="explicit destination folder (e.g. a synced Google Drive folder outside Colab)")
    ap.add_argument("--no-drive", action="store_true", help="do not save to Google Drive")
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

    if os.path.exists(os.path.join(res, "horizontal_results.json")) and os.path.exists(os.path.join(res, "inference_results.json")):
        status["tables"] = run([py, "make_tables.py", "--results", res, "--outdir", os.path.join(res, "latex_tables")], code)
    workbook = build_workbook(res, os.path.join(res, "LR_all_results.xlsx"))
    dest = drive_target(a)
    if dest:
        copy_to_drive(res, workbook, dest)

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
