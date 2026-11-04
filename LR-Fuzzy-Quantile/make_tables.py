"""Build the LaTeX tables of the article (Sections 5.4 and 5.6) from the result files, so the paper and the code agree.

Reads  horizontal_results.json  and  inference_results.json  (made by lr_fuzzy_horizontal.py / lr_fuzzy_inference.py,
FULL runs, seed 2026) and writes  t_hor.tex, t_recover.tex, t_gof.tex, t_mle.tex  into --outdir.
Usage:  python make_tables.py --results results --outdir tables
"""
import argparse
import json
import os


def _w(path, header, rows, colspec):
    with open(path, "w") as f:
        f.write("\\begin{tabular}{%s}\n\\toprule\n%s \\\\\n\\midrule\n" % (colspec, header))
        for r in rows:
            f.write(" & ".join(r) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--outdir", default="tables")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    H = json.load(open(os.path.join(a.results, "horizontal_results.json")))
    I = json.load(open(os.path.join(a.results, "inference_results.json")))
    if H.get("reps", 0) < 1000:
        print("WARNING: horizontal results come from a quick run; use the full run for the article tables.")

    rows = [[r["scenario"], str(r["k"]), "%.3f" % r["median"], "%.3f" % r["q95"], "%.3f" % r["median_bound"],
             "%.3f" % r["coverage_of_bound"]]
            for r in H["horizontal"] if r["k"] in (16, 64, 256, 1000)]
    _w(os.path.join(a.outdir, "t_hor.tex"), "Sc. & $k$ & Median & 95\\% q. & Bound & Cov.", rows, "llrrrr")

    rows = [[str(r["n"]), "%.4f" % r["err_fuzzy"], "%.3f" % r["n_times_err_fuzzy"], "%.4f" % r["exact_expectation"],
             "%.4f" % r["err_crisp_edf"], "%.3f" % r["sqrt_n_times_err_crisp"], "%.1f" % r["ratio"]] for r in I["recover"]]
    _w(os.path.join(a.outdir, "t_recover.tex"),
       "$n$ & Fuzzy err. & $n\\times$err & Exact E & EDF err. & $\\sqrt n\\times$err & Ratio", rows, "rrrrrrr")

    rows = [[str(r["n"]), "%.3f" % r["H0_rejection"], "%.3f" % r["triangular_rejection"],
             "%.3f" % r["misspecified_rejection"]] for r in I["gof"]]
    _w(os.path.join(a.outdir, "t_gof.tex"), "$n$ & Size & Power (tri.) & Power (mism.)", rows, "rrrr")

    rows = [[str(r["n"]), str(r["k"]), "%.3f" % r["efficiency_gain"], "%.3f" % r["sd"], "%.3f" % r["sd_theory"],
             "%.3f" % r["mean"]] for r in I["mle"]]
    _w(os.path.join(a.outdir, "t_mle.tex"), "$n$ & $k$ & $I_k/I_s$ & Sim. s.d. & Theor. s.d. & Mean", rows, "rrrrrr")

    # data files for the pgfplots figures
    fd = os.path.join(a.outdir, "figdata"); os.makedirs(fd, exist_ok=True)
    with open(os.path.join(fd, "gof.dat"), "w") as f:
        f.write("n H0 tri mis\n")
        for r in I["gof"]:
            f.write("%d %g %g %g\n" % (r["n"], r["H0_rejection"], r["triangular_rejection"], r["misspecified_rejection"]))
    with open(os.path.join(fd, "mle.dat"), "w") as f:
        f.write("n k sd sdt\n")
        for r in I["mle"]:
            f.write("%d %d %g %g\n" % (r["n"], r["k"], r["sd"], r["sd_theory"]))
    for n in sorted({r["n"] for r in I["mle"]}):
        with open(os.path.join(fd, "mle%d.dat" % n), "w") as f:
            f.write("k sd sdt\n")
            for r in I["mle"]:
                if r["n"] == n:
                    f.write("%d %g %g\n" % (r["k"], r["sd"], r["sd_theory"]))
    print("wrote figure data in", fd)


if __name__ == "__main__":
    main()
