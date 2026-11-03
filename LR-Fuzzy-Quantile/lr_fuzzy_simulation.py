#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulation study for the proposed contribution to
  Parchami, Grzegorzewski, Romaniuk (2024), "Statistical simulations with LR random fuzzy numbers".

Contribution (Quantile Representation):
    X~_alpha = [ C - Q_l((1-alpha) V_l) ,  C + Q_r((1-alpha) V_r) ],   V_l, V_r ~ U(0,1)
    (where V = F_S(S)), plus optional kernel width W (fuzzy interval) and
    a Gaussian copula linking (C, V_l, V_r).

Experiments (each one is written to a sheet of the output Excel file):
  1. Equivalence   : paper's Algorithm 2 (truncated cdf)  vs  quantile representation (Alg. 2')
  2. Convergence   : sup-distance between k-knot polygon (paper's Algorithm 1) and its LR limit,
                     compared with the conjectured DKW-type bound  sqrt(ln(4/delta)/(2k)) + 2/(k+1)
  3. Aumann        : closed-form / numerical expected alpha-cuts vs Monte-Carlo
  4. Dependence    : Gaussian copula between left and right fuzziness (rho_LR) -> effect on width
  5. Timing        : Algorithm 1 vs Algorithm 2 vs Algorithm 2' (alpha-cuts on a grid)
  6. Fuzzy_sample  : a simulated fuzzy random sample (alpha-cut endpoints) for plotting / later use
  7. Fuzzy_interval: kernel of positive width W (trapezoid-like numbers)

Requirements:  pip install numpy scipy pandas openpyxl
Run:           python lr_fuzzy_simulation.py                (default settings)
               python lr_fuzzy_simulation.py --reps 2000 --seed 7 --out my_results.xlsx
               python lr_fuzzy_simulation.py --quick        (fast test run)
"""
import argparse
import time
import numpy as np
import pandas as pd
from scipy import stats, integrate
from openpyxl import load_workbook
from openpyxl.chart import LineChart, Reference
from openpyxl.styles import Font, PatternFill, Alignment

# ----------------------------------------------------------------------------
# 0. SCENARIOS  (edit here to try other distributions)
#    C: distribution of the core (any scipy frozen dist)
#    L / R: distributions of S^l, S^r  (must be continuous on (0, inf) or (0, xi))
# ----------------------------------------------------------------------------
SCENARIOS = {
    "A: N(0,1), Exp(3), Beta(3,1)": dict(
        C=stats.norm(0, 1), L=stats.expon(scale=1 / 3), R=stats.beta(3, 1)),
    "B: N(0,1), Chi2(8), Chi2(8)": dict(
        C=stats.norm(0, 1), L=stats.chi2(8), R=stats.chi2(8)),
    "C: N(0,1), Gamma(2,1), U(0,3)": dict(
        C=stats.norm(0, 1), L=stats.gamma(2, scale=1), R=stats.uniform(0, 3)),
}
ALPHAS = np.round(np.linspace(0, 1, 11), 2)          # alpha grid for tables


# ----------------------------------------------------------------------------
# 1. GENERATORS
# ----------------------------------------------------------------------------
def trunc_sample(dist, s, k, rng):
    """k iid draws from `dist` truncated to [0, s] (inverse-cdf method)."""
    u = rng.random(k) * dist.cdf(s)
    return dist.ppf(u)


def alg1_arm(dist, s, k, rng):
    """Paper's Algorithm 1 for one arm: returns (offsets, memberships) of the polygon.
    offset t = distance from the core (C - t on the left arm, C + t on the right arm)."""
    l = np.sort(trunc_sample(dist, s, k, rng))
    t = np.concatenate(([0.0], l, [s]))
    mu = np.concatenate(([1.0], (k - np.arange(1, k + 1) + 1) / (k + 1), [0.0]))
    return t, mu


def alg2_membership(dist, s, t):
    """Paper's Algorithm 2 / Corollary 1 : mu(t) = 1 - F_{S|s}(t)  for 0<=t<=s."""
    t = np.clip(t, 0, s)
    return 1.0 - dist.cdf(t) / dist.cdf(s)


def quantile_offset(dist, v, alpha):
    """Proposed Algorithm 2': alpha-cut offset  Q((1-alpha) * V),  V = F(S)."""
    return dist.ppf((1.0 - alpha) * v)


def gaussian_copula_uniforms(n, corr, rng):
    z = rng.multivariate_normal(np.zeros(corr.shape[0]), corr, size=n)
    return stats.norm.cdf(z)


# ----------------------------------------------------------------------------
# 2. EXPERIMENTS
# ----------------------------------------------------------------------------
def exp_equivalence(reps, rng):
    """Alg.2 membership at the endpoint produced by the quantile formula must equal alpha."""
    rows = []
    for name, sc in SCENARIOS.items():
        for arm in ("L", "R"):
            d = sc[arm]
            s = d.rvs(reps, random_state=rng)
            v = d.cdf(s)
            max_err = 0.0
            for a in ALPHAS[1:-1]:
                t = quantile_offset(d, v, a)
                mu = np.array([alg2_membership(d, si, ti) for si, ti in zip(s, t)])
                max_err = max(max_err, float(np.max(np.abs(mu - a))))
            rows.append(dict(Scenario=name, Arm=arm, Replications=reps,
                             Max_abs_error_Alg2_vs_quantile=max_err))
    return pd.DataFrame(rows)


def exp_convergence(ks, reps, delta, rng, grid=1500):
    rows = []
    for name, sc in SCENARIOS.items():
        for k in ks:
            errs = np.empty(reps)
            for r in range(reps):
                e = 0.0
                for arm in ("L", "R"):
                    d = sc[arm]
                    s = d.rvs(random_state=rng)
                    t_pts, mu_pts = alg1_arm(d, s, k, rng)
                    tg = np.linspace(0, s, grid)
                    approx = np.interp(tg, t_pts, mu_pts)
                    exact = alg2_membership(d, s, tg)
                    e = max(e, float(np.max(np.abs(approx - exact))))
                errs[r] = e
            bound = np.sqrt(np.log(4 / delta) / (2 * k)) + 2 / (k + 1)
            rows.append(dict(Scenario=name, k=k, Mean_sup_error=errs.mean(),
                             Median=np.median(errs), Q95=np.quantile(errs, 0.95),
                             Max=errs.max(), Bound_conjectured=bound,
                             Coverage_of_bound=float(np.mean(errs <= bound)),
                             Mean_err_times_sqrt_k=errs.mean() * np.sqrt(k)))
    return pd.DataFrame(rows)


def exp_aumann(reps, rng):
    """Expected left/right offset of the alpha-cut:  E = (1/(1-a)) * int_0^{1-a} Q(u) du."""
    rows = []
    for name, sc in SCENARIOS.items():
        for arm in ("L", "R"):
            d = sc[arm]
            v = d.cdf(d.rvs(reps, random_state=rng))        # V = F(S) ~ U(0,1)
            for a in ALPHAS[:-1]:
                theo, _ = integrate.quad(lambda u: d.ppf(u), 0, 1 - a, limit=200)
                theo /= (1 - a)
                mc = quantile_offset(d, v, a)
                # exact closed forms for two special cases
                closed = np.nan
                if d.dist.name == "expon":
                    lam = 1 / d.kwds.get("scale", 1.0)
                    closed = (1 / lam) * (1 + (a * np.log(a) / (1 - a) if a > 0 else 0.0))
                elif d.dist.name == "uniform" and d.args == (0, 3) or \
                        (d.dist.name == "uniform" and d.kwds.get("scale", None) is not None):
                    closed = None
                if d.dist.name == "uniform":
                    width = d.args[1] if len(d.args) > 1 else 1.0
                    closed = width * (1 - a) / 2
                rows.append(dict(Scenario=name, Arm=arm, alpha=a,
                                 Theory_numeric=theo, Closed_form=closed,
                                 MonteCarlo_mean=mc.mean(),
                                 MC_std_error=mc.std(ddof=1) / np.sqrt(reps),
                                 Abs_diff=abs(theo - mc.mean())))
    return pd.DataFrame(rows)


def exp_dependence(n, rng, rhos=(-0.9, -0.5, 0.0, 0.5, 0.9)):
    """Gaussian copula on (C, V_l, V_r): correlation rho_LR between V_l and V_r
    (and rho_CS between C and both arms) -> effect on width and on the arms."""
    name, sc = list(SCENARIOS.items())[0]
    rows = []
    for rho_cs in (0.0, 0.5):
        for rho in rhos:
            R = np.array([[1, rho_cs, rho_cs], [rho_cs, 1, rho], [rho_cs, rho, 1]])
            if np.linalg.eigvalsh(R).min() < 1e-9:      # not a valid correlation matrix
                continue
            u = gaussian_copula_uniforms(n, R, rng)
            C = sc["C"].ppf(u[:, 0])
            sl = sc["L"].ppf(u[:, 1])
            sr = sc["R"].ppf(u[:, 2])
            width = sl + sr
            rows.append(dict(Scenario=name, rho_C_with_arms=rho_cs, rho_L_R=rho,
                             Mean_width=width.mean(), Var_width=width.var(ddof=1),
                             Pearson_Sl_Sr=np.corrcoef(sl, sr)[0, 1],
                             Spearman_Sl_Sr=stats.spearmanr(sl, sr)[0],
                             Pearson_C_width=np.corrcoef(C, width)[0, 1]))
    return pd.DataFrame(rows)


def exp_timing(nfuzzy, k, rng):
    name, sc = list(SCENARIOS.items())[0]
    al = np.linspace(0, 1, 101)
    dl, dr = sc["L"], sc["R"]
    out = []

    t0 = time.perf_counter()                                     # Alg 1 (k-knot)
    for _ in range(nfuzzy):
        for d in (dl, dr):
            s = d.rvs(random_state=rng)
            t_pts, mu_pts = alg1_arm(d, s, k, rng)
            np.interp(al, mu_pts[::-1], t_pts[::-1])
    out.append(("Algorithm 1 (k-knot polygon, k=%d)" % k, time.perf_counter() - t0))

    t0 = time.perf_counter()                                     # Alg 2 (truncated cdf + inversion on grid)
    for _ in range(nfuzzy):
        for d in (dl, dr):
            s = d.rvs(random_state=rng)
            tg = np.linspace(0, s, 400)
            mu = alg2_membership(d, s, tg)
            np.interp(al, mu[::-1], tg[::-1])
    out.append(("Algorithm 2 (truncated cdf, grid inversion)", time.perf_counter() - t0))

    t0 = time.perf_counter()                                     # Alg 2' (quantile)
    for _ in range(nfuzzy):
        for d in (dl, dr):
            v = rng.random()
            quantile_offset(d, v, al)
    out.append(("Algorithm 2' (quantile representation)", time.perf_counter() - t0))

    df = pd.DataFrame(out, columns=["Method", "Seconds_total"])
    df["Fuzzy_numbers"] = nfuzzy
    df["Microsec_per_fuzzy_number"] = df["Seconds_total"] / nfuzzy * 1e6
    df["Speedup_vs_Alg1"] = df["Seconds_total"].iloc[0] / df["Seconds_total"]
    return df


def exp_sample(n, rng, width=None):
    """Fuzzy random sample by the quantile representation; optional kernel width W."""
    name, sc = list(SCENARIOS.items())[0]
    C = sc["C"].rvs(n, random_state=rng)
    vl, vr = rng.random(n), rng.random(n)
    W = np.zeros(n) if width is None else width.rvs(n, random_state=rng)
    rows = []
    for a in ALPHAS:
        left = C - quantile_offset(sc["L"], vl, a)
        right = C + W + quantile_offset(sc["R"], vr, a)
        for i in range(n):
            rows.append(dict(Fuzzy_number=i + 1, alpha=a, Left=left[i], Right=right[i]))
    df = pd.DataFrame(rows)
    wide = df.pivot(index="alpha", columns="Fuzzy_number", values=["Left", "Right"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    return df, wide.reset_index()


# ----------------------------------------------------------------------------
# 3. EXCEL OUTPUT
# ----------------------------------------------------------------------------
def style_sheet(ws):
    head = PatternFill("solid", fgColor="1F4E78")
    for c in ws[1]:
        c.font = Font(bold=True, color="FFFFFF")
        c.fill = head
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    for col in ws.columns:
        width = max(len(str(c.value)) if c.value is not None else 0 for c in col[:60])
        ws.column_dimensions[col[0].column_letter].width = min(max(12, width + 2), 45)
    ws.freeze_panes = "A2"


def add_convergence_chart(wb, df):
    ws = wb["Convergence"]
    first = df["Scenario"].iloc[0]
    n_k = int((df["Scenario"] == first).sum())
    ch = LineChart()
    ch.title = f"Sup error vs k  ({first})"
    ch.y_axis.title = "sup |X_k - X|"
    ch.x_axis.title = "k (row order)"
    cols = [list(df.columns).index(c) + 1 for c in ("Mean_sup_error", "Q95", "Bound_conjectured")]
    for ci in cols:
        ch.add_data(Reference(ws, min_col=ci, min_row=1, max_row=1 + n_k), titles_from_data=True)
    ch.set_categories(Reference(ws, min_col=list(df.columns).index("k") + 1,
                                min_row=2, max_row=1 + n_k))
    ch.height, ch.width = 9, 18
    ws.add_chart(ch, "L2")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=500, help="replications per cell")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", default="LR_fuzzy_simulation_results.xlsx")
    ap.add_argument("--quick", action="store_true", help="small run for testing")
    a = ap.parse_args()

    reps = 40 if a.quick else a.reps
    ks = [2, 4, 8, 16, 50, 200] if a.quick else [2, 4, 8, 16, 50, 200, 1000]
    delta = 0.05
    rng = np.random.default_rng(a.seed)

    print("1/7 equivalence ...");  df_eq = exp_equivalence(max(reps, 100), rng)
    print("2/7 convergence ...");  df_cv = exp_convergence(ks, reps, delta, rng)
    print("3/7 Aumann expectation ..."); df_au = exp_aumann(max(reps * 20, 5000), rng)
    print("4/7 dependence ...");   df_dp = exp_dependence(max(reps * 20, 5000), rng)
    print("5/7 timing ...");       df_tm = exp_timing(100 if a.quick else 1000, 200, rng)
    print("6/7 fuzzy sample ...")
    df_sl, df_sw = exp_sample(10, rng)
    print("7/7 fuzzy interval ...")
    df_fi, df_fw = exp_sample(10, rng, width=stats.expon(scale=0.5))

    cfg = pd.DataFrame({
        "Parameter": ["seed", "replications", "k values", "delta (bound)", "alpha grid",
                      "Scenarios", "Note"],
        "Value": [a.seed, reps, str(ks), delta, str(list(ALPHAS)),
                  " | ".join(SCENARIOS.keys()),
                  "Bound = sqrt(ln(4/delta)/(2k)) + 2/(k+1) is a conjecture to be proved; "
                  "coverage column checks it empirically."]})

    with pd.ExcelWriter(a.out, engine="openpyxl") as xw:
        cfg.to_excel(xw, "Config", index=False)
        df_eq.to_excel(xw, "Equivalence", index=False)
        df_cv.to_excel(xw, "Convergence", index=False)
        df_au.to_excel(xw, "Aumann", index=False)
        df_dp.to_excel(xw, "Dependence", index=False)
        df_tm.to_excel(xw, "Timing", index=False)
        df_sl.to_excel(xw, "Fuzzy_sample_long", index=False)
        df_sw.to_excel(xw, "Fuzzy_sample_wide", index=False)
        df_fi.to_excel(xw, "Fuzzy_interval_long", index=False)
        df_fw.to_excel(xw, "Fuzzy_interval_wide", index=False)

    wb = load_workbook(a.out)
    for ws in wb.worksheets:
        style_sheet(ws)
    add_convergence_chart(wb, df_cv)
    wb.save(a.out)
    print("\nDone ->", a.out)
    print(df_cv.groupby("Scenario")[["Mean_sup_error", "Coverage_of_bound"]].agg(["min", "max"]))


if __name__ == "__main__":
    main()
