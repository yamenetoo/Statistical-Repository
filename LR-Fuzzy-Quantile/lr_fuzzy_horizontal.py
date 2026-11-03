#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical verification of the two additions to the article
"Quantile representation of LR random fuzzy numbers":

  (A) Horizontal (Hausdorff / alpha-cut) convergence of the k-knot polygon  [Theorem "horizontal rate",
      Corollaries, Proposition "sharp order"]
        d_H(X_k, X) = sup_alpha max{ |l^k_alpha - l_alpha|, |r^k_alpha - r_alpha| }
                    <= max_j  omega_{Q_j, V_j}( V_j * eps_k(delta) )   with probability >= 1 - delta,
        eps_k(delta) = sqrt(ln(4/delta)/(2k)) + 2/(k+1),
      where omega_{Q,V}(h) = sup{ |Q(x)-Q(y)| : x,y in [0,V], |x-y| <= h } is the modulus of continuity of the
      quantile function on [0,V].  Mechanism: the polygon's left offset at level u = 1-alpha is the linear
      interpolation, at the nodes j/(k+1), of the values H(W_(j)),  H(u) = Q(V u),  W_(j) uniform order statistics.

  (B) Characterisation of the LR class [Theorem "characterisation"]: arms are LR  <=>  multiplicative
      consistency   G_s(t) = G_s(s') G_{s'}(t)  (0 < t <= s' <= s)  AND self-consistency (law of S proportional to F).
      We check the identity numerically for LR arms and for a non-LR family, and recover F from the arms.

Run:  python lr_fuzzy_horizontal.py            (full)      python lr_fuzzy_horizontal.py --quick
      python lr_fuzzy_horizontal.py --out results.json --datadir figdata
Requires: numpy, scipy.
"""
import argparse
import json
import os
import numpy as np
from scipy import stats

SCENARIOS = {
    "A": ("Exp(3) / Beta(3,1)", stats.expon(scale=1 / 3), stats.beta(3, 1)),
    "B": ("Chi2(8) / Chi2(8)", stats.chi2(8), stats.chi2(8)),
    "C": ("Gamma(2,1) / U(0,3)", stats.gamma(2, scale=1), stats.uniform(0, 3)),
}


# ------------------------------------------------------------------ horizontal error of one arm
def arm_horizontal_error(dist, V, k, rng, grid=4001):
    """sup_u |H_k(u) - H(u)|  for one arm, H(u) = Q(V u); polygon = interpolation of H(W_(j)) at j/(k+1)."""
    W = np.sort(rng.random(k))
    nodes = np.arange(k + 2) / (k + 1)
    vals = dist.ppf(V * np.concatenate(([0.0], W, [1.0])))
    u = np.union1d(np.linspace(0, 1, grid), nodes)
    Hk = np.interp(u, nodes, vals)
    H = dist.ppf(V * u)
    return float(np.max(np.abs(Hk - H)))


def modulus(dist, V, h, ngrid=4001):
    """omega_{Q,V}(h) = sup_{x in [0,V]} [Q(min(x+h,V)) - Q(x)]   (Q increasing)."""
    x = np.linspace(0.0, V, ngrid)
    y = np.minimum(x + h, V)
    q = dist.ppf(x)
    return float(np.max(dist.ppf(y) - q))


def eps_k(k, delta, arms=2):
    return np.sqrt(np.log(2 * arms / delta) / (2 * k)) + 2 / (k + 1)


def exp_horizontal(rng, ks, reps, delta=0.05):
    rows = []
    for key, (name, dl, dr) in SCENARIOS.items():
        for k in ks:
            err = np.empty(reps); bnd = np.empty(reps)
            e_k = eps_k(k, delta, arms=2)
            for r in range(reps):
                Vl, Vr = rng.random(2)
                el = arm_horizontal_error(dl, Vl, k, rng)
                er = arm_horizontal_error(dr, Vr, k, rng)
                err[r] = max(el, er)
                bnd[r] = max(modulus(dl, Vl, Vl * e_k), modulus(dr, Vr, Vr * e_k))
            rows.append(dict(scenario=key, laws=name, k=int(k), mean_sup_err=float(err.mean()),
                             median=float(np.median(err)), q95=float(np.quantile(err, 0.95)),
                             sqrt_k_times_median=float(np.sqrt(k) * np.median(err)),
                             coverage_of_bound=float(np.mean(err <= bnd + 1e-12)),
                             median_bound=float(np.median(bnd)),
                             median_ratio_bound_over_err=float(np.median(bnd / err))))
    return rows


def fitted_slopes(rows, kmin=64):
    """log-log slope of the median horizontal error against k (k >= kmin): observed order k^(-slope)."""
    out = {}
    for key in SCENARIOS:
        r = [x for x in rows if x["scenario"] == key and x["k"] >= kmin]
        out[key] = float(-np.polyfit(np.log([x["k"] for x in r]), np.log([x["median"] for x in r]), 1)[0])
    return out


# ------------------------------------------------------------------ sharp order (Proposition)
def exp_sharp(rng, ks, reps):
    """Uniform spread, V = 1: E sup|H_k - H| * sqrt(k) -> E sup|Brownian bridge| = sqrt(pi/2) ln 2."""
    d = stats.uniform(0, 1)
    rows = []
    for k in ks:
        e = np.array([arm_horizontal_error(d, 1.0, k, rng, grid=2 * k + 3) for _ in range(reps)])
        rows.append(dict(k=int(k), mean_err=float(e.mean()), sqrt_k_mean_err=float(np.sqrt(k) * e.mean())))
    return rows, float(np.sqrt(np.pi / 2) * np.log(2))


# ------------------------------------------------------------------ characterisation (Theorem)
def ce_statistic(G, spreads, ngrid=200):
    """max | G_s(t) - G_s(s') G_{s'}(t) |  over s' <= s in `spreads`, t in [0,s']  (zero iff CE holds)."""
    S = np.sort(spreads); best = 0.0
    for a in range(len(S)):
        for b in range(a):
            t = np.linspace(0, S[b], ngrid)
            best = max(best, float(np.max(np.abs(G(S[a], t) - G(S[a], S[b]) * G(S[b], t)))))
    return best


def exp_char(rng, n=15, lam=3.0):
    F = lambda t: 1 - np.exp(-lam * t)
    G_lr = lambda s, t: F(np.minimum(t, s)) / F(s)                              # LR arms
    G_tri = lambda s, t: np.minimum(t, s) / s                                    # triangular = U(0,w)
    G_bad = lambda s, t: (np.minimum(t, s) / s) ** (1.0 + s)                     # exponent depends on s: not LR
    S = rng.exponential(1 / lam, n)
    out = dict(CE_LR_exponential=ce_statistic(G_lr, S), CE_triangular=ce_statistic(G_tri, S),
               CE_non_LR_family=ce_statistic(G_bad, S))
    # identification of F from the arms alone: phi(t) = G_{s0}(t) (t<=s0),  1/G_t(s0) (t>s0),  F = phi/phi(inf)
    s0 = 0.5
    tg = np.linspace(1e-6, 4.0, 4000)
    phi = np.where(tg <= s0, G_lr(s0, tg), 1.0 / G_lr(tg, s0))
    Fhat = phi / phi[-1] * F(tg[-1])           # phi(inf) approximated at the last grid point (F(4)=1-6e-6)
    out["identification_max_err"] = float(np.max(np.abs(Fhat - F(tg))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default="horizontal_results.json")
    ap.add_argument("--datadir", default=None, help="also write whitespace-separated .dat files for pgfplots")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    ks = [4, 16, 64, 256] if a.quick else [4, 8, 16, 32, 64, 128, 256, 512, 1000]
    reps = 100 if a.quick else 1000
    print("(A) horizontal error of the k-knot polygon ...")
    hor = exp_horizontal(rng, ks, reps)
    print("(A') sharp order, uniform spread ...")
    sharp, lim = exp_sharp(rng, [16, 64, 256, 1024] if a.quick else [16, 64, 256, 1024, 4096], 200 if a.quick else 2000)
    print("(B) characterisation ...")
    ch = exp_char(rng)
    slopes = fitted_slopes(hor)
    out = dict(horizontal=hor, fitted_slopes=slopes, sharp=sharp, sharp_limit=lim, characterisation=ch, seed=a.seed, reps=reps)
    json.dump(out, open(a.out, "w"), indent=1)
    for r in hor: print(r)
    print("fitted slopes (observed order k^-slope):", slopes)
    print("sharp (limit %.4f):" % lim); [print(r) for r in sharp]
    print("characterisation:", ch)
    if a.datadir:
        os.makedirs(a.datadir, exist_ok=True)
        for key in SCENARIOS:
            with open(os.path.join(a.datadir, "hor_%s.dat" % key), "w") as f:
                f.write("k err q95 bound\n")
                for r in [x for x in hor if x["scenario"] == key]:
                    f.write("%d %.6f %.6f %.6f\n" % (r["k"], r["median"], r["q95"], r["median_bound"]))
        with open(os.path.join(a.datadir, "sharp.dat"), "w") as f:
            f.write("k val\n")
            for r in sharp: f.write("%d %.6f\n" % (r["k"], r["sqrt_k_mean_err"]))
    print("saved ->", a.out)


if __name__ == "__main__":
    main()
