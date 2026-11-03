#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inference experiments of Section 4 / Section 5.6 of the article
"Quantile representation of LR random fuzzy numbers".

Setting: i.i.d. LR random fuzzy numbers with Exp(lambda) spreads (lambda = 3), observed without noise.

  recover  : Proposition 4.1  - the spread c.d.f. is recovered from the membership function of the fuzzy
             number with the largest spread (error ~ 1/n), compared with the crisp empirical c.d.f. (~ 1/sqrt(n)).
  compat   : Lemma 4.3        - the necessary 'rescaled truncation' identity
             1 - mu_i(t) = (1 - mu_i(S_j)) (1 - mu_j(t)),  S_j <= S_i,  holds for LR arms (and, as expected, for
             triangular arms, which are the uniform-spread case): it is necessary, not sufficient, for H0, which is
             why the goodness-of-fit test is needed in addition.
  gof      : Proposition 4.2  - Kolmogorov-type test D_n = sup |F_n - Phi_n|; size under H0 and power against
             (i) triangular arms and (ii) a misspecified spread law (Gamma(2, 1/6), same mean).
  info/mle : Proposition 4.4  - Fisher information of k knots per arm; MLE standard deviation vs 1/sqrt(n I_k).

Run:  python lr_fuzzy_inference.py                 (full run, several minutes)
      python lr_fuzzy_inference.py --quick         (fast test)
      python lr_fuzzy_inference.py --only gof --out gof.json
Requires: numpy, scipy.
"""
import argparse
import json
import numpy as np
from scipy import integrate, optimize, stats

LAM = 3.0
F = lambda t: 1 - np.exp(-LAM * t)
KS_CRIT_5 = 1.3581          # asymptotic 5% critical value of the Kolmogorov distribution


# ----------------------------------------------------------------------------- recovery (Prop. 4.1)
def exp_recover(rng, ns=(10, 50, 250, 1000), R=2000):
    tg = np.linspace(0, 12 / LAM, 4001)
    rows = []
    for n in ns:
        eA = np.empty(R); eB = np.empty(R)
        for r in range(R):
            s = rng.exponential(1 / LAM, n)
            smax = s.max(); vmax = F(smax)
            mu = np.where(tg <= smax, 1 - F(np.minimum(tg, smax)) / vmax, 0.0)   # membership of the largest-spread number
            Phi = np.where(tg <= smax, (1 - mu) * n / (n + 1), n / (n + 1))
            eA[r] = np.max(np.abs(Phi - F(tg)))
            ss = np.sort(s)
            eB[r] = np.max(np.abs(np.searchsorted(ss, tg, side="right") / n - F(tg)))
        g = lambda v: max(abs(n / (n + 1) - v), 1 / (n + 1)) * n * v ** (n - 1)   # V_(n) ~ Beta(n,1)
        exact = integrate.quad(g, 0, 1, points=[n / (n + 1)], limit=400)[0]
        rows.append(dict(n=n, err_fuzzy=float(eA.mean()), err_crisp_edf=float(eB.mean()),
                         n_times_err_fuzzy=float(eA.mean() * n), sqrt_n_times_err_crisp=float(eB.mean() * np.sqrt(n)),
                         exact_expectation=float(exact), ratio=float(eB.mean() / eA.mean())))
    print("limit of n*E[err] = 1+exp(-2) = %.4f ;  sqrt(pi/2) ln2 = %.4f" % (1 + np.exp(-2), np.sqrt(np.pi / 2) * np.log(2)))
    return rows


# ----------------------------------------------------------------------------- compatibility (Lemma 4.3)
def _T_stat(S, mu):
    S = np.sort(S); best = 0.0
    for a in range(len(S)):
        for b in range(a):
            t = np.linspace(0, S[b], 200)
            lhs = 1 - mu(S[a], t)
            rhs = (1 - mu(S[a], np.array([S[b]]))) * (1 - mu(S[b], t))
            best = max(best, float(np.max(np.abs(lhs - rhs))))
    return best


def exp_compat(rng, n=12):
    S = rng.exponential(1 / LAM, n)
    mu_lr = lambda s, t: np.clip(1 - F(np.minimum(t, s)) / F(s), 0, 1)
    mu_tri = lambda s, t: np.clip(1 - t / s, 0, 1)
    return dict(n=n, T_lr=_T_stat(S, mu_lr), T_triangular=_T_stat(S, mu_tri))


# ----------------------------------------------------------------------------- GOF test (Prop. 4.2)
def _gof_once(n, scen, rng):
    if scen == "H0":      # LR arms from F, spreads Exp(3)
        s = rng.exponential(1 / LAM, n); mu = lambda t, sm: 1 - F(np.minimum(t, sm)) / F(sm)
    elif scen == "triangular":
        s = rng.exponential(1 / LAM, n); mu = lambda t, sm: 1 - np.minimum(t, sm) / sm
    elif scen == "misspecified":   # arms from Exp(3), spreads Gamma(2, scale = 1/(2 lambda))
        s = rng.gamma(2, 1 / (2 * LAM), n); mu = lambda t, sm: 1 - F(np.minimum(t, sm)) / F(sm)
    smax = s.max()
    Phi = lambda t: np.where(t <= smax, (1 - mu(t, smax)) * n / (n + 1), n / (n + 1))
    U = Phi(np.sort(s)); i = np.arange(1, n + 1)
    return max(np.max(i / n - U), np.max(U - (i - 1) / n))


def exp_gof(rng, ns=(20, 50, 100, 250), R=3000):
    rows = []
    for n in ns:
        row = dict(n=n)
        for sc in ("H0", "triangular", "misspecified"):
            D = np.array([_gof_once(n, sc, rng) for _ in range(R)])
            row[sc + "_rejection"] = float(np.mean(np.sqrt(n) * D > KS_CRIT_5))
            row[sc + "_mean_sqrt_n_D"] = float(np.mean(np.sqrt(n) * D))
        rows.append(row)
    return rows


# ----------------------------------------------------------------------------- Fisher information / MLE (Prop. 4.4)
def fisher_info(k):
    """I_k = I_s + k J for Exp(lambda): I_s = 1/lambda^2, J = E Var(l | l < S)."""
    def var_trunc(s):
        e = np.exp(-LAM * s)
        return 1 / LAM ** 2 - s ** 2 * e / (1 - e) ** 2
    J = integrate.quad(lambda s: var_trunc(s) * LAM * np.exp(-LAM * s), 1e-9, 60)[0]
    return 1 / LAM ** 2 + k * J, J


def _mle(s, L):
    n = len(s); k = L.shape[1]
    def nll(l):
        v = -(n * np.log(l) - l * s.sum())
        if k > 0:
            v += -(n * k * np.log(l) - l * L.sum() - k * np.sum(np.log1p(-np.exp(-l * s))))
        return v
    return optimize.minimize_scalar(nll, bounds=(0.05, 50), method="bounded", options=dict(xatol=1e-9)).x


def exp_mle(rng, ns=(20, 100), ks=(0, 2, 5, 20), R=1500):
    rows = []
    for n in ns:
        for k in ks:
            est = np.empty(R)
            for r in range(R):
                s = rng.exponential(1 / LAM, n)
                if k > 0:
                    u = rng.random((n, k)); L = -np.log1p(-u * (1 - np.exp(-LAM * s))[:, None]) / LAM
                else:
                    L = np.empty((n, 0))
                est[r] = _mle(s, L)
            I, J = fisher_info(k)
            rows.append(dict(n=n, k=k, mean=float(est.mean()), sd=float(est.std(ddof=1)),
                             sd_theory=float(1 / np.sqrt(n * I)), efficiency_gain=float(I * LAM ** 2)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--quick", action="store_true", help="small replication numbers (testing)")
    ap.add_argument("--only", choices=["recover", "compat", "gof", "mle"], default=None)
    ap.add_argument("--out", default="inference_results.json")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    q = a.quick
    out = {}
    if a.only in (None, "recover"):
        print("recovery ..."); out["recover"] = exp_recover(rng, ns=(10, 50) if q else (10, 50, 250, 1000), R=100 if q else 2000)
    if a.only in (None, "compat"):
        print("compatibility ..."); out["compat"] = exp_compat(rng)
    if a.only in (None, "gof"):
        print("GOF test ...");      out["gof"] = exp_gof(rng, ns=(20, 50) if q else (20, 50, 100, 250), R=200 if q else 3000)
    if a.only in (None, "mle"):
        print("Fisher information / MLE ...")
        out["mle"] = exp_mle(rng, ns=(20,) if q else (20, 100), ks=(0, 2) if q else (0, 2, 5, 20), R=100 if q else 1500)
        out["fisher_info_per_knot_J"] = fisher_info(1)[1]
    with open(a.out, "w") as f:
        json.dump(out, f, indent=1)
    for k, v in out.items():
        print("\n==", k)
        for row in (v if isinstance(v, list) else [v]):
            print(row)
    print("\nsaved ->", a.out)


if __name__ == "__main__":
    main()
