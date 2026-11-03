#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python implementation of the base paper

    Parchami, A., Grzegorzewski, P., Romaniuk, M. (2024).
    Statistical simulations with LR random fuzzy numbers.  Statistical Papers 65, 3585-3609.

Contents (numbering follows the paper)
--------------------------------------
  Algorithm 1  -> knot_fuzzy_number()        k-knot piecewise linear random fuzzy number [f_C,f_Sl,f_Sr]_k
  Algorithm 2  -> lr_fuzzy_number()           LR random fuzzy number (C,S^l,S^r)_LR (limit of Algorithm 1)
  Algorithm 3  -> fuzzy_random_sample()       simple fuzzy random sample, i.i.d. (f_C,f_Sl,f_Sr)
  Example 1    -> example1()                  2-knot number, the paper's realisation (c=1.717, ...)
  Example 2    -> example2()                  X_k, Y_k, Z_k for k = 2,4,8,16,50,200 and their limits (Fig. 3)
  Example 3    -> example3()                  fuzzy random sample from (N(0,1), U(0,1), Beta(3,1)) (Fig. 4)

Notation.  A realisation is (c, s_l, s_r).  Left arm: mu(c - t) = 1 - F_{S^l|s_l}(t), 0 <= t <= s_l, where
F_{S|s}(t) = F_S(t)/F_S(s) is the c.d.f. truncated to [0, s]  (Corollary 1, eqs. (14)-(15):
L(x) = 1 - F_{S^l|s_l}(x s_l), x in [0,1]; same for R).

Knots (eqs. (6)-(7)), with l_(1)<=...<=l_(k) and r_(1)<=...<=r_(k) the sorted truncated draws:
    left : (c - l_(k-i+1),  i/(k+1)),   right: (c + r_(k-i+1),  i/(k+1)),   i = 1..k.
(The printed eq. (7) has c + r_(i); Example 1 of the paper uses r_(k-i+1), as done here.)

Usage
-----
    python parchami2024_base_paper.py --example all --outdir figures
    python parchami2024_base_paper.py --example 1
Requirements: numpy, scipy, matplotlib.
"""
import argparse
import os
from dataclasses import dataclass

import numpy as np
from scipy import stats


# ------------------------------------------------------------------ truncated distributions
def trunc_sample(dist, s, k, rng):
    """k i.i.d. draws from `dist` truncated to [0, s] (inverse-c.d.f. method)."""
    return dist.ppf(rng.random(k) * dist.cdf(s))


def trunc_cdf(dist, s, t):
    """F_{S|s}(t) = F_S(t)/F_S(s) for 0 <= t <= s (clipped outside)."""
    return dist.cdf(np.clip(t, 0.0, s)) / dist.cdf(s)


# ------------------------------------------------------------------ LR random fuzzy number (Alg. 2)
@dataclass
class LRFuzzyNumber:
    """Realisation x = (c, s_l, s_r)_LR of an LR random fuzzy number (Corollary 1)."""
    c: float
    sl: float
    sr: float
    dl: object      # frozen scipy distribution of S^l
    dr: object      # frozen scipy distribution of S^r

    @property
    def support(self):
        return (self.c - self.sl, self.c + self.sr)

    @property
    def core(self):
        return self.c

    def L(self, x):
        """Left shape function, eq. (14): L(x) = 1 - F_{S^l|s_l}(x s_l), x in [0,1]."""
        return 1.0 - trunc_cdf(self.dl, self.sl, np.asarray(x) * self.sl)

    def R(self, x):
        """Right shape function, eq. (15)."""
        return 1.0 - trunc_cdf(self.dr, self.sr, np.asarray(x) * self.sr)

    def membership(self, x):
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        left = (x > self.c - self.sl) & (x < self.c)
        right = (x > self.c) & (x < self.c + self.sr)
        out[left] = 1.0 - trunc_cdf(self.dl, self.sl, self.c - x[left])
        out[right] = 1.0 - trunc_cdf(self.dr, self.sr, x[right] - self.c)
        out[x == self.c] = 1.0
        return out


def lr_fuzzy_number(dC, dSl, dSr, rng):
    """Algorithm 2: generate C ~ f_C, S^l ~ f_Sl, S^r ~ f_Sr and return (c, s_l, s_r)_LR."""
    return LRFuzzyNumber(float(dC.rvs(random_state=rng)), float(dSl.rvs(random_state=rng)),
                         float(dSr.rvs(random_state=rng)), dSl, dSr)


# ------------------------------------------------------------------ k-knot polygon (Alg. 1)
@dataclass
class KnotFuzzyNumber:
    """k-knot piecewise linear fuzzy number: vertices (x_j, mu_j), sorted by x."""
    x: np.ndarray
    mu: np.ndarray
    c: float
    sl: float
    sr: float
    k: int

    def membership(self, x):
        return np.interp(x, self.x, self.mu, left=0.0, right=0.0)


def knot_from_realisation(c, sl, sr, dSl, dSr, k, rng):
    """Steps of Algorithm 1 for a given realisation (c, s_l, s_r): draw and sort k truncated
    values per arm, build knots (6)-(7) and join them with the core and the support endpoints."""
    l = np.sort(trunc_sample(dSl, sl, k, rng))
    r = np.sort(trunc_sample(dSr, sr, k, rng))
    i = np.arange(1, k + 1)
    left = c - l[::-1]                      # c - l_(k-i+1), i = 1..k   (increasing in i)
    right = c + r[::-1]                     # c + r_(k-i+1)
    x = np.concatenate(([c - sl], left, [c], right[::-1], [c + sr]))
    mu = np.concatenate(([0.0], i / (k + 1), [1.0], (i / (k + 1))[::-1], [0.0]))
    return KnotFuzzyNumber(x, mu, c, sl, sr, k)


def knot_fuzzy_number(k, dC, dSl, dSr, rng):
    """Algorithm 1: simulate a realisation of [f_C, f_Sl, f_Sr]_k."""
    c = float(dC.rvs(random_state=rng))
    sl = float(dSl.rvs(random_state=rng))
    sr = float(dSr.rvs(random_state=rng))
    return knot_from_realisation(c, sl, sr, dSl, dSr, k, rng)


# ------------------------------------------------------------------ fuzzy random sample (Alg. 3)
def fuzzy_random_sample(n, dC, dSl, dSr, rng):
    """Algorithm 3: n i.i.d. LR random fuzzy numbers (Definition 2)."""
    return [lr_fuzzy_number(dC, dSl, dSr, rng) for _ in range(n)]


# ------------------------------------------------------------------ examples and figures
def _plot_mu(ax, fn, lw=1.3, color=None, label=None, n=600):
    lo, hi = fn.support if hasattr(fn, "support") else (fn.c - fn.sl, fn.c + fn.sr)
    xs = np.linspace(lo, hi, n)
    xs = np.sort(np.append(xs, fn.c))
    ax.plot(xs, fn.membership(xs), lw=lw, color=color, label=label)


def example1(outdir, show=False):
    """Example 1: X = [N(5,2^2), Exp(3), U(0,1)]_2 with the realisation of the paper."""
    import matplotlib.pyplot as plt
    c, sl, sr = 1.717, 0.057, 0.186
    # the paper's truncated draws l1=0.028, l2=0.017, r1=0.052, r2=0.156
    k = 2
    l = np.sort([0.028, 0.017]); r = np.sort([0.052, 0.156])
    i = np.arange(1, k + 1)
    left = np.round(c - l[::-1], 3); right = np.round(c + r[::-1], 3)
    print("Example 1: core =", c, " support = [%.3f, %.3f]" % (c - sl, c + sr))
    print("  left knots :", [(float(a), round(float(b), 3)) for a, b in zip(left, i / (k + 1))])
    print("  right knots:", [(float(a), round(float(b), 3)) for a, b in zip(right, i / (k + 1))])
    # paper: (1.689,0.333),(1.700,0.667) and (1.873,0.333),(1.769,0.667)
    assert np.allclose(left, [1.689, 1.700]) and np.allclose(right, [1.873, 1.769])
    x = np.concatenate(([c - sl], c - l[::-1], [c], (c + r[::-1])[::-1], [c + sr]))
    mu = np.concatenate(([0], i / 3, [1], (i / 3)[::-1], [0]))
    fig, ax = plt.subplots(figsize=(5, 3.2))
    ax.plot(x, mu, "o-", lw=1.4)
    ax.set_xlabel("x"); ax.set_ylabel("membership")
    ax.set_title(r"Example 1: $[N(5,2^2),Exp(3),U(0,1)]_2$", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig1_example1.png"), dpi=200)
    if show: plt.show()
    plt.close(fig)


# distributions of Example 2.  Gamma(1,6): shape 1, rate 6 (R convention); change here if needed.
EXAMPLE2 = {
    "X: [N(0,1), Chi2(8), Chi2(8)]": (stats.norm(0, 1), stats.chi2(8), stats.chi2(8)),
    "Y: [Poiss(4), Gamma(1,6), Beta(1,2)]": (stats.poisson(4), stats.gamma(1, scale=1 / 6), stats.beta(1, 2)),
    "Z: [Cauchy(0,1), U(0,1), U(0,3)]": (stats.cauchy(0, 1), stats.uniform(0, 1), stats.uniform(0, 3)),
}


def example2(outdir, seed=2024, ks=(2, 4, 8, 16, 50, 200), show=False):
    """Example 2 / Fig. 3: for each of X, Y, Z one realisation (c, s_l, s_r); rows = k-knot
    polygons for increasing k, last row = the LR limit (Algorithm 2) in red."""
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(seed)
    fig, axes = plt.subplots(len(ks) + 1, 3, figsize=(10, 12), sharex="col")
    for col, (name, (dC, dSl, dSr)) in enumerate(EXAMPLE2.items()):
        c, sl, sr = (float(dC.rvs(random_state=rng)), float(dSl.rvs(random_state=rng)),
                     float(dSr.rvs(random_state=rng)))
        for row, k in enumerate(ks):
            kn = knot_from_realisation(c, sl, sr, dSl, dSr, k, rng)
            ax = axes[row, col]
            ax.plot(kn.x, kn.mu, lw=1.1, color="C0")
            ax.text(0.02, 0.8, "k=%d" % k, transform=ax.transAxes, fontsize=8)
        lim = LRFuzzyNumber(c, sl, sr, dSl, dSr)
        _plot_mu(axes[-1, col], lim, color="red")
        axes[0, col].set_title(name, fontsize=9)
        axes[-1, col].text(0.02, 0.8, "limit", transform=axes[-1, col].transAxes, fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig3_example2.png"), dpi=170)
    if show: plt.show()
    plt.close(fig)


def example3(outdir, seed=2024, n=10, show=False):
    """Example 3 / Fig. 4: fuzzy random sample from (N(0,1), U(0,1), Beta(3,1)).
    Left shape 1 - x (linear); right shape 1 - x^3, as derived in the paper."""
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(seed)
    sample = fuzzy_random_sample(n, stats.norm(0, 1), stats.uniform(0, 1), stats.beta(3, 1), rng)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    for fn in sample:
        _plot_mu(ax, fn, lw=1.1)
    ax.set_xlabel("x"); ax.set_ylabel("membership")
    ax.set_title("Example 3: fuzzy random sample (N(0,1), U(0,1), Beta(3,1))", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig4_example3.png"), dpi=200)
    # check the closed forms of the paper: L(x)=1-x, R(x)=1-x^3
    x = np.linspace(0, 1, 11)
    assert np.allclose(sample[0].L(x), 1 - x) and np.allclose(sample[0].R(x), 1 - x ** 3)
    print("Example 3: L(x)=1-x and R(x)=1-x^3 verified on a grid.")
    if show: plt.show()
    plt.close(fig)


def convergence_check(seed=2024, ks=(2, 4, 8, 16, 50, 200, 1000), reps=200, grid=1000):
    """Mean sup-distance between the k-knot polygon and its LR limit (Theorem 1)."""
    rng = np.random.default_rng(seed)
    dC, dSl, dSr = EXAMPLE2["X: [N(0,1), Chi2(8), Chi2(8)]"]
    print("k      mean sup|X_k - X|")
    for k in ks:
        errs = []
        for _ in range(reps):
            c = float(dC.rvs(random_state=rng)); sl = float(dSl.rvs(random_state=rng)); sr = float(dSr.rvs(random_state=rng))
            kn = knot_from_realisation(c, sl, sr, dSl, dSr, k, rng)
            lim = LRFuzzyNumber(c, sl, sr, dSl, dSr)
            xs = np.linspace(c - sl, c + sr, grid)
            errs.append(np.max(np.abs(kn.membership(xs) - lim.membership(xs))))
        print("%-6d %.4f" % (k, np.mean(errs)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--example", default="all", choices=["1", "2", "3", "convergence", "all"])
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--seed", type=int, default=2024)
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    os.makedirs(a.outdir, exist_ok=True)
    if a.example in ("1", "all"): example1(a.outdir)
    if a.example in ("2", "all"): example2(a.outdir, a.seed)
    if a.example in ("3", "all"): example3(a.outdir, a.seed)
    if a.example in ("convergence", "all"): convergence_check(a.seed)
    print("figures in", a.outdir)


if __name__ == "__main__":
    main()
