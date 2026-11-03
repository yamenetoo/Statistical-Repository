# Base paper: Parchami, Grzegorzewski & Romaniuk (2024)

Python implementation of *Statistical simulations with LR random fuzzy numbers* (Statistical Papers 65, 3585-3609).

| Paper | Function in `parchami2024_base_paper.py` |
|---|---|
| Algorithm 1 (k-knot polygon) | `knot_fuzzy_number`, `knot_from_realisation` |
| Algorithm 2 (LR random fuzzy number) | `lr_fuzzy_number`, class `LRFuzzyNumber` (shapes L, R of eqs. (14)-(15)) |
| Algorithm 3 (fuzzy random sample) | `fuzzy_random_sample` |
| Example 1 / Fig. 1 | `example1` (reproduces the knots (1.689,0.333), (1.700,0.667), (1.873,0.333), (1.769,0.667)) |
| Example 2 / Fig. 3 | `example2` |
| Example 3 / Fig. 4 | `example3` (checks L(x)=1-x, R(x)=1-x^3) |

```bash
pip install numpy scipy matplotlib
python parchami2024_base_paper.py --example all --outdir figures
```

Notes: the right-arm knots use `c + r_(k-i+1)` as in Example 1 of the paper (the printed eq. (7) has `r_(i)`);
Gamma(1,6) in Example 2 is read as shape 1, rate 6 (change `EXAMPLE2` if the scale convention is different).
The quantile-based generator proposed in our article is in `../lr_fuzzy_simulation.py`.
