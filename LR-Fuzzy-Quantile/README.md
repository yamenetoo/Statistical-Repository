# LR-Fuzzy-Quantile

Simulation code for the article **"Quantile representation of LR random fuzzy numbers"**
(a contribution to Parchami, Grzegorzewski & Romaniuk (2024), *Statistical simulations with LR random fuzzy numbers*).

Alpha-cuts of the limiting LR random fuzzy number: `[C - Q_l((1-a)V_l), C + Q_r((1-a)V_r)]`, with `V_j ~ U(0,1)`.

## Files
- `lr_fuzzy_simulation.py` - all experiments (equivalence, convergence rate, Aumann expectation, dependence, speed, fuzzy samples); writes an Excel file.
- `lr_fuzzy_simulation.ipynb` - the same code as a notebook (open in Google Colab: *Runtime -> Run all*).

## Run
```bash
pip install -r requirements.txt
python lr_fuzzy_simulation.py                 # default: 500 replications, seed 2026
python lr_fuzzy_simulation.py --quick         # fast test
```
