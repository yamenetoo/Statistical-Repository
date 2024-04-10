This repository is dedicated to the LG BootStrap Control Chart. [Prof Dr Ali Akbar Heydari](https://scholar.google.com/citations?user=68RAHCoAAAAJ&hl=en).
# Function Hyperparameters

- `u`: Quantile value used in the analysis.
- `alpha`: Significance level or Type I error rate, controlling parameter in Control Chart analysis. (Also known as "significance level" or "alpha level".)
- `theta`, `p`: Parameters of the probability distribution being analyzed.
- `n`: Sample size used in the analysis. Default is 5.
- `BootStrap_iterations`: Number of iterations for Bootstrap resampling. Default is 1E+4.
- `MontiCarlo_iter`: Number of iterations for Monte Carlo simulation. Default is 30.
- `ARL_repetition`: Number of repetitions/generations for calculating Average Run Length (ARL). Default is 100.
- `verbos`: Boolean parameter controlling whether to display results during fitting. Default is FALSE.

# Function Output

Output is a list with the following components:
- `MLCL`: Mean of the Lower Control Limit (LCL).
- `MUCL`: Mean of the Upper Control Limit (UCL).
- `SDLCL`: Standard Deviation of the Lower Control Limit (LCL).
- `SDUCL`: Standard Deviation of the Upper Control Limit (UCL).
- `MARL`: Mean Average Run Length (ARL) across generations.
- `SDARL`: Standard Deviation of Average Run Length (ARL) across generations.

# Citation

Professor Ali Akbar Heydari, University of Tabriz · Department of Statistics.
