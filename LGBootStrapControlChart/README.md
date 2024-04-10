This repository is dedicated to the LG BootStrap Control Chart. 
# LG BootStrapControlChart
## Function Hyperparameters

- `u`: Quantile value used in the analysis.
- `alpha`: Significance level or Type I error rate, controlling parameter in Control Chart analysis. (Also known as "significance level" or "alpha level".)
- `theta`, `p`: Parameters of the probability distribution being analyzed.
- `n`: Sample size used in the analysis. Default is 5.
- `BootStrap_iterations`: Number of iterations for Bootstrap resampling. Default is 1E+4.
- `MontiCarlo_iter`: Number of iterations for Monte Carlo simulation. Default is 30.
- `ARL_repetition`: Number of repetitions/generations for calculating Average Run Length (ARL). Default is 100.
- `verbos`: Boolean parameter controlling whether to display results during fitting. Default is FALSE.

## Function Output

Output is a list with the following components:
- `MLCL`: Mean of the Lower Control Limit (LCL).
- `MUCL`: Mean of the Upper Control Limit (UCL).
- `SDLCL`: Standard Deviation of the Lower Control Limit (LCL).
- `SDUCL`: Standard Deviation of the Upper Control Limit (UCL).
- `MARL`: Mean Average Run Length (ARL) across generations.
- `SDARL`: Standard Deviation of Average Run Length (ARL) across generations.
# MontiCarloSimulation
This R script is conducting a Monte Carlo simulation to evaluate the performance of a control chart using the `LGBootStrapControlChart` function from the sourced file `LGBootStrapControlChart.R`.

this Monte Carlo simulation works as follows:

1. **Setup and Parameters**: 
   - `n_values`, `u_values`, `p_values`, and `theta_values` are vectors defining the values of different parameters to be tested.
   - `LGBootStrapControlChart` function requires parameters such as `u`, `alpha`, `theta`, `p`, `n`, `BootStrap_iterations`, `MontiCarlo_iter`, `ARL_rep`, and `verbos`. These parameters control various aspects of the control chart analysis.

2. **Nested Loops**:
   - Four nested loops iterate through different combinations of `n`, `u`, `p`, and `theta` values. 
   - For each combination, the `LGBootStrapControlChart` function is called with the specified parameters.
   - Inside the innermost loop, the function's result is captured, and information about the parameters and results is printed (`cat()` function).
   
3. **Data Collection**:
   - For each iteration, the results of the control chart analysis, including `MLCL`, `MUCL`, `SDLCL`, `SDUCL`, `MARL`, and `SDARL`, are stored in a data frame called `results_df`.
   - This data frame accumulates results from all iterations.

4. **Writing Results to File**:
   - After all iterations are completed, the `results_df` data frame is written to a CSV file named "res.csv" using the `write.csv()` function.

5. **Output**:
   - The output CSV file contains columns representing the different parameters (`u`, `p`, `theta`) and the corresponding control chart analysis results (`MLCL`, `MUCL`, `SDLCL`, `SDUCL`, `MARL`, `SDARL`).

In summary, this script systematically explores various combinations of parameters for the control chart analysis using a Monte Carlo simulation approach. It captures the results of each simulation run and aggregates them into a CSV file for further analysis and interpretation.


## Example
Load the LGBootStrapControlChart function from the sourced file
source("LGBootStrapControlChart.R")

Example of using LGBootStrapControlChart function:
Perform a control chart analysis with the following parameters:
- Quantile value (u) = 0.99
- Significance level (alpha) = 0.07
- Parameter of the probability distribution (theta) = 2
- Another parameter of the probability distribution (p) = 0.8
- Sample size (n) = 5
- Number of iterations for Bootstrap resampling = 10000
- Number of iterations for Monte Carlo simulation = 20
- Number of repetitions for calculating Average Run Length (ARL) = 1000
- Verbose mode enabled (verbos = TRUE), to display results during fitting
result <- LGBootStrapControlChart(
  u = 0.99, alpha = 0.07,
  theta = 2, p = 0.8, n = 5,
  BootStrap_iterations = 10000,
  MontiCarlo_iter = 20,
  ARL_rep = 1000, verbos = TRUE
)



# Citation

[Prof Dr Ali Akbar Heydari](https://scholar.google.com/citations?user=68RAHCoAAAAJ&hl=en),
University of Tabriz ·
Department of Statistics.
