# DistFit Framework

DistFit is a MATLAB-based toolkit for fitting and ranking probability distributions against observed data. It streamlines exploratory analysis across a broad library of continuous and discrete models, computes goodness-of-fit metrics, and visualises results with ready-made plots and demos.

## Key Features
- Fit and rank 30+ probability distributions with Kolmogorov-Smirnov, Anderson-Darling, Chi-square, log-likelihood, AIC, and BIC metrics.
- Automatic model selection for continuous or discrete samples.
- Robust parameter estimation with fallback optimisers to avoid convergence warnings.
- Plotting helpers for histograms, empirical CDF overlays, and Q-Q diagnostics.
- Scenario-driven demos covering exploratory analysis, reliability, finance, discrete demand, simulation, extreme value (Gumbel, GEV), and threshold exceedances (GPD).

## Requirements
- MATLAB R2017b or newer (Statistics and Machine Learning Toolbox recommended but optional).
- The `matlab` directory must be on the MATLAB path (`addpath`).

## Quick Start
1. Open MATLAB and change directory to the repository root.
2. Add the MATLAB sources to your path:
   ```matlab
   addpath(fullfile(pwd, 'matlab'));
   addpath(fullfile(pwd, 'matlab', 'utils'));
   addpath(fullfile(pwd, 'matlab', 'gof'));
   addpath(fullfile(pwd, 'matlab', 'plot'));
   ```
3. Run the standard regression check:
   ```matlab
   run_demo
   ```
   This fits the bundled datasets in `data/` and prints the top-ranked distributions per file.

## Demo Gallery
Each demo script is located in `matlab/demo/`. Run any script directly in MATLAB after adding the repository to your path. Every demo prints a ranked table, exports key variables to the base workspace, and opens informative plots.

| Demo | Scenario | Highlights |
| --- | --- | --- |
| `demo_exploratory` | Continuous exploratory analysis | Mixed normal samples, histogram + CDF overlay |
| `demo_reliability` | Accelerated life test data | Survival curve comparison, Weibull/Gamma fits |
| `demo_finance` | Heavy-tailed returns | Student's t vs logistic diagnostics, Q-Q plot |
| `demo_demand` | Discrete demand counts | PMF overlay with Poisson/negative binomial |
| `demo_simulation` | Service time simulation | Synthetic sampling with best-fit RNG |
| `demo_gumbel` | Block maxima (Type I EVT) | Density and Q-Q checks |
| `demo_weibull` | Weibull life analysis | Empirical CDF vs fitted CDF |
| `demo_gev` | Generalized extreme value | Density + CDF overlays for block maxima |
| `demo_gpd` | Threshold exceedances | Tail survival comparison |

All demos call `setup_demo_paths` to configure the MATLAB path, so no additional setup is required beyond loading the repository.

## Distribution Library
Distributions are defined under `matlab/+dists/`. Each model exposes:
- `name`, `type` (`cont` or `disc`), parameter metadata, and support limits.
- `start` (initial parameter guesses) and `fit` functions with robust optimisation.
- Probability functions: `pdf`, `cdf`, and optional `rnd` for sampling.

The registry (`matlab/+dists/registry.m`) enumerates every model used by `fitAllDistributions`. Newly added GEV and Generalized Pareto implementations include parameter transforms and non-toolbox fallbacks for wider compatibility.

## Adding a New Distribution
1. Create a function in `matlab/+dists/` that returns a struct following the existing pattern.
2. Implement `start`, `fit`, and probability functions, ensuring parameters stay within their valid domain.
3. Register the distribution in `matlab/+dists/registry.m`.
4. (Optional) Add a demo under `matlab/demo/` that showcases the new model.

## Utilities
- `fitAllDistributions`: core engine that filters models by data type (continuous vs discrete), fits parameters, computes GOF metrics, and returns ranked results.
- Plot helpers in `matlab/plot/` (histograms with PDFs, ECDF overlays, Q-Q plots).
- Demo utilities in `matlab/demo/` for consistent table printing, model lookup, and path setup.

## Author
Ahmad Faisal Mohamad Ayob, VSG Labs, 2025.

