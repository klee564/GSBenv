# GSBenv

Group Sparse Bayesian Envelope (GSBenv) — an R package for Bayesian envelope model estimation with group sparsity structure. This package implements MCMC-based Gibbs sampling methods for envelope dimension reduction combined with structured sparsity selection via spike-and-slab priors.

## Overview

Envelope methods reduce the dimension of multivariate regression by identifying the minimal reducing subspace of the response that is material to the regression. **GSBenv** extends this framework by incorporating:

- **Bayesian estimation** via Gibbs sampling (MCMC) for the envelope subspace
- **Group sparsity** through spike-and-slab priors with Beta-Bernoulli group structure
- **Automatic dimension selection** using BIC over a grid of candidate envelope dimensions (`u`) and spike variance parameters (`tau0`)

The package also provides a frequentist sparse envelope method (`spenv`) based on adaptive lasso for comparison.

## Installation

```r
# Install from GitHub using devtools
devtools::install_github("klee564/GSBenv")
```

## Main Functions

| Function | Description |
|---|---|
| `envGS()` | Core Gibbs sampler for the Bayesian sparse envelope model |
| `envGS_selu()` | Wrapper that selects optimal `u` and `tau0` via BIC over a grid search |
| `spenv()` | Frequentist sparse envelope estimation using adaptive lasso |
| `senv_selu()` | Wrapper that selects optimal `u` via BIC for the sparse envelope method |
| `generate_data()` | Generates simulated data from the envelope model |
| `generate_qspar()` / `generate_GSpar()` | Generates true sparse parameters for simulation studies |

## Quick Start

```r
library(GSBenv)

# --- Generate simulation parameters ---
r <- 10   # number of response variables
p <- 20   # number of predictors
u <- 2    # true envelope dimension
q <- 4    # number of sparsity groups

params <- generate_qspar(r, p, u, q)

# --- Simulate data ---
data <- generate_data(params$param, n = 100)

# --- Bayesian Group Sparse Envelope (GSBenv) ---
G <- max(params$groupind)
ex_Betashapes <- matrix(c(1e-6, 1e-6), G + 1, 2, byrow = TRUE)
ex_group_prior <- list(
  ex_groupind = c(rep(G + 1, u), params$groupind),
  ex_Betashapes = ex_Betashapes
)
ex_hyperlist <- list(
  M = diag(1e-6, p),
  B0 = matrix(0, r, p),
  Psi = diag(1, r),
  Psi0 = diag(1, r),
  kappa = 1e-6
)

result_bayes <- envGS_selu(
  Y = data$Y, X = data$X,
  ex_group_prior = ex_group_prior,
  uvec = 1:3, tau0vec = c(0.01, 0.05, 0.1),
  tau1 = 10, ex_hyperlist = ex_hyperlist
)

# --- Frequentist Sparse Envelope ---
result_freq <- senv_selu(Y = data$Y, X = data$X, uvec = 1:3)

# --- Results ---
result_bayes$u     # selected envelope dimension (Bayesian)
result_freq$u      # selected envelope dimension (frequentist)
```

## Running Simulations

`run_simul.R` in the repository root runs a full simulation study comparing the Bayesian group sparse envelope method against the frequentist sparse envelope method across multiple sample sizes.

```r
source("run_simul.R")
# By default, runs simulations with r=10, p=20, q=4
# Sample sizes: n = 50, 100, 200 (repeated 50 times each)
# Results are saved to simulres10_20_4.rds
```

The simulation compares:
- **GSBenv (Bsenv)**: Bayesian group sparse envelope with MCMC
- **Sparse envelope (senv)**: Frequentist sparse envelope with adaptive lasso

and evaluates dimension selection (`u`) and spike variance selection (`tau0`) across different sample sizes.

## Dependencies

- `Renvlp` — envelope methods
- `purrr`, `furrr`, `future` — functional programming and parallel computation
- `caret` — data partitioning utilities
- `magrittr` — pipe operator

## Package Structure

```
GSBenv/
├── R/                   # Source functions
│   ├── envGS.R          # Core Gibbs sampler
│   ├── envGS_MCMC.R     # MCMC routines
│   ├── sel_u.R          # BIC-based model selection (envGS_selu, senv_selu)
│   ├── spenv.R          # Sparse envelope (adaptive lasso)
│   ├── simulsetting.R   # Simulation parameter generation
│   ├── functions.R      # Matrix utilities and parameter transformations
│   └── ...              # Additional helper functions
├── man/                 # Documentation
├── DESCRIPTION
├── NAMESPACE
└── run_simul.R          # Simulation study script
```
