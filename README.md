# Importance Sampling with Variable Selection for Reliability Analysis

This repository contains the code to reproduce the numerical experiments in the paper.

## Repository Structure

```
IS_VS/
├── R_methods/                      # R implementation (all methods except ECL)
│   ├── main.R                      # Main runner script
│   ├── install_packages.R          # Package installation script
│   ├── AMISE_bivariate.R           # Bandwidth optimization (IS methods)
│   ├── AMISE_univariate.R          # Bandwidth optimization (IS methods)
│   ├── AMISE_univariate_OptiTree.R # Bandwidth optimization (OptiTree)
│   ├── Additive_S.R                # Additive model functions
│   ├── function_OptiTreeStrata_independent.R
│   ├── function_OptiTreeStrata_dependent.R
│   ├── function_OptiTreeStrata_highdim.R
│   ├── Ex1/                        # Example 1 (D=4)
│   │   ├── functions_ex1.R
│   │   └── methods/
│   │       ├── ex1_IS_VS.R
│   │       ├── ex1_IS_CE.R
│   │       ├── ex1_IS_Pareto.R
│   │       ├── ex1_WAMK_SIS.R
│   │       ├── ex1_Lasso.R
│   │       ├── ex1_SpAM.R
│   │       ├── ex1_RF_RFE.R
│   │       └── ex1_OptiTreeStrata.R
│   ├── Ex2/                        # Example 2 (D=5, correlated)
│   ├── Ex3/                        # Example 3 (D=10)
│   └── Ex4/                        # Example 4 (D=50)
│
├── Python_ECL/                     # Python implementation (ECL baseline)
│   ├── ecl_example1.ipynb          # ECL for Example 1
│   ├── ecl_example2.ipynb          # ECL for Example 2
│   ├── ecl_example3.ipynb          # ECL for Example 3
│   ├── ecl_method.py               # ECL method implementation
│   ├── eclGP.py                    # Gaussian Process for ECL
│   ├── MFIS_Simulation.py          # MFIS simulation functions
│   ├── psis.py                     # Pareto Smoothed Importance Sampling
│   ├── data/                       # Data folder
│   └── mfis/                       # MFIS helper functions
│
└── README.md                       # This file
```

---

## Requirements

### R Dependencies

R version: 4.3.0 or later

| Package | Version | Description |
|---------|---------|-------------|
| doParallel | 1.0.17 | Parallel backend for foreach |
| foreach | 1.5.2 | Parallel looping |
| truncnorm | 1.0-9 | Truncated normal distribution |
| rootSolve | 1.8.2.4 | Nonlinear root finding |
| rmutil | 1.1.10 | Utilities for regression models |
| pracma | 2.4.4 | Numerical math functions |
| mvtnorm | 1.2-4 | Multivariate normal distribution |
| cubature | 2.1.0 | Adaptive multivariate integration |
| caret | 6.0-94 | Machine learning utilities |
| lhs | 1.1.6 | Latin Hypercube Sampling |
| LHD | 1.3.3 | Orthogonal Array LHS |
| R.utils | 2.12.3 | Programming utilities |
| crch | 1.1-2 | Censored regression |
| loo | 2.6.0 | Leave-one-out cross-validation |
| GauPro | 0.2.12 | Gaussian process fitting |
| DiceKriging | 1.6.0 | Kriging methods |
| SAM | 1.1.2 | Sparse Additive Models |
| glmnet | 4.1-8 | Lasso regression |
| randomForest | 4.7-1.1 | Random forest |

### Python Dependencies

Python version: 3.10 or later

| Package | Version | Description |
|---------|---------|-------------|
| numpy | 1.24.3 | Numerical computing |
| pandas | 2.0.3 | Data manipulation |
| pyDOE | 0.3.8 | Design of Experiments (LHS) |
| scipy | 1.11.1 | Scientific computing |
| scikit-learn | 1.3.0 | Machine learning (Gaussian Process) |
| matplotlib | 3.7.2 | Plotting |
| dill | 0.3.7 | Extended pickling |

---

## Installation

### R Setup

```r
# Navigate to R_methods folder
setwd("path/to/R_methods")

# Run the installation script
source("install_packages.R")
```

Or install manually:

```r
install.packages(c(
  "doParallel", "foreach", "truncnorm", "rootSolve", "rmutil",
  "pracma", "mvtnorm", "cubature", "caret", "lhs", "LHD",
  "R.utils", "crch", "loo", "GauPro", "DiceKriging",
  "SAM", "glmnet", "randomForest"
))
```

### Python Setup

```bash
# Navigate to Python_ECL folder
cd path/to/Python_ECL

# Install dependencies
pip install numpy pandas pyDOE scipy scikit-learn matplotlib dill
```

---

## Reproducing Results

### Table 1: R Methods (All methods except ECL)

```r
# Set working directory to R_methods folder
setwd("path/to/R_methods")

# Source main script
source("main.R")

# Run all methods on all examples (100 replications, 25 cores)
all_results <- run_all_examples(n_reps = 100, ncores = 25)

# Or run individual examples/methods:
result_ex1 <- run_method(1, "IS-VS", n_reps = 100, ncores = 25)
result_ex2 <- run_method(2, "OptiTree", n_reps = 100, ncores = 25)
```

**Available methods:** IS-VS, IS-CE, IS-Pareto, WAMK-SIS, Lasso, SpAM, RF-RFE, OptiTree

**Estimated run time:** ~14 hours total (25 cores, 100 replications)

### Table 1: ECL Baseline (Python)

Open and run the Jupyter notebooks in the `Python_ECL` folder:

1. `ecl_example1.ipynb` - ECL for Example 1 (D=4)
2. `ecl_example2.ipynb` - ECL for Example 2 (D=5)
3. `ecl_example3.ipynb` - ECL for Example 3 (D=10)

Each notebook is self-contained and can be run independently.

**Estimated run time:** ~2 hours per example

---

## Numerical Examples

| Example | Dimension (D) | Threshold (l) | Description |
|---------|---------------|---------------|-------------|
| Example 1 | 4 | 18.99 | Independent inputs |
| Example 2 | 5 | 24.2 | Correlated inputs |
| Example 3 | 10 | 18.74 | Higher dimension |
| Example 4 | 50 | 18.42 | High dimension |

---

## Method Availability

| Method | Ex1 | Ex2 | Ex3 | Ex4 |
|--------|-----|-----|-----|-----|
| IS-VS | ✓ | ✓ | ✓ | ✓ |
| IS-CE | ✓ | ✓ | ✓ | ✓ |
| IS-Pareto | ✓ | ✓ | ✓ | ✓ |
| WAMK-SIS | ✓ | ✓ | ✓ | ✓ |
| Lasso | ✓ | ✓ | ✗ | ✓ |
| SpAM | ✓ | ✓ | ✓ | ✓ |
| RF-RFE | ✓ | ✓ | ✓ | ✓ |
| OptiTree | ✓ | ✓ | ✓ | ✓ |
| ECL | ✓ | ✓ | ✓ | ✗ |

---

## Notes

- Random seeds are set for reproducibility. Default seed_base = 123.
- Numerical results may have minor variations (< 1% relative difference) due to parallel execution order.
- Run times are estimated based on a 25-core machine. Adjust `ncores` parameter based on your system.

---

## License

MIT License

Copyright (c) 2025

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
