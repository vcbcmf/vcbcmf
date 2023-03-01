# Code usage

The code in this directory reproduces the results of the paper Estimating
Heterogeneous Causal Mediation Effects with Bayesian Decision Tree Ensembles. 
This directory is organized with function definitions in
the `lib/` folder and scripts to reproduce figures in the `src/` directory. The
files in `src/` are organized such that one simply needs to run the files in
alphabetical order to reproduce all of the results in the paper. So, for
example, we should first run `source("src/01_fit_bcmf.R")` to obtain an initial
fit of the varying coefficient BCMF model.

The files `07_01_rlearner_setup.R` to `09_make_fig7_fig8.R` reproduce the
simulation experiments, which are very time consuming. If one just wants to
reproduce the results for the MEPS dataset then they need not run these files.

The directories in the folder are as follows:

- `cache/` contains model fits and other large objects.
- `data/` contains the MEPS data and some other quantities computed by the scripts.
- `figures/` contains pdf figures obtained from some scripts.
- `lib/` contains function definitions.
- `Simulation XYZ` are a cache of simulation results for simulation setting XYZ.
- `src/` contains the scripts for replicating the results.

# Dependencies

To reproduce our code, the following packages should be installed from CRAN using the `install.packages()` function:

- `broom`
- `cowplot`
- `estimatr`
- `extrafont`
- `ggdist`
- `glmnet`
- `latex2exp`
- `matrixStats`
- `Metrics`
- `mgcv`
- `progress`
- `rlearner` from the following GitHub repository: https://github.com/xnie/rlearner
- `rpart`
- `rpart.plot`
- `SoftBart`
- `tikzDevice`
- `tidybayes`
- `tidyverse`
