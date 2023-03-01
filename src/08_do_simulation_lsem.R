## Load ----

library(tidyverse)

source('lib/lsem_mediate.R')
source('lib/lsem_boot.R')
source('lib/do_simulation_lsem.R')

NUM_SEED <- 200
shuffle <- sample(1:NUM_SEED)

## MEPS train and test set ----

meps <- readRDS('data/meps.rds')
n <- nrow(meps)
set.seed(12345)
i_train <- sample(1:nrow(meps), floor(n/2))
i_test <- c(1:nrow(meps))[-i_train]

## True effects with LSEM ----

formula_m_lsem <- phealth ~ smoke * (age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt)
formula_y_lsem <- logY ~ (smoke + phealth) * (age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt)
formula_X_lsem <- phealth ~ age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt

fit_m_lsem <- lm(formula_m_lsem, data = meps)
fit_y_lsem <- lm(formula_y_lsem, data = meps)

out_lsem <- lsem_mediate(
  fit_m = fit_m_lsem, 
  fit_y = fit_y_lsem, 
  formula_X = formula_X_lsem, 
  data_test = meps, 
  mediator_name = "phealth", 
  treat_name = "smoke"
)

zeta_hat_lsem <- out_lsem$zeta
delta_hat_lsem <- out_lsem$delta
sigma_m_hat_lsem <- out_lsem$sigma_m
sigma_y_hat_lsem <- out_lsem$sigma_y
rm(out_lsem)
gc()

## True effects with BART ----

out_bart <- readRDS('cache/01_fit_bcmf_out_bart.rds')
mu_y_hat_bart    <- colMeans(out_bart$mu_y_samples)
zeta_hat_bart    <- colMeans(out_bart$zeta_samples)
d_hat_bart       <- colMeans(out_bart$d_samples)
mu_m_hat_bart    <- colMeans(out_bart$mu_m_samples)
tau_hat_bart     <- colMeans(out_bart$tau_samples)
sigma_y_hat_bart <- mean(out_bart$sigma_y_samples)
sigma_m_hat_bart <- mean(out_bart$sigma_m_samples)
rm(out_bart)
gc()

## Run Simulation ----

set.seed(1)
seeds <- sample.int(10e6, NUM_SEED)
simulation_lsem <- do_simulation_lsem(
  data             = meps,
  i_train          = i_train,
  i_test           = i_test,
  fit_m_lsem       = fit_m_lsem,
  fit_y_lsem       = fit_y_lsem,
  formula_X_lsem   = formula_X_lsem,
  formula_m_lsem   = formula_m_lsem,
  formula_y_lsem   = formula_y_lsem,
  mediator_name    = 'phealth',
  outcome_name     = 'logY',
  treat_name       = 'smoke',
  zeta_hat_lsem    = zeta_hat_lsem,
  delta_hat_lsem   = delta_hat_lsem,
  sigma_m_hat_lsem = sigma_m_hat_lsem,
  sigma_y_hat_lsem = sigma_y_hat_lsem,
  mu_y_hat_bart    = mu_y_hat_bart,
  zeta_hat_bart    = zeta_hat_bart,
  d_hat_bart       = d_hat_bart,
  mu_m_hat_bart    = mu_m_hat_bart,
  tau_hat_bart     = tau_hat_bart,
  sigma_y_hat_bart = sigma_y_hat_bart,
  sigma_m_hat_bart = sigma_m_hat_bart,
  n_reps           = NUM_SEED,
  seeds            = seeds[shuffle]
) 
