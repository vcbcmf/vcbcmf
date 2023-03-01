## Load ----

library(SoftBart)
library(tidyverse)
library(matrixStats)
library(rpart)
library(progress)
library(glmnet)
library(rlearner)

source('lib/get_clever_cov.R')
source('lib/get_ps.R')
source('lib/bart_mediate.R')
source('lib/lsem_mediate.R')
source('lib/do_simulation_bart.R')

meps <- readRDS('data/meps.rds')

## Pooling ----

meps_pi <- meps %>%
  mutate(race = ifelse(race == "PacificIslander", "Asian", race))

## Getting design matrix ----

formula <- ~age + bmi + edu + log(income + 1000) +
  log(povlev + 100) + region + race +
  sex + marital + seatbelt - 1

X_meps <- model.matrix(formula, data = meps_pi)
# X_meps_2 <- model.matrix(~age + bmi + edu + log(income + 1000) + 
#                            log(povlev + 100) + region + race + 
#                            sex + marital + seatbelt, data = meps_pi)
# 

## Getting rlasso simulation objects ----

set.seed(112233)
rlearn_m <- rlasso(x = X_meps, w = meps$smoke, y = meps$phealth, k_folds = 30)

set.seed(112233)
rlearn_total <- rlasso(
  x = cbind(X_meps),
  y = meps$logY,
  w = meps$smoke,
  k_folds = 30
)

set.seed(112233)
rlearn_indirect <- rlasso(
  x = cbind(X_meps, meps$smoke),
  y = meps$logY,
  w = meps$phealth,
  k_folds = 30
)

set.seed(112233)
rlearn_direct <- rlasso(
  x = cbind(X_meps, meps$phealth),
  y = meps$logY,
  w = meps$smoke,
  k_folds = 30
)

set.seed(112233)
rlearn_test <- rlasso(
  x = cbind(X_meps),
  y = meps$logY,
  w = meps$phealth,
  k_folds = 30
)

rlearn_fits <- list(
  fit_m = rlearn_m,
  fit_tau = rlearn_total,
  fit_delta = rlearn_indirect,
  w = meps$smoke,
  tau = rlearn_total$tau_hat %>% as.numeric(),
  delta = (rlearn_indirect$tau_hat * rlearn_m$tau_hat) %>% as.numeric(),
  sigma_y = 1.5,
  sigma_m = 1
)

rlearn_fits$zeta <- rlearn_fits$tau - rlearn_fits$delta
saveRDS(object = rlearn_fits, file = "cache/07_01_rlearn_fits.rds")
