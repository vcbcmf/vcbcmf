## Load ----

library(SoftBart)
library(tidyverse)
library(matrixStats)
library(rpart)
library(progress)

NUM_SEED <- 200
shuffle <- sample(1:NUM_SEED)

source('lib/get_clever_cov.R')
source('lib/get_ps.R')
source('lib/bart_mediate.R')
source('lib/lsem_mediate.R')
source('lib/do_simulation_bart.R')
source('lib/simulate_rlearn_med.R')
source("lib/do_simulation_rlearn.R")
source('lib/rlasso_mediate.R')

set.seed(12345)
meps          <- readRDS('data/meps.rds')
rlearner_fits <- readRDS("cache/07_01_rlearn_fits.rds")
n             <- nrow(meps)
i_train       <- sample(1:nrow(meps), floor(n/2))
i_test        <- c(1:nrow(meps))[-i_train]

## A Bunch of Formulas

formula_m_bart  <- phealth ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
formula_y_bart  <- logY ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt + phealth
formula_ps_bart <- smoke ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt

formula_m_lsem <- phealth ~ smoke * (age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt)
formula_y_lsem <- logY ~ (smoke + phealth) * (age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt)
formula_X_lsem <- phealth ~ age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt

## Do Sim ----

set.seed(1)
seeds <- sample.int(10e6, NUM_SEED)

tmp <- do_simulation_rlearn(
  data          = meps,
  rlearner_fits = rlearner_fits,
  i_train       = i_train,
  i_test        = i_test,
  formula_m     = formula_m_bart,
  formula_y     = formula_y_bart,
  formula_ps    = formula_ps_bart,
  formula_y_lm  = formula_y_lsem,
  formula_m_lm  = formula_m_lsem,
  mediator_name = "phealth", 
  outcome_name  = "logY", 
  treat_name    = "smoke", 
  n_iter        = 8000, 
  burnin        = 4000, 
  thin          = 4, 
  seeds         = seeds[shuffle]
)
