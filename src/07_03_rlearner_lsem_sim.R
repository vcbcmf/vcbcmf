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
source("lib/do_simulation_rlearn_lsem.R")
source('lib/lsem_boot.R')

set.seed(12345)
meps          <- readRDS('data/meps.rds')
rlearner_fits <- readRDS("cache/07_01_rlearn_fits.rds")
n             <- nrow(meps)
i_train       <- sample(1:nrow(meps), floor(n/2))
i_test        <- c(1:nrow(meps))[-i_train]

## A Bunch of Formulas

formula_m_lsem <- phealth ~ smoke * (age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt)
formula_y_lsem <- logY ~ (smoke + phealth) * (age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt)
formula_X_lsem <- phealth ~ age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt

## Do Sim ----

set.seed(1)
seeds <- sample.int(10e6, NUM_SEED)

tmp <- do_simulation_rlearn_lsem(
  data           = meps,
  rlearner_fits  = rlearner_fits,
  i_train        = i_train,
  i_test         = i_test,
  formula_m_lsem = formula_m_lsem,
  formula_y_lsem = formula_y_lsem,
  formula_X_lsem = formula_X_lsem,
  mediator_name  = "phealth",
  outcome_name   = "logY",
  treat_name     = "smoke",
  seeds          = seeds[shuffle]
)
