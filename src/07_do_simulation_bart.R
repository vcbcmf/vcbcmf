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
source('lib/rlasso_mediate.R')

## Divide MEPS into training and test sets

meps <- readRDS('data/meps.rds')
n <- nrow(meps)
set.seed(12345)
i_train <- sample(1:nrow(meps), floor(n/2))
i_test <- c(1:nrow(meps))[-i_train]

## Compute the true effects under BART model ----

formula_m_bart <- phealth ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
formula_y_bart <- logY ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt + phealth
formula_ps_bart <- smoke ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt

if(!file.exists("cache/07_do_simulation_bart_bart.rds")) {
  out_bart         <- readRDS('cache/01_fit_bcmf_out_bart.rds')
  idx_chains       <- which(out_bart$chain != 2)
  mu_y_hat_bart    <- colMeans(out_bart$mu_y_samples[idx_chains,])
  zeta_hat_bart    <- colMeans(out_bart$zeta_samples[idx_chains,])
  d_hat_bart       <- colMeans(out_bart$d_samples[idx_chains,])
  mu_m_hat_bart    <- colMeans(out_bart$mu_m_samples[idx_chains,])
  tau_hat_bart     <- colMeans(out_bart$tau_samples[idx_chains,])
  sigma_y_hat_bart <- mean(out_bart$sigma_y_samples[idx_chains])
  sigma_m_hat_bart <- mean(out_bart$sigma_m_samples[idx_chains])
  rm(out_bart)
  gc()
  
  bart <- list(mu_y_hat_bart, zeta_hat_bart, d_hat_bart, mu_m_hat_bart, 
               tau_hat_bart, sigma_y_hat_bart, sigma_m_hat_bart)
  saveRDS(bart, file = "cache/07_do_simulation_bart_bart.rds")
}

bart <- readRDS("cache/07_do_simulation_bart_bart.rds")
mu_y_hat_bart <- bart[[1]]; zeta_hat_bart <- bart[[2]];  d_hat_bart <- bart[[3]]
mu_m_hat_bart <- bart[[4]]; tau_hat_bart <- bart[[5]]; 
sigma_y_hat_bart <- bart[[6]];  sigma_m_hat_bart <- bart[[7]]


## Compute the true effects with LSEM ----

formula_m_lsem <- phealth ~ smoke * (age + bmi + edu + log(income + 1000) + log(povlev + 100) + region + sex + marital + race + seatbelt)
formula_y_lsem <- logY ~ (smoke + phealth) * (age + bmi + edu + log(income + 1000) + log(povlev + 100) + region + sex + marital + race + seatbelt)
formula_X_lsem <- phealth ~ age + bmi + edu + log(income + 1000) + log(povlev + 100) + region + sex + marital + race + seatbelt

fit_m_lsem <- lm(formula_m_lsem, data = meps)
fit_y_lsem <- lm(formula_y_lsem, data = meps)
# coefs_m <- readRDS("data/coefs_m.rds")
# fit_m_lsem$coefficients <- coefs_m$beta

out_lsem <- lsem_mediate(fit_m_lsem, fit_y_lsem, formula_X_lsem, meps, 'phealth', 'smoke')
saveRDS(out_lsem, 'data/out_lsem.rds')
zeta_hat_lsem <- out_lsem$zeta
delta_hat_lsem <- out_lsem$delta
sigma_m_hat_lsem <- out_lsem$sigma_m
sigma_y_hat_lsem <- out_lsem$sigma_y
rm(out_lsem)
gc()

## Checking for reasonbleness of stuff ----

# library(rpart)
# library(rpart.plot)
# meps_X <- meps %>% select(-logY, -phealth)
# 
# tau_pos_lsem <- ifelse(out_lsem$tau_hat > 0, "Pos", "Neg")
# delt_pos_lsem <- ifelse(out_lsem$delta > 0, "Pos", "Neg")
# 
# rpart.plot(rpart(delt_pos_lsem ~ ., data = meps_X))
# rpart.plot(rpart(tau_pos_lsem ~ ., data = meps_X))
# 
# hist(out_lsem$delta / out_lsem$tau_hat)
# hist(d_hat_bart)

## Running Simulation ----

set.seed(1)
seeds <- sample.int(10e6, NUM_SEED)
simulation_bart <- do_simulation_bart(meps, i_train, i_test, formula_m_bart,
                                      formula_y_bart, formula_ps_bart,
                                      fit_m_lsem, fit_y_lsem, formula_X_lsem,
                                      'phealth', 'logY', 'smoke',
                                      mu_y_hat_bart, zeta_hat_bart, d_hat_bart,
                                      mu_m_hat_bart, tau_hat_bart,
                                      sigma_y_hat_bart, sigma_m_hat_bart,
                                      zeta_hat_lsem, delta_hat_lsem,
                                      sigma_y_hat_lsem, sigma_m_hat_lsem,
                                      8000, 4000, 4, seeds[shuffle])
