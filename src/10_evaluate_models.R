library(SoftBart)
library(tidyverse)
library(Metrics)
library(progress)

source('lib/get_clever_cov.R')
source('lib/get_ps.R')
source('lib/bart_mediate.R')
source('lib/rlasso_mediate.R')

# MEPS train and test set
meps <- readRDS('data/meps.rds')
n <- nrow(meps)
set.seed(12345)
i_train <- sample(1:nrow(meps), floor(n/2))
i_test <- c(1:nrow(meps))[-i_train]
meps_train <- meps[i_train,]
meps_test <- meps[i_test,]

# BART Fit
formula_m_bart <- phealth ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
formula_y_bart <- logY ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt + phealth
formula_ps_bart <- smoke ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt

clever_cov <- get_clever_cov(meps_train, meps, formula_m_bart, 'phealth', 'logY', 'smoke')
m0_hat <- clever_cov$m0_hat
m1_hat <- clever_cov$m1_hat
rm(clever_cov)
gc()

pi_hat <- get_ps(meps_train, meps, formula_ps_bart)

out_bart <- bart_mediate(meps_train, meps_test, formula_m_bart, formula_y_bart,
                         pi_hat[i_train], pi_hat[i_test],
                         m0_hat[i_train], m0_hat[i_test],
                         m1_hat[i_train], m1_hat[i_test],
                         'phealth', 'logY', 'smoke',
                         8000, 4000)

mu_y_hat_bart    <- colMeans(out_bart$mu_y_samples)
zeta_hat_bart    <- colMeans(out_bart$zeta_samples)
d_hat_bart       <- colMeans(out_bart$d_samples)
mu_m_hat_bart    <- colMeans(out_bart$mu_m_samples)
tau_hat_bart     <- colMeans(out_bart$tau_samples)
sigma_y_hat_bart <- mean(out_bart$sigma_y_samples)
sigma_m_hat_bart <- mean(out_bart$sigma_m_samples)
rm(out_bart)
gc()

A <- meps$smoke
m <- meps$phealth
m_pred_bart <- mu_m_hat_bart + A[i_test] * tau_hat_bart
y_pred_bart <- mu_y_hat_bart + A[i_test] * zeta_hat_bart + m[i_test] * d_hat_bart

# LSEM Fit
formula_m_lsem <- phealth ~ smoke * (age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt)
formula_y_lsem <- logY ~ (smoke + phealth) * (age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt)

fit_m_lsem <- lm(formula_m_lsem, data = meps_train)
fit_y_lsem <- lm(formula_y_lsem, data = meps_train)

m_pred_lsem <- predict(fit_m_lsem, newdata = meps_test)
y_pred_lsem <- predict(fit_y_lsem, newdata = meps_test)

# Heldout Correlation
cor(meps_test$phealth, m_pred_bart)
cor(meps_test$phealth, m_pred_lsem)
cor(meps_test$logY, y_pred_bart)
cor(meps_test$logY, y_pred_lsem)

# Wilcoxon Test
bart_sq_err <- (meps_test$logY - y_pred_bart)^2
lsem_sq_err <- (meps_test$logY - y_pred_lsem)^2
wilcox.test(bart_sq_err, lsem_sq_err, paired = TRUE)

bart_m_err <- (meps_test$phealth - m_pred_bart)^2
lsem_m_err <- (meps_test$phealth - m_pred_lsem)^2
wilcox.test(bart_m_err, lsem_m_err, paired = TRUE)

# Linear tests

library(broom)
library(estimatr)
lm_robust(logY ~ y_pred_lsem + y_pred_bart, data = meps_test) %>% 
  tidy() %>% 
  knitr::kable(digits = 4, booktabs = TRUE, format = 'latex')
lm_robust(phealth ~ m_pred_lsem + m_pred_bart, data = meps_test) %>% 
  tidy() %>% 
  knitr::kable(digits = 4, booktabs = TRUE, format = "latex")

