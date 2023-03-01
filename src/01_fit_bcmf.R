library(SoftBart)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(progress)

source('lib/bart_mediate2.R')
source('lib/get_clever_cov.R')
source('lib/get_ps.R')
source('lib/tikzprint.R')
source('lib/rlasso_mediate.R')

## Load MEPS data ----
meps <- readRDS('data/meps.rds')

## Define formulas ---- 
formula_m <- phealth ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
formula_y <- logY ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt + phealth
formula_ps <- smoke ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt

## Get clever covariates and propensity score ----
if(!file.exists("cache/01_fit_bcmf_clever_cov.rds")) {
  set.seed(123)
  clever_cov <- get_clever_cov(meps, meps, formula_m,
                               'phealth', 'logY', 'smoke')
  saveRDS(clever_cov, file = "cache/01_fit_bcmf_clever_cov.rds")
} 

clever_cov <- readRDS("cache/01_fit_bcmf_clever_cov.rds")
m0_hat <- clever_cov$m0_hat
m1_hat <- clever_cov$m1_hat

if(!file.exists("cache/01_fit_bcmf_pi_hat.rds")) {
  set.seed(456)
  pi_hat <- get_ps(meps, meps, formula_ps)
  saveRDS(pi_hat, file = "cache/01_fit_bcmf_pi_hat.rds")
}

pi_hat <- readRDS("cache/01_fit_bcmf_pi_hat.rds")

## Fit BCMF ----
if(!file.exists("cache/01_fit_bcmf_out_bart.rds")) {
  set.seed(789)
  out_bart <- bart_mediate2(
    data_train = meps, 
    data_test = meps, 
    formula_m = formula_m, 
    formula_y = formula_y, 
    pi_hat_train = pi_hat, 
    pi_hat_test = pi_hat,
    m0_hat_train = m0_hat, 
    m0_hat_test = m0_hat, 
    m1_hat_train = m1_hat, 
    m1_hat_test = m1_hat,
    mediator_name = 'phealth', 
    outcome_name = 'logY', 
    treat_name = 'smoke',
    n_iter = 5000, 
    burnin = 4000, 
    chain = 4
  )
  
  saveRDS(out_bart, 'cache/01_fit_bcmf_out_bart.rds')
}

out_bart <- readRDS("cache/01_fit_bcmf_out_bart.rds")

## Chain Stacking ----

niter <- 4000
Y <- meps$logY
M <- meps$phealth
A <- meps$smoke
Ymat <- t(matrix(Y, nrow = length(Y), ncol = niter))
Mmat <- t(matrix(M, nrow = length(M), ncol = niter))

mu_y <- out_bart %>% with(mu_y_samples + zeta_samples * A + d_samples * M)
sd_y <- matrix(out_bart$sigma_y_samples, nrow = nrow(mu_y), ncol = ncol(mu_y))

mu_m <- out_bart %>% with(mu_m_samples + tau_samples * A)
sd_m <- matrix(out_bart$sigma_m_samples, nrow = nrow(mu_y), ncol = ncol(mu_y))

loglik_y <- dnorm(Ymat, mu_y, sd_y, log = TRUE)
loglik_m <- dnorm(Mmat, mu_m, sd_m, log = TRUE)

loglik_agg <- loglik_y + loglik_m

loo_chain <- function(i) {
  loo::loo(loglik_agg[1:1000 + (i - 1) * 1000,])
}

moo <- lapply(1:4, loo_chain)
ham <- exp(sapply(1:4, function(i) exp(moo[[i]]$pointwise[,1])))

library(rstan)
stan_code <- "

data {
  int<lower=1> N;
  int<lower = 1> K;
  matrix[N,K] elpd;
  vector[K] lambda;
}

parameters {
  simplex[K] w;
}

model {
  w ~ dirichlet(lambda);
  target += sum(log(elpd * w));
}
"

sm <- stan_model(model_code = stan_code)

smsamps <- sampling(
  object = sm, 
  data = list(N = 16113, K = 4, elpd = ham, lambda = rep(1/4, 4)), 
  seed = 101112
)

chain_weights <- (smsamps %>% summary() %>% pluck("summary"))[1:4,1]
saveRDS(chain_weights, file = "01_fit_bcmf_chain_weights.rds")
traceplot(smsamps)

## Trace Plots: Chain 2 removed on basis of stacking ----

chain <- out_bart$chain
indv_zeta <- out_bart$zeta_samples
indv_delta <- out_bart$tau_samples * out_bart$d_samples
avg_zeta <- data.frame(i = rep(1:1000, times = 4), avg_zeta = rowMeans(indv_zeta), chain = chain)
avg_delta <- data.frame(i = rep(1:1000, times = 4), avg_delta = rowMeans(indv_delta), chain = chain)
sigma_m <- data.frame(i = rep(1:1000, times = 4), sigma_m = out_bart$sigma_m_samples, chain = chain)
sigma_y <- data.frame(i = rep(1:1000, times = 4), sigma_y = out_bart$sigma_y_samples, chain = chain)

traceplot_avg_zeta <- ggplot(avg_zeta %>% filter(chain != 2), aes(x = i, y = avg_zeta)) +
  geom_line(aes(color = chain), alpha = 0.7) +
  xlab('Iteration') + 
  ylab('$\\bar\\zeta$') + 
  # facet_wrap(~'zeta_avg') + 
  theme_bw()

legend <- get_legend(traceplot_avg_zeta)
traceplot_avg_zeta <- traceplot_avg_zeta + theme(legend.position = 'none')
tikzprint(
  traceplot_avg_zeta, 
  "01_fit_bcmf_zeta_traceplot", 
  "figures", 
  width = 7, 
  height = 4
)

traceplot_avg_delta <- 
  ggplot(avg_delta %>% filter(chain != 2), aes(x = i, y = avg_delta)) +
  geom_line(aes(color = chain), alpha = 0.7) +
  theme_bw() + 
  xlab('Iteration') + 
  ylab('$\\bar\\delta$') + 
  theme_bw() + 
  theme(legend.position = "none")

tikzprint(
  fig = traceplot_avg_delta,
  file_name = "01_fit_bcmf_delta_traceplot",
  folder = "figures",
  width = 7,
  height = 4
)

traceplot_sigma_m <- ggplot(sigma_m %>% filter(chain != 2), aes(x = i, y = sigma_m)) + 
  geom_line(aes(color = chain), alpha = 0.7) +
  theme_bw() + 
  xlab('Iteration') + 
  ylab('$\\sigma_m$') + 
  theme_bw() + 
  theme(legend.position = "none")

traceplot_sigma_y <- ggplot(sigma_y %>% filter(chain != 2), aes(x = i, y = sigma_y)) + 
  geom_line(aes(color = chain), alpha = 0.7) +
  theme_bw() + 
  xlab('Iteration') + 
  ylab('$\\sigma_y$') + 
  theme_bw() + 
  theme(legend.position = "none")

grid <- plot_grid(
  traceplot_avg_zeta, 
  traceplot_avg_delta, 
  traceplot_sigma_m, 
  traceplot_sigma_y, 
  align = 'vh'
)


plotted_grid <- plot_grid(grid, ncol = 2, rel_widths = c(1, .1))
tikzprint(
  fig = plotted_grid, 
  file_name = "01_fit_bcmf_plotted_grid", 
  folder = "figures",
  height = 4, 
  width = 7
)

## Individual Level Traceplots ----

set.seed(131415)
i <- sample.int(nrow(meps), 100)

get_df <- function(j) {
  z <- data.frame(subj_id = as.factor(rep(i, each = 1000))) %>% 
    mutate(iter = rep(1:1000, length(i)), chain = j) %>%
    mutate(zeta = as.vector(indv_zeta[1:1000 + (j - 1) * 1000, i]))
}

get_dfd <- function(j) {
  z <- data.frame(subj_id = as.factor(rep(i, each = 1000))) %>% 
    mutate(iter = rep(1:1000, length(i)), chain = j) %>%
    mutate(delta = as.vector(indv_delta[1:1000 + (j - 1) * 1000, i]))
}

zeta_subset <- map_df(1:4, get_df)

ess_zeta <- zeta_subset %>% filter(chain != 2) %>% 
  group_by(subj_id, chain) %>% 
  summarise(ess = coda::effectiveSize(zeta)) %>% 
  summarise(ess = sum(ess)) %>% 
  mutate(Parameter = "$\\zeta(x)$")


delta_subset <- map_df(1:4, get_dfd)

ess_delta <- delta_subset %>% filter(chain != 2) %>% 
  group_by(subj_id, chain) %>% 
  summarise(ess = coda::effectiveSize(delta)) %>% 
  summarise(ess = sum(ess)) %>%
  mutate(Parameter = "$\\delta(x)$")

ess_df <- rbind(ess_zeta, ess_delta)

ess_plot <- ggplot(ess_df, aes(y = ess, x = Parameter, fill = Parameter)) + 
  geom_boxplot(alpha = 0.7) + 
  ylim(0, 1700) + 
  xlab("") + 
  ylab("ESS") +
  theme_bw() + 
  theme(legend.position = "none")

tikzprint(
  ess_plot, 
  file_name = "01_fit_bcmf_ess_plot", 
  "figures", 
  height = 3, 
  width = 7
)

