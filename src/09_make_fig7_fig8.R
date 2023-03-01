## Load ----

library(tidyverse)
library(ggplot2)
library(latex2exp)
library(cowplot)
library(Metrics)
library(extrafont)
# font_import()

source("lib/load_sim.R")
source("lib/lsem_mediate.R")
source("lib/tikzprint.R")

## MEPS data ----

meps <- readRDS('data/meps.rds')
n <- nrow(meps)
set.seed(12345)
i_train <- sample(1:nrow(meps), floor(n/2))
i_test <- c(1:nrow(meps))[-i_train]

subgroup_labels <- c('white, age < 34', 'non-white, age < 34', 'age ≥ 67',
                     'white, 34 ≤ age < 67', 'non-white, 34 ≤ age < 67')
subgroup1 <- which(meps$age < 34 & meps$race == 'White')
subgroup2 <- which(meps$age < 34 & meps$race != 'White')
subgroup3 <- which(meps$age >= 67)
subgroup4 <- which(meps$age >= 34 & meps$age < 67 & meps$race == 'White')
subgroup5 <- which(meps$age >= 34 & meps$age < 67 & meps$race != 'White')

## Simulation Results ----

simulation_bart_bart <- load_sim(folder = "Simulation BARTBART", setting = "BART", method = "BART")
simulation_lsem_bart <- load_sim("Simulation LSEMBART", "LSEM", "BART")
simulation_bart_lsem <- load_sim("Simulation BARTLSEM", "BART", "LSEM")
simulation_lsem_lsem <- load_sim("Simulation LSEMLSEM", "LSEM", "LSEM")
simulation_rlearn_bart <- load_sim("Simulation RLEARN", "RLEARN", "BART")
simulation_rlearn_lsem <- load_sim("Simulation RLEARN LSEM", "RLEARN", "LSEM")

## Get True Values ----

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

out_bart <- readRDS('cache/01_fit_bcmf_out_bart.rds')
zeta_hat_bart    <- colMeans(out_bart$zeta_samples)
delta_hat_bart   <- colMeans(out_bart$d_samples) * colMeans(out_bart$tau_samples)
rm(out_bart)
gc()

## Make Truth/Fit data frames ----

# bb_avg <- simulation_bart_bart$avg_results %>%
#   mutate(group = 'average', .after = seed) %>%
#   mutate(zeta_true = rep(mean(zeta_hat_bart), 200),
#          delta_true = rep(mean(delta_hat_bart), 200))

bb_avg <- simulation_bart_bart$avg_results %>% mutate(group = 'average', .after = seed)
bb_indv <- simulation_bart_bart$indv_results
bb_subgroups <- simulation_bart_bart$subgroup_results
bb_tree_subgroups <- simulation_bart_bart$tree_results
lb_avg <- simulation_lsem_bart$avg_results %>% mutate(group = 'average', .after = seed)
lb_indv <- simulation_lsem_bart$indv_results
ll_avg <- simulation_lsem_lsem$avg_results %>% mutate(group = 'average', .after = seed)
ll_indv <- simulation_lsem_lsem$indv_results
bl_avg <- simulation_bart_lsem$avg_results %>% mutate(group = 'average', .after = seed)
bl_indv <- simulation_bart_lsem$indv_results
rb_avg <- simulation_rlearn_bart$avg_results %>% mutate(group = 'average', .after = seed)
rb_indv <- simulation_rlearn_bart$indv_results
rl_avg <- simulation_rlearn_lsem$avg_results %>% mutate(group = 'average', .after = seed)
rl_indv <- simulation_rlearn_lsem$indv_results

indv_df       <- rbind(bb_indv, lb_indv, ll_indv, bl_indv, rb_indv, rl_indv)
avg_subgroups_df <- rbind(bb_avg, bb_subgroups, bb_tree_subgroups, lb_avg, ll_avg, bl_avg, rb_avg, rl_avg)

## Coverage, RMSE, bias, interval length for average zeta and delta ----

avg_df_zeta <- avg_subgroups_df %>%
  filter(group == 'average') %>%
  group_by(Setting, Method) %>%
  summarize(cov_zeta = mean(zeta_catch),
            rmse_zeta = rmse(zeta_true, zeta_mean),
            bias_zeta = abs(bias(zeta_true, zeta_mean)),
            len_zeta = mean(zeta_len)) %>%
  kableExtra::kable(format = "latex", booktabs = TRUE, digits = 2)

avg_df_delta <- avg_subgroups_df %>%
  filter(group == 'average') %>%
  group_by(Setting, Method) %>%
  summarize(cov_delta = mean(delta_catch),
            rmse_delta = rmse(delta_true, delta_mean),
            bias_delta = abs(bias(delta_true, delta_mean)),
            len_delta = mean(delta_len)) %>%
  kableExtra::kable(format = "latex", booktabs = TRUE, digits = 2)

## Coverage, RMSE, bias, interval length for subgroup zeta and delta ----

sub_df_zeta <- avg_subgroups_df %>%
  filter(group %in% subgroup_labels) %>%
  group_by(Setting, Method, group) %>%
  summarize(cov_zeta = mean(zeta_catch),
            rmse_zeta = rmse(zeta_true, zeta_mean),
            bias_zeta = abs(bias(zeta_true, zeta_mean)),
            len_zeta = mean(zeta_len)) %>%
  kableExtra::kable(format = "latex", booktabs = TRUE, digits = 2)

sub_df_delta <- avg_subgroups_df %>%
  filter(group %in% subgroup_labels) %>%
  group_by(Setting, Method, group) %>%
  summarize(cov_delta = mean(delta_catch),
            rmse_delta = rmse(delta_true, delta_mean),
            bias_delta = abs(bias(delta_true, delta_mean)),
            len_delta = mean(delta_len)) %>%
  kableExtra::kable(format = "latex", booktabs = TRUE, digits = 2)



## Coverage, RMSE, bias, interval length for individual zeta & delta ----

indv_df_zeta <- indv_df %>%
  group_by(Setting, Method, subj_id) %>%
  summarize(cov_zeta = mean(zeta_catch),
            rmse_zeta = rmse(zeta_true, zeta_mean),
            zeta_true = mean(zeta_true),
            bias_zeta = abs(bias(zeta_true, zeta_mean)),
            len_zeta = mean(zeta_len))

indv_df_zeta %>% 
  summarise_all(mean) %>% 
  kableExtra::kable(x = ., format = "latex", booktabs = TRUE, digits = 2)

indv_df_delta <- indv_df %>%
  group_by(Setting, Method, subj_id) %>%
  summarize(cov_delta = mean(delta_catch),
            rmse_delta = rmse(delta_true, delta_mean),
            bias_delta = abs(bias(delta_true, delta_mean)),
            delta_true = mean(delta_true),
            len_delta = mean(delta_len))

indv_df_delta %>% 
  summarise_all(mean) %>% 
  kableExtra::kable(x = ., format = "latex", booktabs = TRUE, digits = 2)

# # Coverage, RMSE, bias, interval length for subgroup zeta & delta
# avg_subgroups_df %>%
#   filter(group %in% subgroup_labels) %>%
#   group_by(Setting, Method, group) %>%
#   summarize(cov_zeta = mean(zeta_catch),
#             rmse_zeta = rmse(zeta_true, zeta_mean),
#             bias_zeta = abs(bias(zeta_true, zeta_mean)),
#             len_zeta = mean(zeta_len))
# 
# avg_subgroups_df %>%
#   filter(group %in% subgroup_labels) %>%
#   group_by(Setting, Method, group) %>%
#   summarize(cov_delta = mean(delta_catch),
#             rmse_delta = rmse(delta_true, delta_mean),
#             bias_delta = abs(bias(delta_true, delta_mean)),
#             len_delta = mean(delta_len))
# 
# avg_subgroups_df %>%
#   group_by(Setting, Method) %>%
#   filter(!(group %in% subgroup_labels) & group != 'average') %>%
#   summarize(cov_zeta = mean(zeta_catch),
#             len_zeta = mean(zeta_len))
# 
# avg_subgroups_df %>%
#   group_by(Setting, Method) %>%
#   filter(!(group %in% subgroup_labels) & group != 'average') %>%
#   summarize(cov_delta = mean(delta_catch),
#             len_delta = mean(delta_len))
# 

## Zeta ----

f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  # r['middle'] <- mean(x)
  r
}

plot_cov_zeta <- 
  ggplot(indv_df_zeta, aes(
    x = Setting, 
    y = cov_zeta,
    fill = Method)) +
  stat_summary(alpha = 0.5, fun.data = f, geom = 'boxplot', outlier.shape = NA, position = 'dodge') + 
  ylab('Coverage') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw() + 
  geom_hline(yintercept = 0.95, lty = 2, color = 'darkgray')
legend        <- get_legend(plot_cov_zeta)
plot_cov_zeta <- plot_cov_zeta + theme(legend.position = 'none')
print(plot_cov_zeta)

# plot_cov_zeta <- 
#   ggplot(indv_df_zeta, aes(x = Setting, y = cov_zeta, fill = Method)) +
#   geom_boxplot(aes(
#     ymin = quantile(cov_zeta, 0.025), 
#     ymax = quantile(cov_zeta, 0.975)), 
#     outlier.shape = NA, alpha = 0.5) +
#   ylab('Coverage') + xlab('Truth') + labs(fill = 'Fit') +
#   theme_bw() +
#   geom_hline(yintercept = 0.95, lty = 2, color = 'darkgray')
# legend <- get_legend(plot_cov_zeta)
# plot_cov_zeta <- plot_cov_zeta + theme(legend.position = 'none')
# print(plot_cov_zeta)

plot_rmse_zeta <- 
  ggplot(indv_df_zeta, aes(x = Setting, y = rmse_zeta, fill = Method)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('RMSE') + xlab('Truth') + 
  labs(fill = 'Fit') +
  theme_bw()

# print(plot_rmse_zeta)


plot_bias_zeta <- 
  ggplot(indv_df_zeta, aes(x = Setting, y = bias_zeta, fill = Method)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('Bias') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw()

plot_len_zeta <- 
  ggplot(indv_df_zeta, aes(x = Setting, y = len_zeta, fill = Method)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('Interval Length') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw()

grid_zeta <- plot_grid(plot_cov_zeta, plot_rmse_zeta, plot_bias_zeta, plot_len_zeta)
plots_zeta <- plot_grid(grid_zeta, legend, ncol = 2, rel_widths = c(1, .1))
pz <- ggdraw(add_sub(plots_zeta, '$\\zeta(x)$', size = 12,
               vpadding = grid::unit(0, "lines"), y = 5, x = 0.48, vjust = 4.5))

tikzprint(
  fig = pz,
  file_name = "09_plot_zeta",
  folder = "figures",
  clean = TRUE,
  height = 5,
  width = 8
)

## Delta ----

f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  # r['middle'] <- mean(x)
  r
}

plot_cov_delta <- 
  ggplot(indv_df_delta, aes(
    x = Setting, 
    y = cov_delta,
    fill = Method)) +
  stat_summary(alpha = 0.5, fun.data = f, geom = 'boxplot', outlier.shape = NA, position = 'dodge') + 
  ylab('Coverage') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw() + 
  geom_hline(yintercept = 0.95, lty = 2, color = 'darkgray')

legend        <- get_legend(plot_cov_delta)
plot_cov_delta <- plot_cov_delta + theme(legend.position = 'none')

# plot_cov_delta <- 
#   ggplot(indv_df_delta, aes(
#     x = Setting, 
#     y = cov_delta,
#     ymin = quantile(cov_delta, 0.1),
#     fill = Method)) +
#   geom_boxplot(alpha = 0.5, show.legend = FALSE, outlier.shape = NA) +
#   ylab('Coverage') + xlab('Truth') + labs(fill = 'Fit') +
#   theme_bw() +
#   geom_hline(yintercept = 0.95, lty = 2, color = 'darkgray')

print(plot_cov_delta)

plot_rmse_delta <- ggplot(indv_df_delta, aes(x = Setting, y = rmse_delta, fill = Method)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('RMSE') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw()

plot_bias_delta <- ggplot(indv_df_delta, aes(x = Setting, y = bias_delta, fill = Method)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('Bias') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw()

plot_len_delta <- ggplot(indv_df_delta, aes(x = Setting, y = len_delta, fill = Method)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('Interval Length') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw()

grid_delta <- plot_grid(plot_cov_delta, plot_rmse_delta, plot_bias_delta, plot_len_delta)
plots_delta <- plot_grid(grid_delta, legend, ncol = 2, rel_widths = c(1, .1))
pd <-
  ggdraw(
    add_sub(
      plots_delta,
      '$\\delta(x)$',
      size = 12,
      vpadding = grid::unit(0, "lines"),
      y = 5,
      x = 0.48,
      vjust = 4.5
    )
  )

tikzprint(
  fig = pd,
  file_name = "09_plot_delta",
  folder = "figures",
  clean = TRUE,
  height = 5,
  width = 8
)
