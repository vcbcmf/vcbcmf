## Load ----

library(dplyr)
library(tidyverse)
library(rpart)
library(rpart.plot)
library(tikzDevice)
library(tidybayes)

meps <- readRDS('data/meps.rds')
source("lib/tikzprint.R")

## Load Results ----

out_bart <- readRDS('cache/01_fit_bcmf_out_bart.rds')
idx_chains <- which(out_bart$chain != 2)

indv_zeta <- out_bart$zeta_samples[idx_chains,]
indv_delta <- out_bart$tau_samples[idx_chains,] * out_bart$d_samples[idx_chains,]

meps_post <- meps %>%
  mutate(zeta_hat = colMeans(indv_zeta),
         delta_hat = colMeans(indv_delta)) %>%
  select(age, bmi, edu, income, povlev, region, sex, marital, race, seatbelt, delta_hat, zeta_hat)

subgroup1 <- meps_post %>% filter(age < 34 & race == 'White')
subgroup2 <- meps_post %>% filter(age < 34 & race != 'White')
subgroup3 <- meps_post %>% filter(age >= 67)
subgroup4 <- meps_post %>% filter(age >= 34 & age < 67 & race == 'White')
subgroup5 <- meps_post %>% filter(age >= 34 & age < 67 & race != 'White')

# subgroup_labels <- c('white, age < 34', 'non-white, age < 34', 'age ≥ 67',
#                      'white, 34 ≤ age < 67', 'non-white, 34 ≤ age < 67')
subgroup_labels <- c('white, age $<$ 34', 'non-white, age $<$ 34', 'age $\\ge$ 67',
                     'white, 34 $\\le$ age $<$ 67', 'non-white, 34 $\\le$ age $<$ 67')
df_subgroups = data.frame(
  subgroup = c(rep(subgroup_labels[1], nrow(subgroup1)),
               rep(subgroup_labels[2], nrow(subgroup2)),
               rep(subgroup_labels[3], nrow(subgroup3)),
               rep(subgroup_labels[4], nrow(subgroup4)),
               rep(subgroup_labels[5], nrow(subgroup5))),
  value = c(subgroup1$delta_hat, subgroup2$delta_hat, subgroup3$delta_hat,
            subgroup4$delta_hat, subgroup5$delta_hat)
)

# subgroup_order <- c('white, 34 ≤ age < 67', 'age ≥ 67', 'white, age < 34',
#                     'non-white, 34 ≤ age < 67', 'non-white, age < 34')
subgroup_order <- c('white, 34 $\\le$ age $<$ 67', 'age $\\ge$ 67', 'white, age $<$ 34',
                    'non-white, 34 $\\le$ age $<$ 67', 'non-white, age $<$ 34')
subgroup_plot <- df_subgroups %>%
  ggplot(aes(x = value, 
             y = factor(subgroup, level = subgroup_order), 
             fill = subgroup)) +
  stat_halfeye(alpha = 1, show.legend = FALSE) + 
  xlab('$\\delta(G_i)$') +
  ylab("") + 
  scale_fill_manual(
    values = c("#3182BD", "#9ECAE1", "#C6DBEF", "#08519C", "#6BAED6")
  ) +
  theme_bw() + 
  # theme(text = element_text(family = "Times New Roman")) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))
  )

tikzprint(
  fig = subgroup_plot,
  file_name = "04_make_fig3",
  folder = "figures",
  height = 5
)