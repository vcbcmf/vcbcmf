## Load ----

library(ggplot2)
library(latex2exp)
library(cowplot)
library(tidyverse)
library(tidybayes)

source("lib/tikzprint.R")

out_bart <- readRDS('cache/01_fit_bcmf_out_bart.rds')
idx_chains <- which(out_bart$chain != 2)

## Average direct and indirect effects ----

avg_zeta <- rowMeans(out_bart$zeta_samples[idx_chains,])
avg_delta <- rowMeans(out_bart$tau_samples[idx_chains] * out_bart$d_samples[idx_chains,])
avg_tau <- avg_zeta + avg_delta

big_df <- data.frame(
    zeta = avg_zeta,
    delta = avg_delta,
    tau = avg_tau,
    iteration = 1:length(avg_zeta)) %>% 
  pivot_longer(c(zeta, delta, tau), 
               names_to = "Parameter", 
               values_to = "Value") %>%
  mutate(Parameter = case_when(
    Parameter == "zeta" ~ "$\\bar\\zeta$",
    Parameter == "delta" ~ "$\\bar\\delta$",
    Parameter == "tau" ~ "$\\bar\\tau$",
  ))

## Make plots ----

plot_avg_zeta <- ggplot(data.frame(avg_zeta = avg_zeta), aes(x = avg_zeta)) +
  geom_histogram(bins = 40, color = 'white', fill = 'cadetblue3') +
  xlab("$\\bar\\zeta$") + 
  ylab("Frequency") + 
  theme_bw()

plot_avg_delta <- ggplot(data.frame(avg_delta = avg_delta), aes(x = avg_delta)) +
  geom_histogram(bins = 40, color = 'white', fill = 'coral1') +
  xlab("$\\bar\\delta$") + 
  ylab("") + 
  theme_bw()

plot_all_three <- 
  ggplot(big_df, aes(x = Value, y = Parameter, fill = Parameter)) + 
  geom_halfeyeh() + 
  scale_colour_identity() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4,"Blues")[c(3,2,1)+1]) + 
  theme_classic() + theme(legend.position = "none") + 
  geom_vline(xintercept = 0, linetype = 2) + 
  xlab("")

tikzprint(
  fig = plot_grid(plot_avg_zeta, plot_avg_delta, align = 'v'), 
  file_name = "02_make_fig1",
  folder = "figures",
  height = 3
)

tikzprint(
  fig       = plot_all_three, 
  file_name = "02_make_fig1_v2", 
  folder    = "figures", 
  height    = 4*.8, width = 7*.8
)
