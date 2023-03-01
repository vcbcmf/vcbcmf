## Load ----

library(ggplot2)
library(rpart)
library(mgcv)
library(cowplot)
library(latex2exp)

source('lib/projection_tree.R')
source('lib/projection_gam.R')
source('lib/r_sq.R')
source('lib/tikzprint.R')

## Load out ----

out_bart <- readRDS('cache/01_fit_bcmf_out_bart.rds')
meps     <- readRDS('data/meps.rds')

idx_chain <- which(out_bart$chain != 2)

indv_indirect <- out_bart$d_samples[idx_chain,] * out_bart$tau_samples[idx_chain,]
formula_y <- logY ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt + phealth

## Fit projections if needed ----

if(!file.exists("cache/06_make_fig6.rds")) {
  proj_tree_delta <- projection_tree(meps, formula_y, indv_indirect)
  r_sq_tree <- r_sq(indv_indirect, proj_tree_delta)
  proj_gam_delta <- projection_gam(meps, indv_indirect)
  r_sq_gam <- r_sq(indv_indirect, proj_gam_delta)
  make_fig6 <- list(proj_tree_delta, r_sq_tree, proj_gam_delta, r_sq_gam)
  saveRDS(make_fig6, "cache/06_make_fig6.rds")
}

make_fig6 <- readRDS("cache/06_make_fig6.rds")
proj_tree_delta <- make_fig6[[1]]
r_sq_tree <- make_fig6[[2]]
proj_gam_delta <- make_fig6[[3]]
r_sq_gam <- make_fig6[[4]]

rsqs <- r_sq_gam_tree(data = meps, formula_y = formula_y, samples = indv_indirect)

## Make Figures ----

plot_r2_tree <- ggplot(data.frame(r_sq = r_sq_tree$r2), aes(x = r_sq)) +
  geom_histogram(bins = 30, color = 'white', fill = 'aquamarine3') + 
  xlab("") + 
  ylab("Frequency") + 
  facet_wrap(~'Tree') +
  geom_vline(xintercept = rsqs$r_sq_tree) + 
  theme_bw() +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

plot_r2_gam <- ggplot(data.frame(r_sq = r_sq_gam$r2), aes(x = r_sq)) +
  geom_histogram(bins = 30, color = 'white', fill = 'aquamarine3') +
  xlab("") + 
  ylab("") + 
  facet_wrap(~'GAM') +
  geom_vline(xintercept = rsqs$r_sq_gam) + 
  theme_bw()

plot_r2 <- plot_grid(plot_r2_tree, plot_r2_gam, align = 'v')

plot_fig <- ggdraw(plot = add_sub(
  plot_r2, 
  ("$R^2$"), 
  # fontfamily = "Times New Roman", 
  size = 11,
  vpadding = grid::unit(0, "lines"), 
  y = 6, 
  x = 0.54, 
  vjust = 4.5
))

tikzprint(plot_fig, "06_make_fig6", "figures", height = 3.5)