## Load ----

source("src/09_make_fig7_fig8.R")

## Make Figure ----

p_delta <- 
  ggplot(
    data = indv_df_delta, 
    aes(x = abs(delta_true - mean(delta_true)), y = cov_delta)) + 
  # geom_point() + 
  # geom_density2d_filled() + 
  geom_hex(bins = 10) + 
  stat_binhex(aes(label=round(..count../sum(..count..)*6, 2)), geom="text", bins=10, colour="white", size = 3) +
  facet_grid(Setting~Method) + 
  theme_bw() + 
  # scale_fill_viridis_c() +
  xlab("$|\\delta(X_i) - \\bar \\delta|$") + 
  ylab("Coverage") + 
  geom_hline(yintercept = 0.95, lty = 2) + 
  theme(legend.position = "none")

p_zeta <- 
  ggplot(
    data = indv_df_zeta, 
    aes(x = abs(zeta_true - mean(zeta_true)), y = cov_zeta)) + 
  # geom_point() + 
  # geom_density2d_filled() + 
  geom_hex(bins = 10) + 
  stat_binhex(aes(label=round(..count../sum(..count..)*6, 2)), geom="text", bins=10, colour="white", size = 3) +
  facet_grid(Setting~Method) + 
  theme_bw() + 
  # scale_fill_viridis_c() +
  xlab("$|\\zeta(X_i) - \\bar \\zeta|$") + 
  ylab("Coverage") + 
  geom_hline(yintercept = 0.95, lty = 2) + 
  theme(legend.position = "none")

p_both <- gridExtra::grid.arrange(p_delta, p_zeta, nrow = 2)

tikzprint(p_both, "12_coverage_by_size", "figures", height = 12, width = 10)
tikzprint(p_zeta, "12_coverage_by_size_zeta", "figures")
tikzprint(p_delta, "12_coverage_by_size_delta", "figures")
