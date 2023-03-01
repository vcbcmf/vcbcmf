## Load ----

library(dplyr)
library(ggplot2)
library(cowplot)
library(ggdist)
library(mgcv)
source("lib/tikzprint.R")

meps <- readRDS('data/meps.rds')

## Load fit ----

out_bart <- readRDS('cache/01_fit_bcmf_out_bart.rds')
idx_chains <- which(out_bart$chain != 2)

## Individual direct and indirect effects ----

indv_zeta <- out_bart$zeta_samples[idx_chains,]
indv_delta <- out_bart$tau_samples[idx_chains,] * out_bart$d_samples[idx_chains,]

meps_post <- meps %>%
  mutate(zeta_hat = colMeans(indv_zeta),
         delta_hat = colMeans(indv_delta)) %>%
  select(age, bmi, edu, income, povlev, region, sex, marital, race, seatbelt, delta_hat, zeta_hat)

# Posterior projection with GAM
smooth_age <- list()
smooth_bmi <- list()
smooth_edu <- list()
smooth_income <- list()
smooth_povlev <- list()
region_list <- list()
sex_list <- list()
marital_list <- list()
race_list <- list()
seatbelt_list <- list()

if(!file.exists("cache/05_make_fig4_fig5.rds")) {
  pdf(NULL)
  for(i in 1:nrow(indv_zeta)) {
    print(i)
    project_gam <- gam(delta_hat ~ s(age, k = 10) + s(bmi, k = 10) + s(edu, k = 10) +
                         s(log(income + 1000), k = 10) + s(povlev, k = 10) + region + sex + 
                         marital + race + seatbelt,
                       data = meps_post %>% mutate(delta_hat = indv_delta[i,]))
    
    # Constraint for binary/categorical variables to sum to 0
    coefs <- coef(project_gam)
    coefs_region0 <- c(regionMidwest = 0, coefs[grep('region', names(coefs))])
    coefs_region <- coefs_region0 - mean(coefs_region0)
    region_list[[i]] <- data.frame(coefs = coefs_region, 
                                   category = names(coefs_region), iteration = i)
    
    coefs_sex0 <- c(sexFemale = 0, coefs[grep('sex', names(coefs))])
    coefs_sex <- coefs_sex0 - mean(coefs_sex0)
    sex_list[[i]] <- data.frame(coefs = coefs_sex,
                                category = names(coefs_sex), iteration = i)
    
    coefs_marital0 <- c(maritalDivorced = 0, coefs[grep('marital', names(coefs))])
    coefs_marital <- coefs_marital0 - mean(coefs_marital0)
    marital_list[[i]] <- data.frame(coefs = coefs_marital,
                                    category = names(coefs_marital), iteration = i)
    
    coefs_race0 <- c(raceWhite = 0, coefs[grep('race', names(coefs))])
    coefs_race <- coefs_race0 - mean(coefs_race0)
    race_list[[i]] <- data.frame(coefs = coefs_race, 
                                 category = names(coefs_race), iteration = i)
    
    coefs_seatbelt0 <- c(seatbeltAlmostAlways = 0, coefs[grep('seatbelt', names(coefs))])
    coefs_seatbelt <- coefs_seatbelt0 - mean(coefs_seatbelt0)
    seatbelt_list[[i]] <- data.frame(coefs = coefs_seatbelt, 
                                     category = names(coefs_seatbelt), iteration = i)
    
    p <- plot(project_gam)
    smooth_age[[i]] <- data.frame(x = p[[1]]$x, y = p[[1]]$fit, iteration = i)
    smooth_bmi[[i]] <- data.frame(x = p[[2]]$x, y = p[[2]]$fit, iteration = i)
    smooth_edu[[i]] <- data.frame(x = p[[3]]$x, y = p[[3]]$fit, iteration = i)
    smooth_income[[i]] <- data.frame(x = p[[4]]$x, y = p[[4]]$fit, iteration = i)
    smooth_povlev[[i]] <- data.frame(x = p[[5]]$x, y = p[[5]]$fit, iteration = i)
  }
  dev.off()
  
  smooth_age_df <- do.call(rbind, smooth_age)
  smooth_bmi_df <- do.call(rbind, smooth_bmi)
  smooth_edu_df <- do.call(rbind, smooth_edu)
  smooth_income_df <- do.call(rbind, smooth_income)
  smooth_povlev_df <- do.call(rbind, smooth_povlev)
  region_df <- do.call(rbind, region_list)
  sex_df <- do.call(rbind, sex_list)
  marital_df <- do.call(rbind, marital_list)
  race_df <- do.call(rbind, race_list)
  seatbelt_df <- do.call(rbind, seatbelt_list)
  
  make_fig4_fig5 <- list(
    smooth_age_df, 
    smooth_bmi_df, 
    smooth_edu_df, 
    smooth_income_df, 
    smooth_povlev_df, 
    region_df, 
    sex_df, 
    marital_df, 
    race_df, 
    seatbelt_df
  )
  
  saveRDS(make_fig4_fig5, file = "cache/05_make_fig4_fig5.rds")
  
}

make_fig4_fig5 <- readRDS("cache/05_make_fig4_fig5.rds")
smooth_age_df <- make_fig4_fig5[[1]] 
smooth_bmi_df <- make_fig4_fig5[[2]]
smooth_edu_df <- make_fig4_fig5[[3]]
smooth_income_df <- make_fig4_fig5[[4]]
smooth_povlev_df <- make_fig4_fig5[[5]]
region_df <- make_fig4_fig5[[6]]
sex_df <- make_fig4_fig5[[7]]
marital_df <- make_fig4_fig5[[8]]
race_df <- make_fig4_fig5[[9]]
seatbelt_df <- make_fig4_fig5[[10]]

## Fixing race_df bug ----

race_df[1 + 6 * (1:3000 - 1), 2] <- "raceAsian"


## Plot continuous covariates from GAM ----
smooth_age_summary <- smooth_age_df %>%
  group_by(x) %>% summarize(mu = mean(y), LCL = quantile(y, 0.025),
                            UCL = quantile(y, 0.975))
smooth_bmi_summary <- smooth_bmi_df %>%
  group_by(x) %>% summarize(mu = mean(y), LCL = quantile(y, 0.025),
                            UCL = quantile(y, 0.975))
smooth_edu_summary <- smooth_edu_df %>%
  group_by(x) %>% summarize(mu = mean(y), LCL = quantile(y, 0.025),
                            UCL = quantile(y, 0.975))
smooth_income_summary <- smooth_income_df %>%
  group_by(x) %>% summarize(mu = mean(y), LCL = quantile(y, 0.025),
                            UCL = quantile(y, 0.975))
smooth_povlev_summary <- smooth_povlev_df %>%
  group_by(x) %>% summarize(mu = mean(y), LCL = quantile(y, 0.025),
                            UCL = quantile(y, 0.975))

plot_age <- ggplot(smooth_age_summary, aes(x = x, y = mu, ymin = LCL, ymax = UCL)) +
  geom_line(size = 1, lty = 2) + geom_ribbon(alpha = 0.5, fill = 'cornflowerblue') +
  geom_hline(yintercept = 0) + xlab('age') + ylab('s(age)') + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

plot_bmi <- ggplot(smooth_bmi_summary, aes(x = x, y = mu, ymin = LCL, ymax = UCL)) +
  geom_line(size = 1, lty = 2) + geom_ribbon(alpha = 0.5, fill = 'cornflowerblue') +
  geom_hline(yintercept = 0) + xlab('bmi') + ylab('s(bmi)') + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

plot_edu <- ggplot(smooth_edu_summary, aes(x = x, y = mu, ymin = LCL, ymax = UCL)) +
  geom_line(size = 1, lty = 2) + geom_ribbon(alpha = 0.5, fill = 'cornflowerblue') +
  geom_hline(yintercept = 0) + xlab('edu') + ylab('s(edu)') + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

plot_income <- ggplot(smooth_income_summary, aes(x = x, y = mu, ymin = LCL, ymax = UCL)) +
  geom_line(size = 1, lty = 2) + geom_ribbon(alpha = 0.5, fill = 'cornflowerblue') +
  geom_hline(yintercept = 0) + xlab('log(income)') + ylab('s(log(income))') + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

plot_povlev <- ggplot(smooth_povlev_summary, aes(x = x, y = mu, ymin = LCL, ymax = UCL)) +
  geom_line(size = 1, lty = 2) + geom_ribbon(alpha = 0.5, fill = 'cornflowerblue') +
  geom_hline(yintercept = 0) + xlab('poverty level') + ylab('s(poverty level)') + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

gam_plots <- plot_grid(plot_age, plot_bmi, plot_edu, plot_income, align = 'vh')
tikzprint(fig = gam_plots, file_name = "05_make_fig4_fig5", folder = "figures", height = 5)

## Plot binary/categorical covariates from GAM ----

level_order <- c('seatbeltAlways', 'seatbeltNever', 'seatbeltAlmostAlways', 'seatbeltNoCar', 'seatbeltSeldom', 'seatbeltSometimes',
                 'raceWhite', 'raceBlack', 'raceAsian',  'raceIndig', 'raceMulti', 'racePacificIslander',
                 'maritalMarried', 'maritalSeparated', 'maritalDivorced', 'maritalWidowed',
                 'regionNortheast', 'regionMidwest', 'regionWest', 'regionSouth',
                 'sexMale', 'sexFemale')

plot_coefs <- rbind(cbind(sex_df, id = 1), cbind(region_df, id = 2),
                       cbind(marital_df, id = 3), cbind(race_df, id = 4),
                       cbind(seatbelt_df, id = 5))

plot_disc <- 
  ggplot(plot_coefs, aes(
    y = factor(category, level = level_order), 
    x = coefs, 
    fill = as.factor(id))
  ) +
  stat_halfeye(show.legend = FALSE, alpha = 0.6) +
  scale_fill_brewer(palette = 'Set3') +
  xlab('Coefficients') +
  ylab('Category') +
  # xlim(-0.1, 0.25) +
  coord_cartesian(xlim = c(-0.08, 0.08)) +
  theme_classic() +
  scale_colour_identity() +
  geom_hline(
    yintercept = c(6.5, 12.5, 16.6, 20.7), 
    lty = 2, 
    col = 'darkgray', 
    size = 0.4
  )

tikzprint(plot_disc, "05_make_fig4_fig5_plot_disc", "figures", height = 7)
