## Load ----

library(tidyverse)
source('lib/get_clever_cov.R')
source('lib/get_ps.R')
source('lib/bart_mediate.R')
source('lib/lsem_mediate.R')
source('lib/tikzprint.R')
meps <- readRDS("data/meps.rds")

## Cache ----

out_bart <- readRDS("cache/01_fit_bcmf_out_bart.rds")

## LSEM Mediate ----

formula_m_lsem <- phealth ~ smoke *
  (age + bmi + edu + log(income + 1000) + 
     povlev + region + sex + marital + race + seatbelt)

formula_y_lsem <- logY ~ (smoke + phealth) * 
  (age + bmi + edu + log(income + 1000) + povlev + region + 
     sex + marital + race + seatbelt)

formula_X_lsem <- 
  phealth ~ age + bmi + edu + log(income + 1000) + povlev + region + 
  sex + marital + race + seatbelt

fit_m <- lm(formula_m_lsem, data = meps)
fit_y <- lm(formula_y_lsem, data = meps)

mediated_lsem <- lsem_mediate(
  fit_m = fit_m, 
  fit_y = fit_y, 
  formula_X = formula_X_lsem, 
  data_test = meps, 
  mediator_name = "phealth", 
  treat_name = "smoke"
)

## Dataframe ----

plot_data_bart <- tibble(
  delta = colMeans(out_bart$d_samples * out_bart$tau_samples), 
  zeta = colMeans(out_bart$zeta_samples),
  Method = "BCMF",
  Race = meps$race,
) %>% mutate(subjid = 1:length(delta))

plot_data_lsem <- tibble(
  delta = mediated_lsem$delta, 
  zeta = mediated_lsem$zeta,
  Method = "LSEM",
  Race = meps$race,
) %>% mutate(subjid = 1:length(delta))

plot_data <- rbind(plot_data_bart, plot_data_lsem)
plot_longer <- plot_data %>% 
  pivot_longer(
    cols = delta:zeta, 
    names_to = "Parameter", 
    values_to = "Value"
) %>% mutate(Parameter = ifelse(Parameter == "delta", "$\\delta$", "$\\zeta$"))


## Plot ----

compare_effect <- ggplot(plot_longer, aes(x = Method, y = Value, fill = Method)) + 
  geom_boxplot(alpha = 0.3) + 
  xlab("Method") + 
  ylab("Individual Effect") + 
  facet_wrap(~Parameter, scales = "free") + 
  theme_bw()

compare_effect_2 <- ggplot(plot_longer, aes(x = Race, y = Value, fill = Method)) + 
  geom_boxplot(alpha = 0.3, outlier.shape = NA) + 
  xlab("Method") + 
  ylab("Individual Effect") + 
  facet_wrap(~Parameter, scales = "free") +
  geom_hline(yintercept = 0, linetype = 2) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, vjust = 1))

print(compare_effect_2)
  
tikzprint(
  fig = compare_effect,
  file_name = "11_compare_effect",
  folder = "figures",
  clean = TRUE,
  height = 3.5
)

tikzprint(
  fig = compare_effect_2,
  file_name = "11_compare_effect_2",
  folder = "figures",
  clean = TRUE,
  height = 3.5
)