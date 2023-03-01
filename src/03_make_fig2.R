## Load ----

library(dplyr)
library(rpart)
library(rpart.plot)
library(tikzDevice)

meps <- readRDS('data/meps.rds')

## Load Results ----

out_bart <- readRDS('cache/01_fit_bcmf_out_bart.rds')
idx_chains <- which(out_bart$chain != 2)

## Individual direct and indirect effects ----

indv_zeta <- out_bart$zeta_samples[idx_chains,]
indv_delta <- out_bart$tau_samples[idx_chains,] * out_bart$d_samples[idx_chains,]

meps_post <- meps %>%
  mutate(zeta_hat = colMeans(indv_zeta),
         delta_hat = colMeans(indv_delta)) %>%
  select(age, bmi, edu, income, povlev, region, sex, marital, race, seatbelt, delta_hat, zeta_hat)

rpart_fit <- rpart(delta_hat ~ . - zeta_hat, 
                   data = meps_post, 
                   control = rpart.control(cp = 0.03))

pdf(file = "figures/03_make_fig2.pdf", height = 5, width = 8)
rpart.plot(
  x = rpart_fit, 
  box.palette = "GnBu", # color scheme
  branch.lty = 3, # dotted branch lines
  shadow.col = "gray", # shadows under the node boxes
  nn = TRUE
)
dev.off()





