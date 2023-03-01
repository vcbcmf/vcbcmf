do_simulation_rlearn <- function(data,
                                 rlearner_fits,
                                 i_train,
                                 i_test,
                                 formula_m,
                                 formula_y,
                                 formula_ps,
                                 formula_y_lm,
                                 formula_m_lm,
                                 mediator_name,
                                 outcome_name,
                                 treat_name,
                                 n_iter,
                                 burnin,
                                 thin,
                                 seeds) {
  
  
  delta <- delta_hat <- delta_lower <- delta_true <- delta_upper <- leaf <-
    zeta <- zeta_hat <- zeta_lower <- zeta_true <- zeta_upper <- iteration <-
    fit_m <- fit_tau <- fit_delta <- w <- sigma_y <- sigma_m <- 
    NULL
  

  # Get training and testing set
  data_train <- data[i_train,]
  data_test <- data[i_test,]
  

  # # Store results
  # colnames_indv      <- c('seed', 'subj_id', 'zeta_mean', 'zeta_lower',
  #                         'zeta_upper', 'zeta_len', 'zeta_true', 'zeta_catch',
  #                         'delta_mean', 'delta_lower', 'delta_upper',
  #                         'delta_len', 'delta_true', 'delta_catch')
  # colnames_subgroup <- c('seed', 'group', 'zeta_mean', 'zeta_lower',
  #                        'zeta_upper', 'zeta_len', 'zeta_true', 'zeta_catch',
  #                        'delta_mean', 'delta_lower', 'delta_upper',
  #                        'delta_len', 'delta_true', 'delta_catch')
  # colnames_avg       <- c('seed', 'zeta_mean', 'zeta_lower', 'zeta_upper',
  #                         'zeta_len', 'zeta_true', 'zeta_catch',
  #                         'delta_mean', 'delta_lower', 'delta_upper',
  #                         'delta_len', 'delta_true', 'delta_catch')
  # colnames_tree_subgroups <- c('seed', 'group', 'zeta_mean', 'zeta_lower',
  #                              'zeta_upper', 'zeta_len', 'zeta_true',
  #                              'zeta_catch', 'delta_mean', 'delta_lower',
  #                              'delta_upper', 'delta_len', 'delta_true',
  #                              'delta_catch')
  
  n_reps <- length(seeds)
  
  for(i in 1:n_reps) {
    
    file_name_indv_rl <- paste0('Simulation RLEARN/indv_seed', seeds[i], '.rds')
    file_name_subgroup_rl <- paste0('Simulation RLEARN/subgroup_seed', seeds[i], '.rds')
    file_name_avg_rl <- paste0('Simulation RLEARN/avg_seed', seeds[i], '.rds')
    file_name_tree_subgroup_rl <- paste0('Simulation RLEARN/tree_subgroup_seed', seeds[i], '.rds')
    
    if (file.exists(file_name_indv_rl) & file.exists(file_name_subgroup_rl) &
        file.exists(file_name_avg_rl) & file.exists(file_name_tree_subgroup_rl))
    {
      
      indv_df_rl_i <- readRDS(file_name_indv_rl)
      subgroup_df_rl_i <- readRDS(file_name_subgroup_rl)
      avg_df_rl_i <- readRDS(file_name_avg_rl)
      tree_subgroups_list_rl_i <- readRDS(file_name_tree_subgroup_rl)
      
    } else {
      
      set.seed(seeds[i])
      
      # Get true effects
      indv_zeta_true <- rlearner_fits$zeta
      indv_delta_true <- rlearner_fits$delta
      avg_zeta_true <- mean(indv_zeta_true)
      avg_delta_true <- mean(indv_delta_true)
      
      subgroup_labels <- c('white, age < 34', 'non-white, age < 34', 'age ≥ 67',
                           'white, 34 ≤ age < 67', 'non-white, 34 ≤ age < 67')
      i1 <- which(data$age < 34 & data$race == 'White')
      i2 <- which(data$age < 34 & data$race != 'White')
      i3 <- which(data$age >= 67)
      i4 <- which(data$age >= 34 & data$age < 67 & data$race == 'White')
      i5 <- which(data$age >= 34 & data$age < 67 & data$race != 'White')
      
      group1_zeta_true <- mean(indv_zeta_true[i1])
      group2_zeta_true <- mean(indv_zeta_true[i2])
      group3_zeta_true <- mean(indv_zeta_true[i3])
      group4_zeta_true <- mean(indv_zeta_true[i4])
      group5_zeta_true <- mean(indv_zeta_true[i5])
      subgroup_zeta_true <- data.frame(
        group = subgroup_labels,
        direct = c(group1_zeta_true, group2_zeta_true,
                   group3_zeta_true, group4_zeta_true,
                   group5_zeta_true))

      group1_delta_true <- mean(indv_delta_true[i1])
      group2_delta_true <- mean(indv_delta_true[i2])
      group3_delta_true <- mean(indv_delta_true[i3])
      group4_delta_true <- mean(indv_delta_true[i4])
      group5_delta_true <- mean(indv_delta_true[i5])
      subgroup_delta_true <- data.frame(
        group = subgroup_labels,
        indirect = c(group1_delta_true, group2_delta_true,
                     group3_delta_true, group4_delta_true,
                     group5_delta_true))
      
       
      # Make training simulation set for rlearn truth
      new_data <- rlearner_fits %>% with(
        simulate_rlearn_med(fit_m = fit_m, fit_tau = fit_tau, 
                            fit_delta = fit_delta, w = w, sigma_y = sigma_y, 
                            sigma_m = sigma_m)
      )

       
      data_train_sim <- data_train
      data_train_sim[[mediator_name]] <- new_data$M[i_train]
      data_train_sim[[outcome_name]]  <- new_data$Y[i_train]

      
      ## BART Truth / BART Fit ----
      # Clever covariates for training/testing set
      clever_cov_rl <- get_clever_cov(data_train_sim, data, formula_m,
                                      mediator_name, outcome_name, treat_name)
      m0_hat_rl <- clever_cov_rl$m0_hat
      m1_hat_rl <- clever_cov_rl$m1_hat
      m0_hat_train_rl <- m0_hat_rl[i_train]
      m1_hat_train_rl <- m1_hat_rl[i_train]
      rm(clever_cov_rl)
      gc()
      
      # Propensity score for training/testing set
      pi_hat_rl <- get_ps(data, data, formula_ps)
      pi_hat_train_rl <- pi_hat_rl[i_train]
      
      # Get output for simulated dataset
      out_rl <- bart_mediate(data_train_sim, data,
                             formula_m, formula_y,
                             pi_hat_train_rl, pi_hat_rl,
                             m0_hat_train_rl, m0_hat_rl,
                             m1_hat_train_rl, m1_hat_rl,
                             'phealth', 'logY', 'smoke',
                             n_iter, burnin, thin)

                                              
      # Get simulated direct/indirect distributions
      zeta_rl <- out_rl$zeta_samples
      zeta_train_rl <- zeta_rl[,i_train]
      delta_rl <- out_rl$d_samples * out_rl$tau_samples
      delta_train_rl <- delta_rl[,i_train]
    
      rm(out_rl)
      gc()

      avg_zeta_rl    <- rowMeans(zeta_train_rl)
      avg_delta_rl  <- rowMeans(delta_train_rl)
      indv_zeta_rl   <- zeta_rl[,i_test]
      indv_delta_rl <- delta_rl[,i_test]
      
      i1_train <- which(data_train$age < 34 & data_train$race == 'White')
      i2_train <- which(data_train$age < 34 & data_train$race != 'White')
      i3_train <- which(data_train$age >= 67)
      i4_train <- which(data_train$age >= 34 & data_train$age < 67 & data_train$race == 'White')
      i5_train <- which(data_train$age >= 34 & data_train$age < 67 & data_train$race != 'White')
      
      group1_zeta_rl <- rowMeans(zeta_train_rl[,i1_train])
      group2_zeta_rl <- rowMeans(zeta_train_rl[,i2_train])
      group3_zeta_rl <- rowMeans(zeta_train_rl[,i3_train])
      group4_zeta_rl <- rowMeans(zeta_train_rl[,i4_train])
      group5_zeta_rl <- rowMeans(zeta_train_rl[,i5_train])
      subgroup_zeta_rl <- cbind(group1_zeta_rl, group2_zeta_rl, group3_zeta_rl,
                                group4_zeta_rl, group5_zeta_rl)
      
      group1_delta_rl <- rowMeans(delta_train_rl[,i1_train])
      group2_delta_rl <- rowMeans(delta_train_rl[,i2_train])
      group3_delta_rl <- rowMeans(delta_train_rl[,i3_train])
      group4_delta_rl <- rowMeans(delta_train_rl[,i4_train])
      group5_delta_rl <- rowMeans(delta_train_rl[,i5_train])
      subgroup_delta_rl <- cbind(group1_delta_rl, group2_delta_rl, group3_delta_rl,
                                 group4_delta_rl, group5_delta_rl)
      
      # Get direct/indirect credible intervals
      avg_zeta_interval_rl    <- quantile(avg_zeta_rl, probs = c(0.025, 0.975))
      avg_delta_interval_rl   <- quantile(avg_delta_rl, probs = c(0.025, 0.975))
      
      indv_zeta_interval_rl   <- colQuantiles(indv_zeta_rl, probs = c(0.025, 0.975))
      indv_delta_interval_rl  <- colQuantiles(indv_delta_rl, probs = c(0.025, 0.975))
      
      subgroup_zeta_interval_rl <- colQuantiles(subgroup_zeta_rl, probs = c(0.025, 0.975))
      subgroup_delta_interval_rl <- colQuantiles(subgroup_delta_rl, probs = c(0.025, 0.975))
      
      # Mean point estimates
      avg_zeta_mean_rl    <- mean(avg_zeta_rl)
      avg_delta_mean_rl  <- mean(avg_delta_rl)
      
      indv_zeta_mean_rl   <- colMeans(indv_zeta_rl)
      indv_delta_mean_rl <- colMeans(indv_delta_rl)
      
      subgroup_zeta_mean_rl <- colMeans(subgroup_zeta_rl)
      subgroup_delta_mean_rl <- colMeans(subgroup_delta_rl)
      
      # Catch T/F
      avg_zeta_catch_rl <- (avg_zeta_true >= avg_zeta_interval_rl[1]) &
        (avg_zeta_true <= avg_zeta_interval_rl[2])
      avg_delta_catch_rl <- (avg_delta_true >= avg_delta_interval_rl[1]) &
        (avg_delta_true <= avg_delta_interval_rl[2])
      
      indv_zeta_catch_rl <- unlist(lapply(
        X = 1:length(indv_zeta_true[i_test]),
        FUN = function(i) (indv_zeta_true[i_test][i] >= indv_zeta_interval_rl[i,1]) &
          (indv_zeta_true[i_test][i] <= indv_zeta_interval_rl[i,2])
      ))
      indv_delta_catch_rl <- unlist(lapply(
        X = 1:length(indv_delta_true[i_test]),
        FUN = function(i) (indv_delta_true[i_test][i] >= indv_delta_interval_rl[i,1]) &
          (indv_delta_true[i_test][i] <= indv_delta_interval_rl[i,2])
      ))
      
      subgroup_zeta_catch_rl <- unlist(lapply(
        X = 1:nrow(subgroup_zeta_true),
        FUN = function(i) (subgroup_zeta_true[i,2] >= subgroup_zeta_interval_rl[i,1]) &
          (subgroup_zeta_true[i,2] <= subgroup_zeta_interval_rl[i,2])
      ))
      subgroup_delta_catch_rl <- unlist(lapply(
        X = 1:nrow(subgroup_delta_true),
        FUN = function(i) (subgroup_delta_true[i,2] >= subgroup_delta_interval_rl[i,1]) &
          (subgroup_delta_true[i,2] <= subgroup_delta_interval_rl[i,2])
      ))
      
      # Saving
      # Individual
      indv_df_rl_i <- data.frame(
        seed = rep(seeds[i], nrow(data_test)),
        subj_id = 1:nrow(data_test),
        zeta_mean = indv_zeta_mean_rl,
        zeta_lower = indv_zeta_interval_rl[, 1],
        zeta_upper = indv_zeta_interval_rl[, 2],
        zeta_len = indv_zeta_interval_rl[, 2] - indv_zeta_interval_rl[, 1],
        zeta_true = indv_zeta_true[i_test],
        zeta_catch = as.numeric(indv_zeta_catch_rl),
        delta_mean = indv_delta_mean_rl,
        delta_lower = indv_delta_interval_rl[, 1],
        delta_upper = indv_delta_interval_rl[, 2],
        delta_len = indv_delta_interval_rl[, 2] - indv_delta_interval_rl[, 1],
        delta_true = indv_delta_true[i_test],
        delta_catch = as.numeric(indv_delta_catch_rl)
      )
      saveRDS(indv_df_rl_i, file_name_indv_rl)
      
      # Subgroup
      subgroup_df_rl_i <- data.frame(
        seed = rep(seeds[i], length(subgroup_labels)),
        group = subgroup_labels,
        zeta_mean = subgroup_zeta_mean_rl,
        zeta_lower = subgroup_zeta_interval_rl[, 1],
        zeta_upper = subgroup_zeta_interval_rl[, 2],
        zeta_len = subgroup_zeta_interval_rl[, 2] - subgroup_zeta_interval_rl[, 1],
        zeta_true = subgroup_zeta_true[, 2],
        zeta_catch = as.numeric(subgroup_zeta_catch_rl),
        delta_mean = subgroup_delta_mean_rl,
        delta_lower = subgroup_delta_interval_rl[, 1],
        delta_upper = subgroup_delta_interval_rl[, 2],
        delta_len = subgroup_delta_interval_rl[, 2] - subgroup_delta_interval_rl[, 1],
        delta_true = subgroup_delta_true[, 2],
        delta_catch = as.numeric(subgroup_delta_catch_rl)
      )
      saveRDS(subgroup_df_rl_i, file_name_subgroup_rl)
      
      # Average
      avg_df_rl_i <- data.frame(
        seed = seeds[i],
        zeta_mean = avg_zeta_mean_rl,
        zeta_lower = avg_zeta_interval_rl[1],
        zeta_upper = avg_zeta_interval_rl[2],
        zeta_len = avg_zeta_interval_rl[2] - avg_zeta_interval_rl[1],
        zeta_true = avg_zeta_true,
        zeta_catch = as.numeric(avg_zeta_catch_rl),
        delta_mean = avg_delta_mean_rl,
        delta_lower = avg_delta_interval_rl[1],
        delta_upper = avg_delta_interval_rl[2],
        delta_len = avg_delta_interval_rl[2] - avg_delta_interval_rl[1],
        delta_true = avg_delta_true,
        delta_catch = as.numeric(avg_delta_catch_rl)
      )
      saveRDS(avg_df_rl_i, file_name_avg_rl)
      
      # Tree Projection Subgroups
      data_post <- model.frame(formula_m, data = data) %>%
        select(-all_of(mediator_name)) %>%
        mutate(zeta_hat = colMeans(zeta_rl), delta_hat = colMeans(delta_rl))
      
      tree <- rpart(delta_hat ~ . - zeta_hat, data = data_post)
      
      # True subroup values
      true_leaf_df <- data.frame(zeta_true = indv_zeta_true,
                                 delta_true = indv_delta_true,
                                 leaf = tree$where)
      
      true_leaf_means <- true_leaf_df %>% group_by(leaf) %>%
        summarise(zeta_true = mean(zeta_true),
                  delta_true = mean(delta_true))
      
      # Simulated subgroup values
      sim_leaf_df <- data.frame(zeta = c(t(zeta_rl)),
                                delta = c(t(delta_rl)),
                                iteration = rep(1:n_iter, each = nrow(data)),
                                leaf = rep(tree$where, times = n_iter))
      
      sim_leaf_means_df <- sim_leaf_df %>% group_by(iteration, leaf) %>%
        summarize(zeta_hat = mean(zeta),
                  delta_hat = mean(delta))
      
      sim_leaf_quantiles <- sim_leaf_means_df %>% group_by(leaf) %>%
        summarize(zeta_mean = mean(zeta_hat),
                  zeta_lower = quantile(zeta_hat, 0.025),
                  zeta_upper = quantile(zeta_hat, 0.975),
                  zeta_len = zeta_upper - zeta_lower,
                  delta_mean = mean(delta_hat),
                  delta_lower = quantile(delta_hat, 0.025),
                  delta_upper = quantile(delta_hat, 0.975),
                  delta_len = delta_upper - delta_lower)
      
      tree_zeta_catch <- (true_leaf_means$zeta_true >= sim_leaf_quantiles$zeta_lower) &
        (true_leaf_means$zeta_true <= sim_leaf_quantiles$zeta_upper)
      tree_delta_catch <- (true_leaf_means$delta_true >= sim_leaf_quantiles$delta_lower) &
        (true_leaf_means$delta_true <= sim_leaf_quantiles$delta_upper)
      
      tree_subgroups_list_rl_i <- data.frame(seed = seeds[i],
                                             group = sim_leaf_quantiles$leaf,
                                             zeta_mean = sim_leaf_quantiles$zeta_mean,
                                             zeta_lower = sim_leaf_quantiles$zeta_lower,
                                             zeta_upper = sim_leaf_quantiles$zeta_upper,
                                             zeta_len = sim_leaf_quantiles$zeta_len,
                                             zeta_true = true_leaf_means$zeta_true,
                                             zeta_catch = as.numeric(tree_zeta_catch),
                                             delta_mean = sim_leaf_quantiles$delta_mean,
                                             delta_lower = sim_leaf_quantiles$delta_lower,
                                             delta_upper = sim_leaf_quantiles$delta_upper,
                                             delta_len = sim_leaf_quantiles$delta_len,
                                             delta_true = true_leaf_means$delta_true,
                                             delta_catch = as.numeric(tree_delta_catch))
      saveRDS(tree_subgroups_list_rl_i, file_name_tree_subgroup_rl)
      rm(data_post, true_leaf_df, sim_leaf_df, sim_leaf_means_df, sim_leaf_means_df)
      gc()
      
      ## Clean ----
      
      rm(zeta_rl, zeta_train_rl, delta_rl, delta_train_rl,
         indv_zeta_rl, indv_delta_rl)
      
      gc()
      
    }
    
  }
    
}
