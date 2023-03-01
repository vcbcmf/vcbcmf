do_simulation_rlearn_lsem <- function(data, 
                                      rlearner_fits,
                                      i_train,
                                      i_test,
                                      formula_m_lsem,
                                      formula_y_lsem,
                                      formula_X_lsem,
                                      mediator_name,
                                      outcome_name,
                                      treat_name,
                                      seeds
                                      ) {
  
  fit_delta <- fit_m <- fit_tau <- sigma_m <- sigma_y <- w <- NULL
  
  # Get training and testing set
  data_train <- data[i_train,]
  data_test <- data[i_test,]
  
  n_reps <- length(seeds)
  
  for(i in 1:n_reps) {
    
    cat(paste("Starting Seed", seeds[i], "\n"))
    
    file_name_indv_rl <- paste0('Simulation RLEARN LSEM/indv_seed', seeds[i], '.rds')
    file_name_subgroup_rl <- paste0('Simulation RLEARN LSEM/subgroup_seed', seeds[i], '.rds')
    file_name_avg_rl <- paste0('Simulation RLEARN LSEM/avg_seed', seeds[i], '.rds')
    # file_name_tree_subgroup_rl <- paste0('Simulation RLEARN LSEM/tree_subgroup_seed', seeds[i], '.rds')
   
    if (file.exists(file_name_indv_rl) & file.exists(file_name_subgroup_rl) &
        file.exists(file_name_avg_rl)) {
      
      indv_df_ll_i <- readRDS(file_name_indv_rl)
      subgroup_df_ll_i <- readRDS(file_name_subgroup_rl)
      avg_df_ll_i <- readRDS(file_name_avg_rl)
      # tree_subgroups_list_ll_i <- readRDS(file_name_tree_subgroup_rl)
      
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
      
      ## LSEM Truth / LSEM Fit ----
      # Fit LSEM to LSEM simulated data and bootstrap for zeta/delta
      fit_m_sim_ll <- lm(formula_m_lsem, data = data_train_sim)
      fit_y_sim_ll <- lm(formula_y_lsem, data = data_train_sim)
      out_boot_ll  <- replicate(1000, lsem_boot(fit_m_sim_ll, fit_y_sim_ll,
                                                formula_X_lsem,
                                                data_train_sim, data,
                                                mediator_name, treat_name,
                                                outcome_name))
      zeta_boot_ll  <- matrix(unlist(out_boot_ll['zeta',]), ncol = nrow(data), byrow = TRUE)
      delta_boot_ll <- matrix(unlist(out_boot_ll['delta',]), ncol = nrow(data), byrow = TRUE)
      zeta_train_ll <- zeta_boot_ll[,i_train]
      delta_train_ll <- delta_boot_ll[,i_train]
      
      indv_zeta_ll <- zeta_boot_ll[,i_test]
      indv_delta_ll <- delta_boot_ll[,i_test]
      avg_zeta_ll   <- rowMeans(zeta_train_ll)
      avg_delta_ll  <- rowMeans(delta_train_ll)
      
      i1_train <- which(data_train$age < 34 & data_train$race == 'White')
      i2_train <- which(data_train$age < 34 & data_train$race != 'White')
      i3_train <- which(data_train$age >= 67)
      i4_train <- which(data_train$age >= 34 & data_train$age < 67 & data_train$race == 'White')
      i5_train <- which(data_train$age >= 34 & data_train$age < 67 & data_train$race != 'White')
      
      group1_zeta_ll <- rowMeans(zeta_train_ll[,i1_train])
      group2_zeta_ll <- rowMeans(zeta_train_ll[,i2_train])
      group3_zeta_ll <- rowMeans(zeta_train_ll[,i3_train])
      group4_zeta_ll <- rowMeans(zeta_train_ll[,i4_train])
      group5_zeta_ll <- rowMeans(zeta_train_ll[,i5_train])
      subgroup_zeta_ll <- cbind(group1_zeta_ll, group2_zeta_ll,
                                group3_zeta_ll, group4_zeta_ll,
                                group5_zeta_ll)
      
      group1_delta_ll <- rowMeans(delta_train_ll[,i1_train])
      group2_delta_ll <- rowMeans(delta_train_ll[,i2_train])
      group3_delta_ll <- rowMeans(delta_train_ll[,i3_train])
      group4_delta_ll <- rowMeans(delta_train_ll[,i4_train])
      group5_delta_ll <- rowMeans(delta_train_ll[,i5_train])
      subgroup_delta_ll <- cbind(group1_delta_ll, group2_delta_ll,
                                 group3_delta_ll, group4_delta_ll,
                                 group5_delta_ll)
      
      zeta_bar_ll <- colMeans(indv_zeta_ll)
      delta_bar_ll <- colMeans(indv_delta_ll)
      se_zeta_ll <- apply(indv_zeta_ll, 2, sd)
      se_delta_ll <- apply(indv_delta_ll, 2, sd)
      t_stat_train <- qt(0.975, nrow(data_train) - 1)
      t_stat_test <- qt(0.975, nrow(data_test) - 1)
      indv_zeta_interval_ll  <- t(apply(cbind(zeta_bar_ll + t_stat_test * se_zeta_ll, zeta_bar_ll - t_stat_test * se_zeta_ll), 1, sort))
      indv_delta_interval_ll <- t(apply(cbind(delta_bar_ll - t_stat_test * se_delta_ll, delta_bar_ll + t_stat_test * se_delta_ll), 1, sort))
      avg_zeta_interval_ll   <- mean(avg_zeta_ll) + c(-1, 1) * t_stat_train * sd(avg_zeta_ll)
      avg_delta_interval_ll  <- mean(avg_delta_ll) + c(-1, 1) * t_stat_train * sd(avg_delta_ll)
      
      se_zeta_subgroup_ll <- apply(subgroup_zeta_ll, 2, sd)
      se_delta_subgroup_ll <- apply(subgroup_delta_ll, 2, sd)
      subgroup_zeta_interval_ll <- t(apply(cbind(colMeans(subgroup_zeta_ll) + t_stat_train * se_zeta_subgroup_ll, colMeans(subgroup_zeta_ll) - t_stat_train * se_zeta_subgroup_ll), 1, sort))
      subgroup_delta_interval_ll <- t(apply(cbind(colMeans(subgroup_delta_ll) + t_stat_train * se_delta_subgroup_ll, colMeans(subgroup_delta_ll) - t_stat_train * se_delta_subgroup_ll), 1, sort))
      
      # Mean estimates
      indv_zeta_mean_ll  <- colMeans(indv_zeta_ll)
      indv_delta_mean_ll <- colMeans(indv_delta_ll)
      avg_zeta_mean_ll   <- mean(avg_zeta_ll)
      avg_delta_mean_ll  <- mean(avg_delta_ll)
      
      subgroup_zeta_mean_ll <- colMeans(subgroup_zeta_ll)
      subgroup_delta_mean_ll <- colMeans(subgroup_delta_ll)
      
      # Catch?
      indv_zeta_catch_ll <- unlist(lapply(1:length(indv_zeta_true[i_test]),
                                          function(i) (indv_zeta_true[i_test][i] >= indv_zeta_interval_ll[i,1]) &
                                            (indv_zeta_true[i_test][i] <= indv_zeta_interval_ll[i,2])))
      indv_delta_catch_ll <- unlist(lapply(1:length(indv_delta_true[i_test]),
                                           function(i) (indv_delta_true[i_test][i] >= indv_delta_interval_ll[i,1]) &
                                             (indv_delta_true[i_test][i] <= indv_delta_interval_ll[i,2])))
      avg_zeta_catch_ll  <- (avg_zeta_true >= avg_zeta_interval_ll[1]) & (avg_zeta_true <= avg_zeta_interval_ll[2])
      avg_delta_catch_ll <- (avg_delta_true >= avg_delta_interval_ll[1]) & (avg_delta_true <= avg_delta_interval_ll[2])
      
      subgroup_zeta_catch_ll <- unlist(lapply(1:nrow(subgroup_zeta_true),
                                              function(i) (subgroup_zeta_true[i,2] >= subgroup_zeta_interval_ll[i,1]) &
                                                (subgroup_zeta_true[i,2] <= subgroup_zeta_interval_ll[i,2])))
      subgroup_delta_catch_ll <- unlist(lapply(1:nrow(subgroup_delta_true),
                                               function(i) (subgroup_delta_true[i,2] >= subgroup_delta_interval_ll[i,1]) &
                                                 (subgroup_delta_true[i,2] <= subgroup_delta_interval_ll[i,2])))
      
      # Save
      # Individual
      indv_df_ll_i <- data.frame(seed = rep(seeds[i], nrow(data_test)),
                                 subj_id = 1:nrow(data_test),
                                 zeta_mean = indv_zeta_mean_ll,
                                 zeta_lower = indv_zeta_interval_ll[,1],
                                 zeta_upper = indv_zeta_interval_ll[,2],
                                 zeta_len = indv_zeta_interval_ll[,2] - indv_zeta_interval_ll[,1],
                                 zeta_true = indv_zeta_true[i_test],
                                 zeta_catch = as.numeric(indv_zeta_catch_ll),
                                 delta_mean = indv_delta_mean_ll,
                                 delta_lower = indv_delta_interval_ll[,1],
                                 delta_upper = indv_delta_interval_ll[,2],
                                 delta_len = indv_delta_interval_ll[,2] - indv_delta_interval_ll[,1],
                                 delta_true = indv_delta_true[i_test],
                                 delta_catch = as.numeric(indv_delta_catch_ll))
      saveRDS(indv_df_ll_i, file_name_indv_rl)
      
      # Average
      avg_df_ll_i <- data.frame(seed = seeds[i],
                                zeta_mean = avg_zeta_mean_ll,
                                zeta_lower = avg_zeta_interval_ll[1],
                                zeta_upper = avg_zeta_interval_ll[2],
                                zeta_len = avg_zeta_interval_ll[2] - avg_zeta_interval_ll[1],
                                zeta_true = avg_zeta_true,
                                zeta_catch = as.numeric(avg_zeta_catch_ll),
                                delta_mean = avg_delta_mean_ll,
                                delta_lower = avg_delta_interval_ll[1],
                                delta_upper = avg_delta_interval_ll[2],
                                delta_len = avg_delta_interval_ll[2] - avg_delta_interval_ll[1],
                                delta_true = avg_delta_true,
                                delta_catch = as.numeric(avg_delta_catch_ll))
      saveRDS(avg_df_ll_i, file_name_avg_rl)
      
      # Subgroup
      subgroup_df_ll_i <- data.frame(seed = rep(seeds[i], length(subgroup_labels)),
                                     group = subgroup_labels,
                                     zeta_mean = subgroup_zeta_mean_ll,
                                     zeta_lower = subgroup_zeta_interval_ll[,1],
                                     zeta_upper = subgroup_zeta_interval_ll[,2],
                                     zeta_len = subgroup_zeta_interval_ll[,2] - subgroup_zeta_interval_ll[,1],
                                     zeta_true = subgroup_zeta_true[,2],
                                     zeta_catch = as.numeric(subgroup_zeta_catch_ll),
                                     delta_mean = subgroup_delta_mean_ll,
                                     delta_lower = subgroup_delta_interval_ll[,1],
                                     delta_upper = subgroup_delta_interval_ll[,2],
                                     delta_len = subgroup_delta_interval_ll[,2] - subgroup_delta_interval_ll[,1],
                                     delta_true = subgroup_delta_true[,2],
                                     delta_catch = as.numeric(subgroup_delta_catch_ll))
      saveRDS(subgroup_df_ll_i, file_name_subgroup_rl)
      
    }
  }
}

