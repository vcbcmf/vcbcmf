do_simulation_lsem <- function(data, i_train, i_test,
                               fit_m_lsem, fit_y_lsem, formula_X_lsem,
                               formula_m_lsem, formula_y_lsem,
                               mediator_name, outcome_name, treat_name,
                               zeta_hat_lsem, delta_hat_lsem,
                               sigma_m_hat_lsem, sigma_y_hat_lsem,
                               mu_y_hat_bart, zeta_hat_bart, d_hat_bart,
                               mu_m_hat_bart, tau_hat_bart,
                               sigma_y_hat_bart, sigma_m_hat_bart,
                               n_reps, seeds) {
  
  # Get training and testing set
  data_train <- data[i_train,]
  data_test <- data[i_test,]
  
  # Store results
  # colnames_indv      <- c('seed', 'subj_id', 'zeta_mean', 'zeta_lower',
  #                         'zeta_upper', 'zeta_len', 'zeta_true', 'zeta_catch',
  #                         'delta_mean', 'delta_lower', 'delta_upper',
  #                         'delta_len', 'delta_true', 'delta_catch')
  # colnames_avg       <- c('seed', 'zeta_mean', 'zeta_lower', 'zeta_upper',
  #                         'zeta_len', 'zeta_true', 'zeta_catch',
  #                         'delta_mean', 'delta_lower', 'delta_upper',
  #                         'delta_len', 'delta_true', 'delta_catch')
  # colnames_subgroup <- c('seed', 'group', 'zeta_mean', 'zeta_lower',
  #                        'zeta_upper', 'zeta_len', 'zeta_true', 'zeta_catch',
  #                        'delta_mean', 'delta_lower', 'delta_upper',
  #                        'delta_len', 'delta_true', 'delta_catch')
  
  # LSEM/LSEM
  # indv_df_ll      <- data.frame(matrix(NA, nrow = n_reps * nrow(data_test),
  #                                      ncol = length(colnames_indv),
  #                                      dimnames = list(c(), colnames_indv)))
  # avg_df_ll       <- data.frame(matrix(NA, nrow = n_reps,
  #                                      ncol = length(colnames_avg),
  #                                      dimnames = list(c(), colnames_avg)))
  # subgroup_df_ll <- data.frame(matrix(NA, nrow = n_reps * 5,
  #                                     ncol = length(colnames_subgroup),
  #                                     dimnames = list(c(), colnames_subgroup)))
  
  # BART/LSEM
  # indv_df_bl      <- data.frame(matrix(NA, nrow = n_reps * nrow(data_test),
  #                                      ncol = length(colnames_indv),
  #                                      dimnames = list(c(), colnames_indv)))
  # avg_df_bl       <- data.frame(matrix(NA, nrow = n_reps,
  #                                      ncol = length(colnames_avg),
  #                                      dimnames = list(c(), colnames_avg)))
  # subgroup_df_bl <- data.frame(matrix(NA, nrow = n_reps * 5,
  #                                     ncol = length(colnames_subgroup),
  #                                     dimnames = list(c(), colnames_subgroup)))
  
  for (i in 1:n_reps) {
    
    
    file_name_indv_ll          <- paste0('Simulation LSEMLSEM/indv_seed', seeds[i], '.rds')
    file_name_subgroup_ll      <- paste0('Simulation LSEMLSEM/subgroup_seed', seeds[i], '.rds')
    file_name_avg_ll           <- paste0('Simulation LSEMLSEM/avg_seed', seeds[i], '.rds')
    file_name_indv_bl          <- paste0('Simulation BARTLSEM/indv_seed', seeds[i], '.rds')
    file_name_subgroup_bl      <- paste0('Simulation BARTLSEM/subgroup_seed', seeds[i], '.rds')
    file_name_avg_bl           <- paste0('Simulation BARTLSEM/avg_seed', seeds[i], '.rds')
    
    skip_sim <- 
      file.exists(file_name_indv_ll) & 
      file.exists(file_name_subgroup_ll) & 
      file.exists(file_name_avg_ll) & 
      file.exists(file_name_indv_bl) & 
      file.exists(file_name_avg_bl) & 
      file.exists(file_name_subgroup_bl)
    
    if (skip_sim) {
      
      indv_df_ll_i          <- readRDS(file_name_indv_ll)
      subgroup_df_ll_i      <- readRDS(file_name_subgroup_ll)
      avg_df_ll_i           <- readRDS(file_name_avg_ll)
      
      indv_df_bl_i          <- readRDS(file_name_indv_bl)
      subgroup_df_bl_i      <- readRDS(file_name_subgroup_bl)
      avg_df_bl_i           <- readRDS(file_name_avg_bl)
      
      
    } else {
    
    
      set.seed(seeds[i])
      
      cat(paste("Starting seed", seeds[i], "\n"))
      
      # Get true effects for LSEM
      indv_zeta_true_lsem  <- zeta_hat_lsem
      indv_delta_true_lsem <- delta_hat_lsem
      avg_zeta_true_lsem   <- mean(indv_zeta_true_lsem)
      avg_delta_true_lsem  <- mean(indv_delta_true_lsem)
      
      subgroup_labels <- c('white, age < 34', 'non-white, age < 34', 'age ≥ 67',
                           'white, 34 ≤ age < 67', 'non-white, 34 ≤ age < 67')
      i1 <- which(data$age < 34 & data$race == 'White')
      i2 <- which(data$age < 34 & data$race != 'White')
      i3 <- which(data$age >= 67)
      i4 <- which(data$age >= 34 & data$age < 67 & data$race == 'White')
      i5 <- which(data$age >= 34 & data$age < 67 & data$race != 'White')
      
      group1_zeta_true_lsem <- mean(indv_zeta_true_lsem[i1])
      group2_zeta_true_lsem <- mean(indv_zeta_true_lsem[i2])
      group3_zeta_true_lsem <- mean(indv_zeta_true_lsem[i3])
      group4_zeta_true_lsem <- mean(indv_zeta_true_lsem[i4])
      group5_zeta_true_lsem <- mean(indv_zeta_true_lsem[i5])
      subgroup_zeta_true_lsem <- data.frame(group = subgroup_labels,
                                            direct = c(group1_zeta_true_lsem, group2_zeta_true_lsem,
                                                       group3_zeta_true_lsem, group4_zeta_true_lsem,
                                                       group5_zeta_true_lsem))
      
      group1_delta_true_lsem <- mean(indv_delta_true_lsem[i1])
      group2_delta_true_lsem <- mean(indv_delta_true_lsem[i2])
      group3_delta_true_lsem <- mean(indv_delta_true_lsem[i3])
      group4_delta_true_lsem <- mean(indv_delta_true_lsem[i4])
      group5_delta_true_lsem <- mean(indv_delta_true_lsem[i5])
      subgroup_delta_true_lsem <- data.frame(group = subgroup_labels,
                                             indirect = c(group1_delta_true_lsem, group2_delta_true_lsem,
                                                          group3_delta_true_lsem, group4_delta_true_lsem,
                                                          group5_delta_true_lsem))
      
      # Get true effects for BART
      indv_zeta_true_bart   <- zeta_hat_bart
      indv_delta_true_bart  <- d_hat_bart * tau_hat_bart
      avg_zeta_true_bart    <- mean(indv_zeta_true_bart)
      avg_delta_true_bart   <- mean(indv_delta_true_bart)
      
      group1_zeta_true_bart <- mean(indv_zeta_true_bart[i1])
      group2_zeta_true_bart <- mean(indv_zeta_true_bart[i2])
      group3_zeta_true_bart <- mean(indv_zeta_true_bart[i3])
      group4_zeta_true_bart <- mean(indv_zeta_true_bart[i4])
      group5_zeta_true_bart <- mean(indv_zeta_true_bart[i5])
      subgroup_zeta_true_bart <- data.frame(group = subgroup_labels,
                                            direct = c(group1_zeta_true_bart, group2_zeta_true_bart,
                                                       group3_zeta_true_bart, group4_zeta_true_bart,
                                                       group5_zeta_true_bart))
      
      group1_delta_true_bart <- mean(indv_delta_true_bart[i1])
      group2_delta_true_bart <- mean(indv_delta_true_bart[i2])
      group3_delta_true_bart <- mean(indv_delta_true_bart[i3])
      group4_delta_true_bart <- mean(indv_delta_true_bart[i4])
      group5_delta_true_bart <- mean(indv_delta_true_bart[i5])
      subgroup_delta_true_bart <- data.frame(group = subgroup_labels,
                                             indirect = c(group1_delta_true_bart, group2_delta_true_bart,
                                                          group3_delta_true_bart, group4_delta_true_bart,
                                                          group5_delta_true_bart))
      
      # Make training simulation set for LSEM truth
      epsilon_m_lsem <- rnorm(nrow(data_train))
      epsilon_y_lsem <- rnorm(nrow(data_train))
      
      m_sim_lsem <- predict(fit_m_lsem, newdata = data_train) + sigma_m_hat_lsem * epsilon_m_lsem
      Y_sim_lsem <- predict(fit_y_lsem, newdata = data_train %>% mutate(!!mediator_name := m_sim_lsem)) + sigma_y_hat_lsem * epsilon_y_lsem
      
      data_train_sim_lsem <- data_train
      data_train_sim_lsem[[mediator_name]] <- m_sim_lsem
      data_train_sim_lsem[[outcome_name]] <- Y_sim_lsem
      
      # Make training simulation set for BART truth
      mu_y_hat_train_bart <- mu_y_hat_bart[i_train]
      zeta_hat_train_bart <- zeta_hat_bart[i_train]
      d_hat_train_bart    <- d_hat_bart[i_train]
      mu_m_hat_train_bart <- mu_m_hat_bart[i_train]
      tau_hat_train_bart  <- tau_hat_bart[i_train]
      
      epsilon_y_bart <- rnorm(nrow(data_train))
      epsilon_m_bart <- rnorm(nrow(data_train))
      A_train <- data_train[[treat_name]]
      
      m_sim_bart <- mu_m_hat_train_bart + A_train * tau_hat_train_bart + sigma_m_hat_bart * epsilon_m_bart
      Y_sim_bart <- mu_y_hat_train_bart + A_train * zeta_hat_train_bart + m_sim_bart * d_hat_train_bart + sigma_y_hat_bart * epsilon_y_bart
      
      data_train_sim_bart <- data_train
      data_train_sim_bart[[mediator_name]] <- m_sim_bart
      data_train_sim_bart[[outcome_name]]  <- Y_sim_bart
      
      ## LSEM Truth / LSEM Fit ----
      # Fit LSEM to LSEM simulated data and bootstrap for zeta/delta
      fit_m_sim_ll <- lm(formula_m_lsem, data = data_train_sim_lsem)
      fit_y_sim_ll <- lm(formula_y_lsem, data = data_train_sim_lsem)
      out_boot_ll  <- replicate(1000, lsem_boot(fit_m_sim_ll, fit_y_sim_ll,
                                                formula_X_lsem,
                                                data_train_sim_lsem, data,
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
      
      # CI
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
      indv_zeta_catch_ll <- unlist(lapply(1:length(indv_zeta_true_lsem[i_test]),
                                          function(i) (indv_zeta_true_lsem[i_test][i] >= indv_zeta_interval_ll[i,1]) &
                                            (indv_zeta_true_lsem[i_test][i] <= indv_zeta_interval_ll[i,2])))
      indv_delta_catch_ll <- unlist(lapply(1:length(indv_delta_true_lsem[i_test]),
                                           function(i) (indv_delta_true_lsem[i_test][i] >= indv_delta_interval_ll[i,1]) &
                                             (indv_delta_true_lsem[i_test][i] <= indv_delta_interval_ll[i,2])))
      avg_zeta_catch_ll  <- (avg_zeta_true_lsem >= avg_zeta_interval_ll[1]) & (avg_zeta_true_lsem <= avg_zeta_interval_ll[2])
      avg_delta_catch_ll <- (avg_delta_true_lsem >= avg_delta_interval_ll[1]) & (avg_delta_true_lsem <= avg_delta_interval_ll[2])
      
      subgroup_zeta_catch_ll <- unlist(lapply(1:nrow(subgroup_zeta_true_lsem),
                                              function(i) (subgroup_zeta_true_lsem[i,2] >= subgroup_zeta_interval_ll[i,1]) &
                                                (subgroup_zeta_true_lsem[i,2] <= subgroup_zeta_interval_ll[i,2])))
      subgroup_delta_catch_ll <- unlist(lapply(1:nrow(subgroup_delta_true_lsem),
                                               function(i) (subgroup_delta_true_lsem[i,2] >= subgroup_delta_interval_ll[i,1]) &
                                                 (subgroup_delta_true_lsem[i,2] <= subgroup_delta_interval_ll[i,2])))
      
      # Save
      # Individual
      indv_df_ll_i <- data.frame(seed = rep(seeds[i], nrow(data_test)),
                                 subj_id = 1:nrow(data_test),
                                 zeta_mean = indv_zeta_mean_ll,
                                 zeta_lower = indv_zeta_interval_ll[,1],
                                 zeta_upper = indv_zeta_interval_ll[,2],
                                 zeta_len = indv_zeta_interval_ll[,2] - indv_zeta_interval_ll[,1],
                                 zeta_true = indv_zeta_true_lsem[i_test],
                                 zeta_catch = as.numeric(indv_zeta_catch_ll),
                                 delta_mean = indv_delta_mean_ll,
                                 delta_lower = indv_delta_interval_ll[,1],
                                 delta_upper = indv_delta_interval_ll[,2],
                                 delta_len = indv_delta_interval_ll[,2] - indv_delta_interval_ll[,1],
                                 delta_true = indv_delta_true_lsem[i_test],
                                 delta_catch = as.numeric(indv_delta_catch_ll))
      saveRDS(indv_df_ll_i, file = file_name_indv_ll)
      
      # Average
      avg_df_ll_i <- data.frame(seed = seeds[i],
                                zeta_mean = avg_zeta_mean_ll,
                                zeta_lower = avg_zeta_interval_ll[1],
                                zeta_upper = avg_zeta_interval_ll[2],
                                zeta_len = avg_zeta_interval_ll[2] - avg_zeta_interval_ll[1],
                                zeta_true = avg_zeta_true_lsem,
                                zeta_catch = as.numeric(avg_zeta_catch_ll),
                                delta_mean = avg_delta_mean_ll,
                                delta_lower = avg_delta_interval_ll[1],
                                delta_upper = avg_delta_interval_ll[2],
                                delta_len = avg_delta_interval_ll[2] - avg_delta_interval_ll[1],
                                delta_true = avg_delta_true_lsem,
                                delta_catch = as.numeric(avg_delta_catch_ll))
      saveRDS(avg_df_ll_i, file = file_name_avg_ll)
      
      # Subgroup
      subgroup_df_ll_i <- data.frame(seed = rep(seeds[i], length(subgroup_labels)),
                                     group = subgroup_labels,
                                     zeta_mean = subgroup_zeta_mean_ll,
                                     zeta_lower = subgroup_zeta_interval_ll[,1],
                                     zeta_upper = subgroup_zeta_interval_ll[,2],
                                     zeta_len = subgroup_zeta_interval_ll[,2] - subgroup_zeta_interval_ll[,1],
                                     zeta_true = subgroup_zeta_true_lsem[,2],
                                     zeta_catch = as.numeric(subgroup_zeta_catch_ll),
                                     delta_mean = subgroup_delta_mean_ll,
                                     delta_lower = subgroup_delta_interval_ll[,1],
                                     delta_upper = subgroup_delta_interval_ll[,2],
                                     delta_len = subgroup_delta_interval_ll[,2] - subgroup_delta_interval_ll[,1],
                                     delta_true = subgroup_delta_true_lsem[,2],
                                     delta_catch = as.numeric(subgroup_delta_catch_ll))
      saveRDS(object = subgroup_df_ll_i, file_name_subgroup_ll)
      
      # Save replications into a data frame
      # indv_df_ll[1:nrow(data_test) + nrow(data_test) * (i-1),] <- indv_df_ll_i
      # avg_df_ll[i,] <- avg_df_ll_i
      # subgroup_df_ll[1:5 + 5 * (i-1),] <- subgroup_df_ll_i
      
      
      ## BART Truth / LSEM Fit ----
      # Fit LSEM to LSEM simulated data and bootstrap for zeta/delta
      fit_m_sim_bl <- lm(formula_m_lsem, data = data_train_sim_bart)
      fit_y_sim_bl <- lm(formula_y_lsem, data = data_train_sim_bart)
      out_boot_bl  <- replicate(1000, lsem_boot(fit_m_sim_bl, fit_y_sim_bl,
                                                formula_X_lsem,
                                                data_train_sim_bart, data,
                                                mediator_name, treat_name,
                                                outcome_name))
      zeta_boot_bl  <- matrix(unlist(out_boot_bl['zeta',]), ncol = nrow(data), byrow = TRUE)
      delta_boot_bl <- matrix(unlist(out_boot_bl['delta',]), ncol = nrow(data), byrow = TRUE)
      zeta_train_bl <- zeta_boot_bl[,i_train]
      delta_train_bl <- delta_boot_bl[,i_train]
      
      indv_zeta_bl <- zeta_boot_bl[,i_test]
      indv_delta_bl <- delta_boot_bl[,i_test]
      avg_zeta_bl   <- rowMeans(zeta_train_bl)
      avg_delta_bl  <- rowMeans(delta_train_bl)
      
      group1_zeta_bl <- rowMeans(zeta_train_bl[,i1_train])
      group2_zeta_bl <- rowMeans(zeta_train_bl[,i2_train])
      group3_zeta_bl <- rowMeans(zeta_train_bl[,i3_train])
      group4_zeta_bl <- rowMeans(zeta_train_bl[,i4_train])
      group5_zeta_bl <- rowMeans(zeta_train_bl[,i5_train])
      subgroup_zeta_bl <- cbind(group1_zeta_bl, group2_zeta_bl,
                                group3_zeta_bl, group4_zeta_bl,
                                group5_zeta_bl)
      
      group1_delta_bl <- rowMeans(delta_train_bl[,i1_train])
      group2_delta_bl <- rowMeans(delta_train_bl[,i2_train])
      group3_delta_bl <- rowMeans(delta_train_bl[,i3_train])
      group4_delta_bl <- rowMeans(delta_train_bl[,i4_train])
      group5_delta_bl <- rowMeans(delta_train_bl[,i5_train])
      subgroup_delta_bl <- cbind(group1_delta_bl, group2_delta_bl,
                                 group3_delta_bl, group4_delta_bl,
                                 group5_delta_bl)
      
      # CI
      zeta_bar_bl <- colMeans(indv_zeta_bl)
      delta_bar_bl <- colMeans(indv_delta_bl)
      se_zeta_bl <- apply(indv_zeta_bl, 2, sd)
      se_delta_bl <- apply(indv_delta_bl, 2, sd)
      indv_zeta_interval_bl  <- t(apply(cbind(zeta_bar_bl + t_stat_test * se_zeta_bl, zeta_bar_bl - t_stat_test * se_zeta_bl), 1, sort))
      indv_delta_interval_bl <- t(apply(cbind(delta_bar_bl - t_stat_test * se_delta_bl, delta_bar_bl + t_stat_test * se_delta_bl), 1, sort))
      avg_zeta_interval_bl   <- mean(avg_zeta_bl) + c(-1, 1) * t_stat_train * sd(avg_zeta_bl)
      avg_delta_interval_bl  <- mean(avg_delta_bl) + c(-1, 1) * t_stat_train * sd(avg_delta_bl)
      
      se_zeta_subgroup_bl <- apply(subgroup_zeta_bl, 2, sd)
      se_delta_subgroup_bl <- apply(subgroup_delta_bl, 2, sd)
      subgroup_zeta_interval_bl <- t(apply(cbind(colMeans(subgroup_zeta_bl) + t_stat_train * se_zeta_subgroup_bl, colMeans(subgroup_zeta_bl) - t_stat_train * se_zeta_subgroup_bl), 1, sort))
      subgroup_delta_interval_bl <- t(apply(cbind(colMeans(subgroup_delta_bl) + t_stat_train * se_delta_subgroup_bl, colMeans(subgroup_delta_bl) - t_stat_train * se_delta_subgroup_bl), 1, sort))
      
      # Mean estimates
      indv_zeta_mean_bl  <- colMeans(indv_zeta_bl)
      indv_delta_mean_bl <- colMeans(indv_delta_bl)
      avg_zeta_mean_bl   <- mean(avg_zeta_bl)
      avg_delta_mean_bl  <- mean(avg_delta_bl)
      
      subgroup_zeta_mean_bl <- colMeans(subgroup_zeta_bl)
      subgroup_delta_mean_bl <- colMeans(subgroup_delta_bl)
      
      # Catch?
      indv_zeta_catch_bl <- unlist(lapply(1:length(indv_zeta_true_bart[i_test]),
                                          function(i) (indv_zeta_true_bart[i_test][i] >= indv_zeta_interval_bl[i,1]) &
                                            (indv_zeta_true_bart[i_test][i] <= indv_zeta_interval_bl[i,2])))
      indv_delta_catch_bl <- unlist(lapply(1:length(indv_delta_true_bart[i_test]),
                                           function(i) (indv_delta_true_bart[i_test][i] >= indv_delta_interval_bl[i,1]) &
                                             (indv_delta_true_bart[i_test][i] <= indv_delta_interval_bl[i,2])))
      avg_zeta_catch_bl  <- (avg_zeta_true_bart >= avg_zeta_interval_bl[1]) & (avg_zeta_true_bart <= avg_zeta_interval_bl[2])
      avg_delta_catch_bl <- (avg_delta_true_bart >= avg_delta_interval_bl[1]) & (avg_delta_true_bart <= avg_delta_interval_bl[2])
      
      subgroup_zeta_catch_bl <- unlist(lapply(1:nrow(subgroup_zeta_true_bart),
                                              function(i) (subgroup_zeta_true_bart[i,2] >= subgroup_zeta_interval_bl[i,1]) &
                                                (subgroup_zeta_true_bart[i,2] <= subgroup_zeta_interval_bl[i,2])))
      subgroup_delta_catch_bl <- unlist(lapply(1:nrow(subgroup_delta_true_bart),
                                               function(i) (subgroup_delta_true_bart[i,2] >= subgroup_delta_interval_bl[i,1]) &
                                                 (subgroup_delta_true_bart[i,2] <= subgroup_delta_interval_bl[i,2])))
      
      # Save
      # Individual
      indv_df_bl_i <- data.frame(seed = rep(seeds[i], nrow(data_test)),
                                 subj_id = 1:nrow(data_test),
                                 zeta_mean = indv_zeta_mean_bl,
                                 zeta_lower = indv_zeta_interval_bl[,1],
                                 zeta_upper = indv_zeta_interval_bl[,2],
                                 zeta_len = indv_zeta_interval_bl[,2] - indv_zeta_interval_bl[,1],
                                 zeta_true = indv_zeta_true_bart[i_test],
                                 zeta_catch = as.numeric(indv_zeta_catch_bl),
                                 delta_mean = indv_delta_mean_bl,
                                 delta_lower = indv_delta_interval_bl[,1],
                                 delta_upper = indv_delta_interval_bl[,2],
                                 delta_len = indv_delta_interval_bl[,2] - indv_delta_interval_bl[,1],
                                 delta_true = indv_delta_true_bart[i_test],
                                 delta_catch = as.numeric(indv_delta_catch_bl))
      saveRDS(indv_df_bl_i, file_name_indv_bl)
      
      # Average
      avg_df_bl_i <- data.frame(seed = seeds[i],
                                zeta_mean = avg_zeta_mean_bl,
                                zeta_lower = avg_zeta_interval_bl[1],
                                zeta_upper = avg_zeta_interval_bl[2],
                                zeta_len = avg_zeta_interval_bl[2] - avg_zeta_interval_bl[1],
                                zeta_true = avg_zeta_true_bart,
                                zeta_catch = as.numeric(avg_zeta_catch_bl),
                                delta_mean = avg_delta_mean_bl,
                                delta_lower = avg_delta_interval_bl[1],
                                delta_upper = avg_delta_interval_bl[2],
                                delta_len = avg_delta_interval_bl[2] - avg_delta_interval_bl[1],
                                delta_true = avg_delta_true_bart,
                                delta_catch = as.numeric(avg_delta_catch_bl))
      saveRDS(avg_df_bl_i, file_name_avg_bl)
      
      # Subgroup
      subgroup_df_bl_i <- data.frame(seed = rep(seeds[i], length(subgroup_labels)),
                                     group = subgroup_labels,
                                     zeta_mean = subgroup_zeta_mean_bl,
                                     zeta_lower = subgroup_zeta_interval_bl[,1],
                                     zeta_upper = subgroup_zeta_interval_bl[,2],
                                     zeta_len = subgroup_zeta_interval_bl[,2] - subgroup_zeta_interval_bl[,1],
                                     zeta_true = subgroup_zeta_true_bart[i,2],
                                     zeta_catch = as.numeric(subgroup_zeta_catch_bl),
                                     delta_mean = subgroup_delta_mean_bl,
                                     delta_lower = subgroup_delta_interval_bl[,1],
                                     delta_upper = subgroup_delta_interval_bl[,2],
                                     delta_len = subgroup_delta_interval_bl[,2] - subgroup_delta_interval_bl[,1],
                                     delta_true = subgroup_delta_true_bart[,2],
                                     delta_catch = as.numeric(subgroup_delta_catch_bl))
      saveRDS(subgroup_df_bl_i, file_name_subgroup_bl)
      
      # Save replications into a data frame
      # indv_df_bl[1:nrow(data_test) + nrow(data_test) * (i-1),] <- indv_df_bl_i
      # avg_df_bl[i,] <- avg_df_bl_i
      # subgroup_df_bl[1:5 + 5 * (i-1),] <- subgroup_df_bl_i
    
    }
  }
  
  return(TRUE)
  
  # return(list(
  #   individual_lsem_lsem = indv_df_ll,
  #   average_lsem_lsem = avg_df_ll,
  #   subgroup_lsem_lsem = subgroup_df_ll,
  #   individual_bart_lsem = indv_df_bl,
  #   average_bart_lsem = avg_df_bl,
  #   subgroup_bart_lsem = subgroup_df_bl
  # ))
}