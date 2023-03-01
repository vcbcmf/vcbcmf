do_simulation_bart <- function(data, i_train, i_test,
                               formula_m_bart, formula_y_bart, formula_ps_bart,
                               fit_m_lsem, fit_y_lsem, formula_X_lsem,
                               mediator_name, outcome_name, treat_name,
                               mu_y_hat_bart, zeta_hat_bart, d_hat_bart,
                               mu_m_hat_bart, tau_hat_bart,
                               sigma_y_hat_bart, sigma_m_hat_bart,
                               zeta_hat_lsem, delta_hat_lsem,
                               sigma_y_hat_lsem, sigma_m_hat_lsem,
                               n_iter, burnin, thin, seeds) {

  
  delta <- delta_hat <- delta_lower <- delta_true <- delta_upper <- leaf <- 
    zeta <- zeta_hat <- zeta_lower <- zeta_true <- zeta_upper <- iteration <- 
    NULL
  
  # Get training and testing set
  data_train <- data[i_train,]
  data_test <- data[i_test,]

  # Store results
  colnames_indv      <- c('seed', 'subj_id', 'zeta_mean', 'zeta_lower',
                          'zeta_upper', 'zeta_len', 'zeta_true', 'zeta_catch',
                          'delta_mean', 'delta_lower', 'delta_upper',
                          'delta_len', 'delta_true', 'delta_catch')
  colnames_subgroup <- c('seed', 'group', 'zeta_mean', 'zeta_lower',
                         'zeta_upper', 'zeta_len', 'zeta_true', 'zeta_catch',
                         'delta_mean', 'delta_lower', 'delta_upper',
                         'delta_len', 'delta_true', 'delta_catch')
  colnames_avg       <- c('seed', 'zeta_mean', 'zeta_lower', 'zeta_upper',
                          'zeta_len', 'zeta_true', 'zeta_catch',
                          'delta_mean', 'delta_lower', 'delta_upper',
                          'delta_len', 'delta_true', 'delta_catch')
  colnames_tree_subgroups <- c('seed', 'group', 'zeta_mean', 'zeta_lower',
                               'zeta_upper', 'zeta_len', 'zeta_true',
                               'zeta_catch', 'delta_mean', 'delta_lower',
                               'delta_upper', 'delta_len', 'delta_true',
                               'delta_catch')
  
  # # BART/BART
  # indv_df_bb      <- data.frame(matrix(NA, nrow = n_reps * nrow(data_test),
  #                                      ncol = length(colnames_indv),
  #                                      dimnames = list(c(), colnames_indv)))
  # subgroup_df_bb <- data.frame(matrix(NA, nrow = n_reps * 5,
  #                                     ncol = length(colnames_subgroup),
  #                                     dimnames = list(c(), colnames_subgroup)))
  # avg_df_bb       <- data.frame(matrix(NA, nrow = n_reps,
  #                                      ncol = length(colnames_avg),
  #                                      dimnames = list(c(), colnames_avg)))
  # tree_subgroups_list_bb <- list()
  #
  # # LSEM/BART
  # indv_df_lb      <- data.frame(matrix(NA, nrow = n_reps * nrow(data_test),
  #                                      ncol = length(colnames_indv),
  #                                      dimnames = list(c(), colnames_indv)))
  # subgroup_df_lb <- data.frame(matrix(NA, nrow = n_reps * 5,
  #                                     ncol = length(colnames_subgroup),
  #                                     dimnames = list(c(), colnames_subgroup)))
  # avg_df_lb       <- data.frame(matrix(NA, nrow = n_reps,
  #                                      ncol = length(colnames_avg),
  #                                      dimnames = list(c(), colnames_avg)))
  
  n_reps <- length(seeds)
  
  for (i in 1:n_reps) {

    file_name_indv_bb <- paste0('Simulation BARTBART/indv_seed', seeds[i], '.rds')
    file_name_subgroup_bb <- paste0('Simulation BARTBART/subgroup_seed', seeds[i], '.rds')
    file_name_avg_bb <- paste0('Simulation BARTBART/avg_seed', seeds[i], '.rds')
    file_name_tree_subgroup_bb <- paste0('Simulation BARTBART/tree_subgroup_seed', seeds[i], '.rds')

    file_name_indv_lb <- paste0('Simulation LSEMBART/indv_seed', seeds[i], '.rds')
    file_name_subgroup_lb <- paste0('Simulation LSEMBART/subgroup_seed', seeds[i], '.rds')
    file_name_avg_lb <- paste0('Simulation LSEMBART/avg_seed', seeds[i], '.rds')

    if (file.exists(file_name_indv_bb) & file.exists(file_name_subgroup_bb) &
        file.exists(file_name_avg_bb) & file.exists(file_name_tree_subgroup_bb) &
        file.exists(file_name_indv_lb) & file.exists(file_name_subgroup_lb) &
        file.exists(file_name_avg_lb)) {
      # if (file.exists(file_name_indv_lb) & file.exists(file_name_avg_lb)) {
      indv_df_bb_i <- readRDS(file_name_indv_bb)
      subgroup_df_bb_i <- readRDS(file_name_subgroup_bb)
      avg_df_bb_i <- readRDS(file_name_avg_bb)
      tree_subgroups_list_bb_i <- readRDS(file_name_tree_subgroup_bb)

      indv_df_lb_i <- readRDS(file_name_indv_lb)
      subgroup_df_lb_i <- readRDS(file_name_subgroup_lb)
      avg_df_lb_i <- readRDS(file_name_avg_lb)
    }
    else {
      set.seed(seeds[i])

      # Get true effects for BART
      indv_zeta_true_bart   <- zeta_hat_bart
      indv_delta_true_bart  <- d_hat_bart * tau_hat_bart
      avg_zeta_true_bart    <- mean(indv_zeta_true_bart)
      avg_delta_true_bart   <- mean(indv_delta_true_bart)

      subgroup_labels <- c('white, age < 34', 'non-white, age < 34', 'age ≥ 67',
                           'white, 34 ≤ age < 67', 'non-white, 34 ≤ age < 67')
      i1 <- which(data$age < 34 & data$race == 'White')
      i2 <- which(data$age < 34 & data$race != 'White')
      i3 <- which(data$age >= 67)
      i4 <- which(data$age >= 34 & data$age < 67 & data$race == 'White')
      i5 <- which(data$age >= 34 & data$age < 67 & data$race != 'White')

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

      # Get true effects for LSEM
      indv_zeta_true_lsem  <- zeta_hat_lsem
      indv_delta_true_lsem <- delta_hat_lsem
      avg_zeta_true_lsem   <- mean(indv_zeta_true_lsem)
      avg_delta_true_lsem  <- mean(indv_delta_true_lsem)

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

      # Make training simulation set for LSEM truth
      epsilon_m_lsem <- rnorm(nrow(data_train))
      epsilon_y_lsem <- rnorm(nrow(data_train))

      m_sim_lsem <- predict(fit_m_lsem, newdata = data_train) + sigma_m_hat_lsem * epsilon_m_lsem
      Y_sim_lsem <- predict(fit_y_lsem, newdata = data_train %>% mutate(!!mediator_name := m_sim_lsem)) + sigma_y_hat_lsem * epsilon_y_lsem

      data_train_sim_lsem <- data_train
      data_train_sim_lsem[[mediator_name]] <- m_sim_lsem
      data_train_sim_lsem[[outcome_name]] <- Y_sim_lsem

      ## BART Truth / BART Fit ----
      # Clever covariates for training/testing set
      clever_cov_bb <- get_clever_cov(data_train_sim_bart, data, formula_m_bart,
                                      mediator_name, outcome_name, treat_name)
      m0_hat_bb <- clever_cov_bb$m0_hat
      m1_hat_bb <- clever_cov_bb$m1_hat
      m0_hat_train_bb <- m0_hat_bb[i_train]
      m1_hat_train_bb <- m1_hat_bb[i_train]
      rm(clever_cov_bb)
      gc()

      # Propensity score for training/testing set
      pi_hat_bb <- get_ps(data, data, formula_ps_bart)
      pi_hat_train_bb <- pi_hat_bb[i_train]

      # Get output for simulated dataset
      out_bb <- bart_mediate(data_train_sim_bart, data,
                             formula_m_bart, formula_y_bart,
                             pi_hat_train_bb, pi_hat_bb,
                             m0_hat_train_bb, m0_hat_bb,
                             m1_hat_train_bb, m1_hat_bb,
                             'phealth', 'logY', 'smoke',
                             n_iter, burnin, thin)

      # Get simulated direct/indirect distributions
      zeta_bb <- out_bb$zeta_samples
      zeta_train_bb <- zeta_bb[,i_train]
      # zeta_test_bb <- zeta_bb[,i_test]

      # d_bb <- out_bb$d_samples
      # d_train_bb <- d_bb[,i_train]
      # d_test_bb <- d_bb[,i_test]

      # tau_bb <- out_bb$tau_samples
      # tau_train_bb <- tau_bb[,i_train]
      # tau_test_bb <- tau_bb[,i_test]

      delta_bb <- out_bb$d_samples * out_bb$tau_samples
      delta_train_bb <- delta_bb[,i_train]
      # delta_test_bb <- delta_bb[,i_test]

      rm(out_bb)
      gc()

      avg_zeta_bb    <- rowMeans(zeta_train_bb)
      avg_delta_bb  <- rowMeans(delta_train_bb)
      indv_zeta_bb   <- zeta_bb[,i_test]
      indv_delta_bb <- delta_bb[,i_test]

      i1_train <- which(data_train$age < 34 & data_train$race == 'White')
      i2_train <- which(data_train$age < 34 & data_train$race != 'White')
      i3_train <- which(data_train$age >= 67)
      i4_train <- which(data_train$age >= 34 & data_train$age < 67 & data_train$race == 'White')
      i5_train <- which(data_train$age >= 34 & data_train$age < 67 & data_train$race != 'White')

      group1_zeta_bb <- rowMeans(zeta_train_bb[,i1_train])
      group2_zeta_bb <- rowMeans(zeta_train_bb[,i2_train])
      group3_zeta_bb <- rowMeans(zeta_train_bb[,i3_train])
      group4_zeta_bb <- rowMeans(zeta_train_bb[,i4_train])
      group5_zeta_bb <- rowMeans(zeta_train_bb[,i5_train])
      subgroup_zeta_bb <- cbind(group1_zeta_bb, group2_zeta_bb, group3_zeta_bb,
                                group4_zeta_bb, group5_zeta_bb)

      group1_delta_bb <- rowMeans(delta_train_bb[,i1_train])
      group2_delta_bb <- rowMeans(delta_train_bb[,i2_train])
      group3_delta_bb <- rowMeans(delta_train_bb[,i3_train])
      group4_delta_bb <- rowMeans(delta_train_bb[,i4_train])
      group5_delta_bb <- rowMeans(delta_train_bb[,i5_train])
      subgroup_delta_bb <- cbind(group1_delta_bb, group2_delta_bb, group3_delta_bb,
                                 group4_delta_bb, group5_delta_bb)

      # Get direct/indirect credible intervals
      avg_zeta_interval_bb    <- quantile(avg_zeta_bb, probs = c(0.025, 0.975))
      avg_delta_interval_bb   <- quantile(avg_delta_bb, probs = c(0.025, 0.975))

      indv_zeta_interval_bb   <- colQuantiles(indv_zeta_bb, probs = c(0.025, 0.975))
      indv_delta_interval_bb  <- colQuantiles(indv_delta_bb, probs = c(0.025, 0.975))

      subgroup_zeta_interval_bb <- colQuantiles(subgroup_zeta_bb, probs = c(0.025, 0.975))
      subgroup_delta_interval_bb <- colQuantiles(subgroup_delta_bb, probs = c(0.025, 0.975))

      # Mean point estimates
      avg_zeta_mean_bb    <- mean(avg_zeta_bb)
      avg_delta_mean_bb  <- mean(avg_delta_bb)

      indv_zeta_mean_bb   <- colMeans(indv_zeta_bb)
      indv_delta_mean_bb <- colMeans(indv_delta_bb)

      subgroup_zeta_mean_bb <- colMeans(subgroup_zeta_bb)
      subgroup_delta_mean_bb <- colMeans(subgroup_delta_bb)

      # Catch T/F
      avg_zeta_catch_bb <- (avg_zeta_true_bart >= avg_zeta_interval_bb[1]) &
        (avg_zeta_true_bart <= avg_zeta_interval_bb[2])
      avg_delta_catch_bb <- (avg_delta_true_bart >= avg_delta_interval_bb[1]) &
        (avg_delta_true_bart <= avg_delta_interval_bb[2])

      indv_zeta_catch_bb <- unlist(lapply(
        X = 1:length(indv_zeta_true_bart[i_test]),
        FUN = function(i) (indv_zeta_true_bart[i_test][i] >= indv_zeta_interval_bb[i,1]) &
                          (indv_zeta_true_bart[i_test][i] <= indv_zeta_interval_bb[i,2])
      ))
      indv_delta_catch_bb <- unlist(lapply(
        X = 1:length(indv_delta_true_bart[i_test]),
        FUN = function(i) (indv_delta_true_bart[i_test][i] >= indv_delta_interval_bb[i,1]) &
                          (indv_delta_true_bart[i_test][i] <= indv_delta_interval_bb[i,2])
      ))

      subgroup_zeta_catch_bb <- unlist(lapply(
        X = 1:nrow(subgroup_zeta_true_bart),
        FUN = function(i) (subgroup_zeta_true_bart[i,2] >= subgroup_zeta_interval_bb[i,1]) &
                          (subgroup_zeta_true_bart[i,2] <= subgroup_zeta_interval_bb[i,2])
      ))
      subgroup_delta_catch_bb <- unlist(lapply(
        X = 1:nrow(subgroup_delta_true_bart),
        FUN = function(i) (subgroup_delta_true_bart[i,2] >= subgroup_delta_interval_bb[i,1]) &
                          (subgroup_delta_true_bart[i,2] <= subgroup_delta_interval_bb[i,2])
      ))


      # Saving
      # Individual
      indv_df_bb_i <- data.frame(
        seed = rep(seeds[i], nrow(data_test)),
        subj_id = 1:nrow(data_test),
        zeta_mean = indv_zeta_mean_bb,
        zeta_lower = indv_zeta_interval_bb[, 1],
        zeta_upper = indv_zeta_interval_bb[, 2],
        zeta_len = indv_zeta_interval_bb[, 2] - indv_zeta_interval_bb[, 1],
        zeta_true = indv_zeta_true_bart[i_test],
        zeta_catch = as.numeric(indv_zeta_catch_bb),
        delta_mean = indv_delta_mean_bb,
        delta_lower = indv_delta_interval_bb[, 1],
        delta_upper = indv_delta_interval_bb[, 2],
        delta_len = indv_delta_interval_bb[, 2] - indv_delta_interval_bb[, 1],
        delta_true = indv_delta_true_bart[i_test],
        delta_catch = as.numeric(indv_delta_catch_bb)
      )
      saveRDS(indv_df_bb_i, file_name_indv_bb)

      # Subgroup
      subgroup_df_bb_i <- data.frame(
        seed = rep(seeds[i], length(subgroup_labels)),
        group = subgroup_labels,
        zeta_mean = subgroup_zeta_mean_bb,
        zeta_lower = subgroup_zeta_interval_bb[, 1],
        zeta_upper = subgroup_zeta_interval_bb[, 2],
        zeta_len = subgroup_zeta_interval_bb[, 2] - subgroup_zeta_interval_bb[, 1],
        zeta_true = subgroup_zeta_true_bart[, 2],
        zeta_catch = as.numeric(subgroup_zeta_catch_bb),
        delta_mean = subgroup_delta_mean_bb,
        delta_lower = subgroup_delta_interval_bb[, 1],
        delta_upper = subgroup_delta_interval_bb[, 2],
        delta_len = subgroup_delta_interval_bb[, 2] - subgroup_delta_interval_bb[, 1],
        delta_true = subgroup_delta_true_bart[, 2],
        delta_catch = as.numeric(subgroup_delta_catch_bb)
      )
      saveRDS(subgroup_df_bb_i, file_name_subgroup_bb)

      # Average
      avg_df_bb_i <- data.frame(
        seed = seeds[i],
        zeta_mean = avg_zeta_mean_bb,
        zeta_lower = avg_zeta_interval_bb[1],
        zeta_upper = avg_zeta_interval_bb[2],
        zeta_len = avg_zeta_interval_bb[2] - avg_zeta_interval_bb[1],
        zeta_true = avg_zeta_true_bart,
        zeta_catch = as.numeric(avg_zeta_catch_bb),
        delta_mean = avg_delta_mean_bb,
        delta_lower = avg_delta_interval_bb[1],
        delta_upper = avg_delta_interval_bb[2],
        delta_len = avg_delta_interval_bb[2] - avg_delta_interval_bb[1],
        delta_true = avg_delta_true_bart,
        delta_catch = as.numeric(avg_delta_catch_bb)
      )
      saveRDS(avg_df_bb_i, file_name_avg_bb)

      # Tree Projection Subgroups
      data_post <- model.frame(formula_m_bart, data = data) %>%
        select(-mediator_name) %>%
        mutate(zeta_hat = colMeans(zeta_bb), delta_hat = colMeans(delta_bb))

      tree <- rpart(delta_hat ~ . - zeta_hat, data = data_post)

      # True subroup values
      true_leaf_df <- data.frame(zeta_true = indv_zeta_true_bart,
                                 delta_true = indv_delta_true_bart,
                                 leaf = tree$where)

      true_leaf_means <- true_leaf_df %>% group_by(leaf) %>%
        summarise(zeta_true = mean(zeta_true),
                  delta_true = mean(delta_true))

      # Simulated subgroup values
      sim_leaf_df <- data.frame(zeta = c(t(zeta_bb)),
                                delta = c(t(delta_bb)),
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

      tree_subgroups_list_bb_i <- data.frame(seed = seeds[i],
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
      saveRDS(tree_subgroups_list_bb_i, file_name_tree_subgroup_bb)
      rm(data_post, true_leaf_df, sim_leaf_df, sim_leaf_means_df, sim_leaf_means_df)
      gc()

      ## Clean ----
      
      rm(zeta_bb, zeta_train_bb, delta_bb, delta_train_bb, 
         indv_zeta_bb, indv_delta_bb)
      
      gc()
      
      ## LSEM Truth / BART Fit ----
      clever_cov_lb <- get_clever_cov(data_train_sim_lsem,
                                      data,
                                      formula_m_bart,
                                      mediator_name,
                                      outcome_name,
                                      treat_name)
      m0_hat_lb <- clever_cov_lb$m0_hat
      m1_hat_lb <- clever_cov_lb$m1_hat
      m0_hat_train_lb <- m0_hat_lb[i_train]
      m1_hat_train_lb <- m1_hat_lb[i_train]
      rm(clever_cov_lb)
      gc()

      # Propensity score for training/testing set
      pi_hat_lb <- get_ps(data, data, formula_ps_bart)
      pi_hat_train_lb <- pi_hat_lb[i_train]

      # Get output for simulated dataset
      out_lb <- bart_mediate(data_train_sim_lsem, data,
                             formula_m_bart, formula_y_bart,
                             pi_hat_train_lb, pi_hat_lb,
                             m0_hat_train_lb, m0_hat_lb,
                             m1_hat_train_lb, m1_hat_lb,
                             'phealth', 'logY', 'smoke',
                             n_iter, burnin, thin)

      # Get simulated direct/indirect distributions
      zeta_lb <- out_lb$zeta_samples
      zeta_train_lb <- zeta_lb[,i_train]
      # zeta_lb[,i_test]

      # d_lb <- out_lb$d_samples
      # d_train_lb <- d_lb[,i_train]
      # d_test_lb <- d_lb[,i_test]

      # tau_lb <- out_lb$tau_samples
      # tau_train_lb <- tau_lb[,i_train]
      # tau_test_lb <- tau_lb[,i_test]

      delta_lb <- out_lb$d_samples * out_lb$tau_samples
      delta_train_lb <- delta_lb[,i_train]
      # delta_test_lb <- delta_lb[,i_test]

      rm(out_lb)
      gc()

      avg_zeta_lb    <- rowMeans(zeta_train_lb)
      avg_delta_lb  <- rowMeans(delta_train_lb)
      indv_zeta_lb   <- zeta_lb[,i_test]
      indv_delta_lb <- delta_lb[,i_test]

      group1_zeta_lb <- rowMeans(zeta_train_lb[,i1_train])
      group2_zeta_lb <- rowMeans(zeta_train_lb[,i2_train])
      group3_zeta_lb <- rowMeans(zeta_train_lb[,i3_train])
      group4_zeta_lb <- rowMeans(zeta_train_lb[,i4_train])
      group5_zeta_lb <- rowMeans(zeta_train_lb[,i5_train])
      subgroup_zeta_lb <- cbind(group1_zeta_lb, group2_zeta_lb, group3_zeta_lb,
                                group4_zeta_lb, group5_zeta_lb)

      group1_delta_lb <- rowMeans((delta_train_lb)[,i1_train])
      group2_delta_lb <- rowMeans((delta_train_lb)[,i2_train])
      group3_delta_lb <- rowMeans((delta_train_lb)[,i3_train])
      group4_delta_lb <- rowMeans((delta_train_lb)[,i4_train])
      group5_delta_lb <- rowMeans((delta_train_lb)[,i5_train])
      subgroup_delta_lb <- cbind(group1_delta_lb, group2_delta_lb, group3_delta_lb,
                                 group4_delta_lb, group5_delta_lb)

      # Get direct/indirect credible intervals
      avg_zeta_interval_lb    <- quantile(avg_zeta_lb, probs = c(0.025, 0.975))
      avg_delta_interval_lb   <- quantile(avg_delta_lb, probs = c(0.025, 0.975))

      indv_zeta_interval_lb   <- colQuantiles(indv_zeta_lb, probs = c(0.025, 0.975))
      indv_delta_interval_lb  <- colQuantiles(indv_delta_lb, probs = c(0.025, 0.975))

      subgroup_zeta_interval_lb <- colQuantiles(subgroup_zeta_lb, probs = c(0.025, 0.975))
      subgroup_delta_interval_lb <- colQuantiles(subgroup_delta_lb, probs = c(0.025, 0.975))

      # Mean point estimates
      avg_zeta_mean_lb    <- mean(avg_zeta_lb)
      avg_delta_mean_lb  <- mean(avg_delta_lb)

      indv_zeta_mean_lb   <- colMeans(indv_zeta_lb)
      indv_delta_mean_lb <- colMeans(indv_delta_lb)

      subgroup_zeta_mean_lb <- colMeans(subgroup_zeta_lb)
      subgroup_delta_mean_lb <- colMeans(subgroup_delta_lb)

      # Catch T/F
      avg_zeta_catch_lb <- (avg_zeta_true_lsem >= avg_zeta_interval_lb[1]) &
        (avg_zeta_true_lsem <= avg_zeta_interval_lb[2])
      avg_delta_catch_lb <- (avg_delta_true_lsem >= avg_delta_interval_lb[1]) &
        (avg_delta_true_lsem <= avg_delta_interval_lb[2])

      indv_zeta_catch_lb <- unlist(lapply(1:length(indv_zeta_true_lsem[i_test]),
                                          function(i) (indv_zeta_true_lsem[i_test][i] >= indv_zeta_interval_lb[i,1]) &
                                            (indv_zeta_true_lsem[i_test][i] <= indv_zeta_interval_lb[i,2])))
      indv_delta_catch_lb <- unlist(lapply(1:length(indv_delta_true_lsem[i_test]),
                                           function(i) (indv_delta_true_lsem[i_test][i] >= indv_delta_interval_lb[i,1]) &
                                             (indv_delta_true_lsem[i_test][i] <= indv_delta_interval_lb[i,2])))

      subgroup_zeta_catch_lb <- unlist(lapply(1:nrow(subgroup_zeta_true_lsem),
                                              function(i) (subgroup_zeta_true_lsem[i,2] >= subgroup_zeta_interval_lb[i,1]) &
                                                (subgroup_zeta_true_lsem[i,2] <= subgroup_zeta_interval_lb[i,2])))
      subgroup_delta_catch_lb <- unlist(lapply(1:nrow(subgroup_delta_true_lsem),
                                               function(i) (subgroup_delta_true_lsem[i,2] >= subgroup_delta_interval_lb[i,1]) &
                                                 (subgroup_delta_true_lsem[i,2] <= subgroup_delta_interval_lb[i,2])))

      # Saving
      # Individual
      indv_df_lb_i <- data.frame(seed = rep(seeds[i], nrow(data_test)),
                                 subj_id = 1:nrow(data_test),
                                 zeta_mean = indv_zeta_mean_lb,
                                 zeta_lower = indv_zeta_interval_lb[,1],
                                 zeta_upper = indv_zeta_interval_lb[,2],
                                 zeta_len = indv_zeta_interval_lb[,2] - indv_zeta_interval_lb[,1],
                                 zeta_true = indv_zeta_true_lsem[i_test],
                                 zeta_catch = as.numeric(indv_zeta_catch_lb),
                                 delta_mean = indv_delta_mean_lb,
                                 delta_lower = indv_delta_interval_lb[,1],
                                 delta_upper = indv_delta_interval_lb[,2],
                                 delta_len = indv_delta_interval_lb[,2] - indv_delta_interval_lb[,1],
                                 delta_true = indv_delta_true_lsem[i_test],
                                 delta_catch = as.numeric(indv_delta_catch_lb))
      saveRDS(indv_df_lb_i, file_name_indv_lb)

      # Subgroup
      subgroup_df_lb_i <- data.frame(seed = rep(seeds[i], length(subgroup_labels)),
                                     group = subgroup_labels,
                                     zeta_mean = subgroup_zeta_mean_lb,
                                     zeta_lower = subgroup_zeta_interval_lb[,1],
                                     zeta_upper = subgroup_zeta_interval_lb[,2],
                                     zeta_len = subgroup_zeta_interval_lb[,2] - subgroup_zeta_interval_lb[,1],
                                     zeta_true = subgroup_zeta_true_lsem[,2],
                                     zeta_catch = as.numeric(subgroup_zeta_catch_lb),
                                     delta_mean = subgroup_delta_mean_lb,
                                     delta_lower = subgroup_delta_interval_lb[,1],
                                     delta_upper = subgroup_delta_interval_lb[,2],
                                     delta_len = subgroup_delta_interval_lb[,2] - subgroup_delta_interval_lb[,1],
                                     delta_true = subgroup_delta_true_lsem[,2],
                                     delta_catch = as.numeric(subgroup_delta_catch_lb))
      saveRDS(subgroup_df_lb_i, file_name_subgroup_lb)

      # Average
      avg_df_lb_i <- data.frame(seed = seeds[i],
                                zeta_mean = avg_zeta_mean_lb,
                                zeta_lower = avg_zeta_interval_lb[1],
                                zeta_upper = avg_zeta_interval_lb[2],
                                zeta_len = avg_zeta_interval_lb[2] - avg_zeta_interval_lb[1],
                                zeta_true = avg_zeta_true_lsem,
                                zeta_catch = as.numeric(avg_zeta_catch_lb),
                                delta_mean = avg_delta_mean_lb,
                                delta_lower = avg_delta_interval_lb[1],
                                delta_upper = avg_delta_interval_lb[2],
                                delta_len = avg_delta_interval_lb[2] - avg_delta_interval_lb[1],
                                delta_true = avg_delta_true_lsem,
                                delta_catch = as.numeric(avg_delta_catch_lb))
      saveRDS(avg_df_lb_i, file_name_avg_lb)
      
      ## Clean ----
      
      rm(zeta_lb, zeta_train_lb, delta_lb, delta_train_lb, indv_zeta_lb, indv_delta_lb)
      gc()
      
    }

    # # Save replications into a data frame
    # indv_df_bb[1:nrow(data_test) + nrow(data_test) * (i-1),] <- indv_df_bb_i
    # subgroup_df_bb[1:5 + 5 * (i-1),] <- subgroup_df_bb_i
    # avg_df_bb[i,] <- avg_df_bb_i
    # tree_subgroups_list_bb[[i]] <- tree_subgroups_list_bb_i
    #
    # indv_df_lb[1:nrow(data_test) + nrow(data_test) * (i-1),] <- indv_df_lb_i
    # subgroup_df_lb[1:5 + 5 * (i-1),] <- subgroup_df_lb_i
    # avg_df_lb[i,] <- avg_df_lb_i
    rm(indv_df_bb_i, subgroup_df_bb_i, avg_df_bb_i, tree_subgroups_list_bb_i,
       indv_df_lb_i, subgroup_df_lb_i, avg_df_lb_i)
    gc()
  }

  # tree_subgroups_df_bb <- do.call(rbind, tree_subgroups_list_bb)

  # return(list(
  #   individual_bart_bart = indv_df_bb,
  #   subgroups_bart_bart = subgroup_df_bb,
  #   average_bart_bart = avg_df_bb,
  #   tree_subgroups_bart_bart = tree_subgroups_df_bb,
  #   individual_lsem_bart = indv_df_lb,
  #   subgroups_lsem_bart = subgroup_df_lb,
  #   average_lsem_bart = avg_df_lb
  # ))
}