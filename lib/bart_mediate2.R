bart_mediate2 <- function(data_train, data_test, formula_m, formula_y,
                         pi_hat_train, pi_hat_test,
                         m0_hat_train, m0_hat_test,
                         m1_hat_train, m1_hat_test,
                         mediator_name, outcome_name, treat_name,
                         n_iter, burnin, chains = 1) {
  
  ## Process data ----
  X_m_train <- quantile_normalize_bart(
    preprocess_df(
      model.frame(formula_m, data = data_train) %>% 
        select(-all_of(mediator_name))
    )[[1]]
  )
  
  X_y_train <- quantile_normalize_bart(
    preprocess_df(
    cbind(
      model.frame(formula_y, data = data_train) %>% select(-all_of(outcome_name)),
      m0_hat_train, 
      m1_hat_train)
    )[[1]]
  )
  
  X_mu_m_train <- cbind(X_m_train, pi_hat_train)
  X_mu_y_train <- cbind(X_y_train, pi_hat_train)
  
  X_m_test <- quantile_normalize_bart(
    preprocess_df(
      model.frame(formula_m, data = data_test) %>% select(-all_of(mediator_name))
    )[[1]]
  )
  
  X_y_test <- quantile_normalize_bart(
    preprocess_df(
      cbind(model.frame(formula_y, data = data_test) %>% select(-all_of(outcome_name)),
            m0_hat_test, 
            m1_hat_test)
      )[[1]]
    )
  
  X_mu_m_test <- cbind(X_m_test, pi_hat_test)
  X_mu_y_test <- cbind(X_y_test, pi_hat_test)
  
  ## Extracting mediator, outcomes, treatment and scaling ----

  m <- data_train[[mediator_name]]
  m_scale <- (m - mean(m)) / sd(m)

  Y <- data_train[[outcome_name]]
  Y_scale <- (Y - mean(Y)) / sd(Y)
  
  A <- data_train[[treat_name]]
  
  ## Store samples for all chains ----
  
  mu_y_list <- list()
  zeta_list <- list()
  d_list <- list()
  mu_m_list <- list()
  tau_list <- list()
  sigma_y_list <- list()
  sigma_m_list <- list()
  
  for (j in 1:chains) {
    
    bar_format <- paste0(" running chain ", j, " [:bar] :percent eta: :eta")
    
    pb <- progress_bar$new(
      format = bar_format, 
      total = n_iter, 
      clear = FALSE, 
      width= 60
    )

    ## Hypers for mu_m and tau ----
    hypers_tau             <- Hypers(X_m_train, m_scale, normalize_Y = FALSE)
    hypers_mu_m            <- hypers_tau
    opts_mu_m              <- Opts(update_s = FALSE, update_sigma_mu = FALSE)
    opts_tau               <- opts_mu_m
    opts_tau$update_sigma  <- FALSE
  
    ## Hypers for mu_y, zeta, d ----
    hypers_mu_y            <- Hypers(X_y_train, Y_scale, normalize_Y = FALSE)
    hypers_zeta            <- hypers_mu_y
    hypers_d               <- hypers_mu_y
    opts_mu_y              <- Opts(update_s = FALSE, update_sigma_mu = FALSE)
    opts_zeta              <- opts_d <- opts_mu_y
    opts_zeta$update_sigma <- FALSE
    opts_d$update_sigma    <- FALSE

    # Forests
    forest_mu_m <- MakeForest(hypers_mu_m, opts_mu_m)
    forest_tau  <- MakeForest(hypers_tau, opts_tau)
    forest_mu_y <- MakeForest(hypers_mu_y, opts_mu_y)
    forest_zeta <- MakeForest(hypers_zeta, opts_zeta)
    forest_d    <- MakeForest(hypers_d, opts_d)
  
    # Initialize mu_m and tau
    mu_m <- forest_mu_m$do_predict(X_mu_m_train)
    tau <- forest_tau$do_predict(X_m_train)

    # Initialize mu_y, zeta, d
    mu_y <- forest_mu_y$do_predict(X_mu_y_train)
    zeta <- forest_zeta$do_predict(X_y_train)
    d <- forest_d$do_predict(X_y_train)

    # Store samples for one chain
    mu_y_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
    zeta_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
    d_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
    mu_m_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
    tau_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
    sigma_y_samples <- rep(NA, n_iter - burnin)
    sigma_m_samples <- rep(NA, n_iter - burnin)

    for (i in 1:n_iter) {
      
      pb$tick()

      # Mi(a) ----
      
      mu_m  <- forest_mu_m$do_gibbs(X_mu_m_train, m_scale - A * tau, X_mu_m_train, 1)
      mu_m_test <- forest_mu_m$do_predict(X_mu_m_test)
      sigma_m <- forest_mu_m$get_sigma()
      
      forest_tau$set_sigma(sigma_m)
      tau <- forest_tau$do_gibbs(X_m_train[A == 1,], m_scale[A == 1] - mu_m[A == 1], X_m_train, 1)
      tau_test <- forest_tau$do_predict(X_m_test)

      # Yi(a) ----
      
      mu_y  <- forest_mu_y$do_gibbs(X_mu_y_train, Y_scale - A * zeta - m_scale * d, X_mu_y_train, 1)
      mu_y_test <- forest_mu_y$do_predict(X_mu_y_test)
      sigma_y <- forest_mu_y$get_sigma()

      forest_zeta$set_sigma(sigma_y)
      zeta  <- forest_zeta$do_gibbs(
        X_y_train[A == 1,], 
        Y_scale[A == 1] - mu_y[A == 1] - (m_scale * d)[A == 1], 
        X_y_train, 
      1)
      zeta_test <- forest_zeta$do_predict(X_y_test)
      
      R <- (Y_scale - mu_y - A * zeta) / m_scale
      forest_d$set_sigma(sigma_y)
      d <- forest_d$do_gibbs_weighted(
        X_y_train, 
        R, 
        m_scale^2, 
        X_y_train, 
      1)
      d_test <- forest_d$do_predict(X_y_test)

      # Unscale
      mu_y_unscale <- mean(Y) + (mu_y_test * sd(Y)) - mean(m) / sd(m) * d_test * sd(Y)
      zeta_unscale <- zeta_test * sd(Y)
      d_unscale <- d_test * sd(Y) / sd(m)
      mu_m_unscale <- mean(m) + (mu_m_test * sd(m))
      tau_unscale <- tau_test * sd(m)
      sigma_y_unscale <- sigma_y * sd(Y)
      sigma_m_unscale <- sigma_m * sd(m)

      if (i > burnin){
        mu_y_samples[i - burnin,] <- mu_y_unscale
        zeta_samples[i - burnin,] <- zeta_unscale
        d_samples[i - burnin,] <- d_unscale
        mu_m_samples[i - burnin,] <- mu_m_unscale
        tau_samples[i - burnin,] <- tau_unscale
        sigma_y_samples[i - burnin] <- sigma_y_unscale
        sigma_m_samples[i - burnin] <- sigma_m_unscale
      }
    }

    mu_y_list[[j]] <- mu_y_samples
    zeta_list[[j]] <- zeta_samples
    d_list[[j]] <- d_samples
    mu_m_list[[j]] <- mu_m_samples
    tau_list[[j]] <- tau_samples
    sigma_y_list[[j]] <- sigma_y_samples
    sigma_m_list[[j]] <- sigma_m_samples
    
  }

  return(list(
    mu_y_samples = do.call(rbind, mu_y_list),
    zeta_samples = do.call(rbind, zeta_list),
    d_samples = do.call(rbind, d_list),
    mu_m_samples = do.call(rbind, mu_m_list),
    tau_samples = do.call(rbind, tau_list),
    sigma_y_samples = unlist(sigma_y_list),
    sigma_m_samples = unlist(sigma_m_list), 
    chain = factor(rep(1:chains, each = n_iter - burnin))
  ))
}