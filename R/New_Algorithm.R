## ----message = FALSE, purl = TRUE-----------------------------------------------------------------------------------------------------------
library(tibble)
library(ggplot2)
library(magrittr)
library(dplyr) 
library(patchwork)
# loads the functions from the package




## ----analysis helpers, purl = TRUE----------------------------------------------------------------------------------------------------------
noise_fn <- function(x, sd){
  x+rnorm(1, mean = 0, sd = sd)
}



## ----Helper Functions, purl = TRUE----------------------------------------------------------------------------------------------------------
# function inside of parameter distance
lim_func <- function(k_1, k_2, tau_1, tau_2) {
  k_1*exp(-1/tau_1)/(1-exp(-1/tau_1))-k_2*exp(-1/tau_2)/(1-exp(-1/tau_2))
}

# computes "parameter distance"
params_dist <- function(params_1, params_2, lambda_1) {
  # ignore params_1[[1]] since it is just the inital value
  lambda_1 * abs(lim_func(params_1[[2]], params_1[[3]], params_1[[4]], params_1[[5]]) -
        lim_func(params_2[[2]], params_2[[3]], params_2[[4]], params_2[[5]])
      )
}

params_local_diff <- function(params, N, lambda_2) {
  coef_1 <- exp(-1/params[[4]]); coef_2 <- exp(-1/params[[5]])
  exp_sum_1 <- 0; exp_sum_2 <- 0
  for (i in 1:N) {
   exp_sum_1 = exp_sum_1 + coef_1^i
   exp_sum_2 = exp_sum_2 + coef_2^i
  }
  return(lambda_2*(params[[2]]*exp_sum_1 - params[[3]]*exp_sum_2))
}

alg_param_cost_func <- function(x, y, params_settings) {
  training_level <- params_settings$training_level
  lambda_1 <- params_settings$lambda_1
  N <- params_settings$N
  training_level * (params_dist(x, y, lambda_1) +
                      abs(
                        params_local_diff(x, N, lambda_2) -
                          params_local_diff(y, N, lambda_2)
                      ))
}

alg_TE_cost_func <-
  function(P_t,
           N_t,
           tau_1_i,
           tau_2_i,
           tau_1_prev,
           tau_2_prev,
           params_settings) {
    N <- params_settings$N
    lambda_2 <- params_settings$lambda_2
    return(lambda_2 * abs((
      exp(-N / tau_1_i) * P_t -
        exp(-N / tau_2_i) * N_t
    ) -
      (
        exp(-N / tau_1_prev) * P_t -
          exp(-N / tau_2_prev) * N_t
      )))
  }


local_cost <- function(params_i, params_old_i, N, tau_i_sum_mat, lambda) {
  cost = 0
  for (j in 1:N) {
    cost = cost + lambda*(
      abs(params_i[[2]]*tau_i_sum_mat[params_i[[6]], j] - 
                        params_old_i[[2]]*tau_i_sum_mat[params_old_i[[6]], j]) +
      abs(params_i[[3]]*tau_i_sum_mat[params_i[[7]], j] - 
                        params_old_i[[3]]*tau_i_sum_mat[params_old_i[[7]], j])
    )
  }
  return(cost)
}


alg_param_cost_func_2 <- function(params_i, 
                                  min_params_mat,
                                  params_settings,
                                  tau_i_sum_mat) {
  n_obs <- nrow(min_params_mat)
  output_vec <- c(rep(0, n_obs))
  N <- params_settings$N 
  lambda <- params_settings$lambda
  for (i in 1:n_obs) {
    output_vec[[i]] = params_dist(params_i, min_params_mat[i, ], lambda) +
      local_cost(params_i, min_params_mat[i, ], N, tau_i_sum_mat, lambda)
  }
  return(output_vec)
}


# plotting function
plot_perf <- function(observed_performance,
                      modeled_performance) {
  data <- tibble(
  "day" = c(0:length(observed_performance)),
  "pred" = modeled_performance,
  "obs" = c(NA,observed_performance)
  )
  
  plot <- ggplot(data, aes(x = day)) + 
  geom_point(aes(y = pred, color = "pred")) + 
  geom_point(aes(y = obs, color = "obs")) +
  labs(x = "Day",
       y = "Performance",
       title = "") 
  plot
}

## ----update parameters functions, purl = TRUE-----------------------------------------------------------------------------------------------
update_cost <- function(params_i, prev_params, # these two are the only things that change each iteration
                        obs_indexes, training_load, obs_perf,
                        error_func, cost_error_func, param_cost_func, 
                        TE_cost_func, 
                        alpha, window, 
                        min_cost_vec, min_params_mat,
                        params_settings
                        ###### pass through prev term cost thing #####
                     ) {
  # initializations 
  params_cost_i <- param_cost_func(params_i, prev_params, params_settings)
  cuml_error <- 0
  p_0 <- params_i[[1]]
  k_1 <- params_i[[2]]
  k_2 <- params_i[[3]]
  coef_1 <- exp(-1 / params_i[[4]])
  coef_2 <- exp(-1 / params_i[[5]])

  T_1 <- 0
  T_2 <- 0
  lower_index <- 1
  j <- 0
  if (window > 0) { # 0 is the case where no window is specified
    error_tmp <- c(rep(0, window))
  }
  
  for (new_index in obs_indexes) {
    for (t in lower_index:new_index) {
      T_1 <- coef_1*(T_1 + k_1*training_load[[t]])
      T_2 <- coef_2*(T_2 + k_2*training_load[[t]])
    }
    error_i <- error_func(p_0 + T_1 - T_2 - obs_perf[[new_index]])
    cuml_error <- alpha*cuml_error + error_i 
    if(window > 0) {
      cuml_error <- cuml_error - alpha^(window)*error_tmp[[j%%window + 1]]
      error_tmp[[j%%window + 1]] <- error_i 
    }
    j <- j + 1
    lower_index <- new_index + 1
    cost_ij <- cost_error_func(cuml_error, j) + params_cost_i
      TE_cost_func(T_1, T_2, params_i[[4]], params_i[[5]],
                   prev_params[[4]], prev_params[[5]], params_settings)
    if (cost_ij < min_cost_vec[[j]]) {
      min_cost_vec[[j]] <- cost_ij
      min_params_mat[j, ] <- params_i
    }
  }
  return(list("params" = min_params_mat, "cost" = min_cost_vec))
}



## ----purl = TRUE----------------------------------------------------------------------------------------------------------------------------


update_cost_2 <- function(params_i, obs_indexes, training_load, obs_perf,
                          error_func, cost_error_func, param_cost_func, 
                          alpha, window, min_cost_vec, min_params_mat,
                          params_settings, tau_i_sum_mat, TE_mat
                 ) {
  # initializations 
  params_cost_i_vec <- param_cost_func(params_i, min_params_mat, params_settings, tau_i_sum_mat)
  cuml_error <- 0
  T_1 <- 0
  T_2 <- 0
  j <- 0
  if (window > 0) { # 0 is the case where no window is specified
    error_tmp <- c(rep(0, window))
  }
  
  for (new_index in obs_indexes) {
    error_i <- error_func(params_i[[1]] + params_i[[2]]*TE_mat[params_i[[6]], new_index] 
                          - params_i[[3]]*TE_mat[params_i[[7]], new_index] - obs_perf[[new_index]])
    cuml_error <- alpha*cuml_error + error_i 
    if(window > 0) {
      cuml_error <- cuml_error - alpha^(window)*error_tmp[[j%%window + 1]]
      error_tmp[[j%%window + 1]] <- error_i 
    }
    j <- j + 1
    cost_ij <- cost_error_func(cuml_error, j) + params_cost_i_vec[[j]]
    if (cost_ij < min_cost_vec[[j]]) {
      min_cost_vec[[j]] <- cost_ij
      min_params_mat[j, ] <- params_i # includes indexes for taus
    }
  }
  return(list("params" = min_params_mat, "cost" = min_cost_vec))
}
###################
## make sure tau_list works
###################


min_params <- function(params_matrix, init_params, # these two could change
                       training_load, obs_perf, # these won't; they are data
                       alpha, window, # these won't; they are specified by the user
                       error_func, cost_error_func, param_cost_func, # these won't; they are specified by the user
                       TE_cost_func, params_settings, grid_settings, 
                       num_update_cost) {
  obs_indexes <- which(!is.na(obs_perf))
  n <- length(obs_indexes)
  min_cost_vec <- c(rep(Inf, n))
  tau_list <- sort(unique(c(grid_settings$tau_1, grid_settings$tau_2)))
  
  M <- length(tau_list) + 2 # +2 for the initial parameter tau values, makes sense in the next chunk #
  init_params = c(init_params, M-1, M)
  min_params_mat <- matrix(init_params, nrow = n, ncol = 7, byrow = TRUE) 
  
  ##############################################################
  #### Precomputing some values for computational efficency ####
  ##############################################################
  
  N <- params_settings$N
  tau_i_vec <- c(as.numeric(purrr::map(tau_list, ~exp(-1/.x))),
                 exp(-1/init_params[[4]]),
                 exp(-1/init_params[[5]])
  )
  tau_i_sum_mat <- matrix(0, nrow = M, ncol = N)
  tau_i_sum_mat[, 1] <- tau_i_vec
  for (i in 2:N) {
    tau_i_sum_mat[, i] = tau_i_sum_mat[, i-1] + tau_i_vec^i
  }
  
  
  N_days <- length(training_load)
  #!!!!!!!!!! maybe a problem here !!!!!!!!!!!!
  TE_mat <- matrix(0, nrow = M, ncol = N_days)
  TE_mat[, 1] <- tau_i_vec*training_load[[1]]
  for (i in 2:N_days){
    TE_mat[, i] = tau_i_vec*(TE_mat[, i-1] + training_load[[i]])
  }
  
  #############################
  #### Doing the algorithm ####
  #############################
  
  for (i in 1:nrow(params_matrix)) {
    
    if (num_update_cost == 1) {
       res <- update_cost(as.numeric(params_matrix[i, ]), prev_params, obs_indexes, training_load, obs_perf,
                         error_func, cost_error_func, param_cost_func,
                         TE_cost_func,
                         alpha, window,
                         min_cost_vec, min_params_mat,
                         params_settings
                         ) 
    }
    
    if (num_update_cost == 2) {
      
      res <- update_cost_2(as.numeric(params_matrix[i, ]), obs_indexes,
                           training_load, obs_perf,
                           error_func, cost_error_func, param_cost_func,
                           alpha, window, min_cost_vec, min_params_mat,
                           params_settings,
                           tau_i_sum_mat, TE_mat) 
    }
    min_cost_vec <- res$cost
    min_params_mat <- res$params
  }
  return(list("opt_params" = min_params_mat, "cost" = min_cost_vec))
}

# this function creates (or calls, we will see) the parameter matrix and applies 
# it to the previous function. makes things tidier in the algorithm funciton
search_params_mat <- function(grid_settings, init_params, num_update_cost) {
  bounds_type <- grid_settings$type
  if(bounds_type == "test") {
    params_matrix <- expand.grid(p_0 = 500,
                                  k_1 = c(1,2,3,4),
                                  k_2 = c(2,4,6,8),
                                  tau_1 = c(5:35),
                                  tau_2 = c(5:35))
  }
  if(bounds_type == "custom" && num_update_cost == 1) {
    params_matrix <- expand.grid(p_0 = grid_settings$p_0,
                              k_1 = grid_settings$k_1,
                              k_2 = grid_settings$k_2,
                              tau_1 = grid_settings$tau_1,
                              tau_2 = grid_settings$tau_2)
  }
  
  if(bounds_type == "custom" && num_update_cost == 2) {
    params_matrix <- expand.grid(p_0 = grid_settings$p_0,
                              k_1 = grid_settings$k_1,
                              k_2 = grid_settings$k_2,
                              tau_1 = grid_settings$tau_1,
                              tau_2 = grid_settings$tau_2)
    
    tau_list <- sort(unique(c(grid_settings$tau_1, grid_settings$tau_2)))
    
    tau_1_index <- c(rep(0, nrow(params_matrix)))
    tau_2_index <- c(rep(0, nrow(params_matrix)))
    for (i in 1:nrow(params_matrix)) {
      tau_1_index[[i]] <- which(tau_list == params_matrix$tau_1[[i]])
      tau_2_index[[i]] <- which(tau_list == params_matrix$tau_2[[i]])
    }
    params_matrix <- cbind(params_matrix, tau_1_index)
    params_matrix <- cbind(params_matrix, tau_2_index)
    
    # special case for initial parameters
    params_matrix <- rbind(params_matrix, c(init_params, 
                                            length(tau_list) + 1,
                                            length(tau_list) + 2))
  }
  
  init_params <- init_params
  
  params_matrix <- params_matrix[params_matrix$tau_1 > params_matrix$tau_2, ]
  params_matrix <- params_matrix[params_matrix$k_1 < params_matrix$k_2, ]
  return(params_matrix)
}





## ----New algorithm, purl = TRUE-------------------------------------------------------------------------------------------------------------
new_pred_perf <- function(init_params,
                          training_load,
                          obs_perf,
                          error_settings = list("type" = "RMSE", "alpha" = 1, "window" = Inf),
                          params_settings,
                          grid_settings = list("type" = "test")) {
  
  curr_params <- c(init_params, 0,0)
  # curr_params take form c(p_0, k_1, k_2, tau_1, tau_2)
  ###
  matrix_params <- matrix(0, nrow = length(training_load), ncol = 7)
  colnames(matrix_params) <- c("p_0", "k_1", "k_2", "tau_1", "tau_2", "index_1", "index_2")
  days <- length(training_load)
  perf_out <- c(rep(curr_params[[1]], days + 1)) 
  
  ###############################################
  ########## Doing the error Settings ###########
  ###############################################
  
  type_error <- error_settings$type
  
  if (type_error %in% c("AE", "MAE", "AEOS")) {
    error_func <- abs
  } else if (type_error %in% c("SSE", "MSE", "RMSE")) {
    error_func <- function(x) x^2
  } else if (type_error %in% c("SPE", "MPE", "RMPE")) {
    p <- error_settings$p
    error_func <- function(x) x^p
  }
  else stop("change error settings")
  
  cost_error_func <- switch(type_error,
    "AE" = function(error, j) {error},
    "MAE" = function(error, j) {error/j},
    "AEOS" = function(error, j ) {error/sqrt(j)},
    "SSE" = function(error, j) {error},
    "MSE" = function(error, j) {error/j},
    "RMSE" = function(error, j) {sqrt(error/j)},
    "SPE" = function(error, j) {error},
    "MPE" = function(error, j) {error/j},
    "RMPE" = function(error, j) {(error/j)^{1/p}},
  )
  
  ###############################################
  ######### Doing the parameter Settings ########
  ###############################################
  
  if (params_settings$type == "alg_cost") {
    ## See helper function chunk for these functions ##
    
    param_cost_func <- alg_param_cost_func
    TE_cost_func <- alg_TE_cost_func
    num_update_cost <- 1
  }
  
  if (params_settings$type == "alg_cost_2") {
    param_cost_func <- alg_param_cost_func_2
    num_update_cost <- 2
  }
  
  #### checking inputs ####
  alpha <- error_settings$alpha
  if (is.null(alpha) == TRUE) {
    alpha = 1
  } else if (!is.numeric(alpha)){
    stop("bad value for alpha, must be a real number")
  } else if (alpha > 1 || alpha < 0) {
    warning("Alpha is intended to be a real number between 0 and 1") 
  } else {
    #pass 
  }

  window <- error_settings$window
  if(is.null(window)==TRUE || is.infinite(window)) {
    window = 0
  } else if (!is.numeric(window)|| window - floor(window) > 0 || window < 0) {
    stop("Bad Value of alpha; must be a positive integer")
  } else {
    #pass
  }
  
  
  ##############################################
  ########### Doing the Algorithm ##############
  ############################################## 
  
  res_1 <- min_params(search_params_mat(grid_settings, init_params, num_update_cost),
                      init_params, training_load, obs_perf,
                      alpha, window,  
                      error_func, cost_error_func, param_cost_func,
                      TE_cost_func, params_settings, grid_settings,
                      num_update_cost)
  opt_params_mat <- res_1$opt_params 
  cost_vec <- res_1$cost 

  ##############################################
  ######## Computing The performance ###########
  ##############################################
  
  T_1 <- 0
  T_2 <- 0
  j = 1
  for (i in 1:length(obs_perf)) {
    # update params if there is new information
    if (is.na(obs_perf[[i]]) == FALSE) {
      curr_params <- as.numeric(opt_params_mat[j, ])
      j <- j + 1
    } else {
    } # pass
    T_1 <- exp(-1 / curr_params[[4]]) * (T_1 + curr_params[[2]] * training_load[[i]])
    # training load index is  since i starts at 2
    T_2 <- exp(-1 / curr_params[[5]]) * (T_2 + curr_params[[3]] * training_load[[i]])
    perf_out[[i+1]] <- curr_params[[1]] +  T_1 - T_2
    matrix_params[i, ] <- as.numeric(curr_params)
  }
  
  ##############################################
  ########## plotting the results ##############
  ##############################################
  
  params_data <- tibble(
    "day" = 0:length(training_load),
    "k_1" = as.numeric(c(matrix_params[1, "k_1"], matrix_params[, "k_1"])),
    "k_2" = as.numeric(c(matrix_params[1, "k_2"], matrix_params[, "k_2"])),
    "tau_1" = as.numeric(c(matrix_params[1, "tau_1"], matrix_params[, "tau_1"])),
    "tau_2" = as.numeric(c(matrix_params[1, "tau_2"], matrix_params[, "tau_2"])),
    "p_0" = as.numeric(c(matrix_params[1, "p_0"], matrix_params[, "p_0"]))
  )
  k_plot <- ggplot(params_data, aes(x = day)) +
    geom_line(aes(y = k_1, color = "k_1")) +
    geom_line(aes(y = k_2, color = "k_2"))

  tau_plot <- ggplot(params_data, aes(x = day)) +
    geom_line(aes(y = tau_1, color = "tau_1")) +
    geom_line(aes(y = tau_2, color = "tau_2"))
  
  training_plot <- ggplot(tibble("day" = 0:length(training_load),
                                 "train_load" = c(0,training_load)
                                 ),
                          aes(x = day,
                              y= train_load)
                          ) +
    geom_bar(stat = "identity")

  #outputs
  out_list <- list()
  out_list$performance <- perf_out
  out_list$cost_vec <- cost_vec
  out_list$perf_plot <- plot_perf(obs_perf, perf_out)
  out_list$params_plots <- k_plot + tau_plot
  out_list$training_plot <- training_plot
  return(out_list)
} 

