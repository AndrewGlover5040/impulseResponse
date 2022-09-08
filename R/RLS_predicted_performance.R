library(ggplot2)

####!!!! make sure that we get an integer number of elements of R!!!!####
#assumes training load has a zero day
#' Get Performance Matrix
#'
#' Returns a Matrix of the value of T(n), computed for numerical vector of tau.
#' Lets refer to this vector as R. This is a helper function for
#' estimate_RLS_parameters (and ultimately RLS_predicted_performance).
#'
#' @param training_load A numerical vector
#' @param bounds A length-two numerical vector. The first entry is the left
#'   endpoint R, and the second entry is the right endpoint R.
#' @param by A number. This is the increment between values of R.
#'
#' @return A Matrix with dimensions length(R) by length(training_load)
#'
#'
#' @examples
#' # For some training load w, we can get the time-invariant prediction of
#' performance by calling
#' [p_0 + k_1*get_performance_matrix(w,c(tau_1,tau_1),1) - k_2*get_performance_matrix(w,c(tau_2,tau_2),1)]
#' This agrees with predictedPerformance(c(p_0, k_1, k_2, tau_1, tau_2), w)
#'
get_performance_matrix <- function(
    training_load,
    bounds,
    by
){
  #need to check that this is a positive integer
  num_tau=(bounds[[2]] - bounds[[1]])/by + 1
  e_to_the_tau_vec=c(rep(0,num_tau))
  for (i in 1:num_tau){
    e_to_the_tau_vec[[i]] = exp(-1/(bounds[[1]] + (i-1)*by))
   }
  e_to_the_tau_vec = as.matrix(e_to_the_tau_vec)
  max_day = length(training_load)
  out_matrix = matrix(0,num_tau,max_day)
  for (i in 2:max_day){
    out_matrix[,i] = e_to_the_tau_vec*(out_matrix[,i-1]) + training_load[[i]]
  }

  out_matrix
}


####assume T_2 is posititvea
#' RLS Algorithm
#'
#' The implementation of the RLS algorithm with a forgetting factor.
#' This is used as a helper function for `estimate_RLS_parameters()` and
#' ultimately for `RLS_predicted_performance()`
#'
#' @inheritParams RLS_predicted_performance
#' @param indexes_perf A list. Passed from the function that calls it
#' @param fixed_initializations A list. This passes on initalizations of the
#'   Algorithm from the function that calls it. This saves computation time
#'   recreating all of these values on every iteration.
#' @param p_0 A number. The performance on day 0.
#' @param T_1_vec A 1 by length(training_load) matrix. A row of
#'   get_performance_matrix(.).
#' @param T_2_vec Same as above
#'
#' @return A list of 2. The first is a list of values of the estimates of theta.
#'  The second is a matrix of the SSE of the estimates
#'
#'
#' @examples ADD in
RLS_Algorithm <- function(
    training_load,
    performance,
    indexes_perf,
    fixed_initializations,
    p_0,
    T_1_vec,
    T_2_vec,
    alpha=1
){
  #initializations
  ##!!!fix theta to be the right thing
  theta = fixed_initializations[[1]]
  #=matrix(c(.1,0),byrow = TRUE)

  P = fixed_initializations[[2]]
  #=delta*matrix(c(1,0,0,1),nrow=2)

  K = fixed_initializations[[3]]
  #=matrix(0,2,1)

  #make dim two identity matrix
  I = fixed_initializations[[4]]
  #=matrix(c(1,0,0,1),nrow=2)

  S_vector = fixed_initializations[[5]]
  #=numeric(len_indexes_perf)

  #sum
  S = fixed_initializations[[6]]
  #=0

  #counter
  j = fixed_initializations[[7]]
  #=1

  theta_list = fixed_initializations[[8]]
  #=as.list(rep(NA,len_indexes_perf))

  #actual recursion algorithm
  for (i in indexes_perf) {
    x = matrix(c(T_1_vec[[i]], T_2_vec[[i]]), byrow=TRUE)

    #so we don't have to compute this twice
    denom = as.numeric(alpha + t(x)%*%P%*%x)

    #update K (with old P)
    K = P%*%x/denom

    #update P ##fix to make more efficent
    P = (I - P%*%x%*%t(x)/denom)%*%P/alpha

    #update theta
    error = performance[[i]] - as.numeric(p_0 + t(x)%*%theta)
    theta = theta + K*(error)
    S = S*alpha+error^2
    S_vector[[j]] = S
    theta_list[[j]] = theta
    j = j + 1
  }
  list(theta_list, S_vector)
}



#' Estimate RLS Parameters
#'
#' To calculate an obtain the parameter estimates for the Time-invariant Model.
#'
#' @inheritParams RLS_predicted_performance
#'
#' @return A list of 2. The first is a list of length(indexes_perf) that gives
#'  the parameter estimates on the corresponding days. The second is a numerical
#'  vector of the marginal SSE
#'
estimate_RLS_parameters <- function(
    training_load,
    performance,
    p_0,
    alpha=1,
    delta=1,
    bounds_T_1,
    by_T_1,
    bounds_T_2,
    by_T_2
){
  #setting up the RLS algorithm.
  T_1_matrix <- get_performance_matrix(training_load, bounds_T_1, by_T_1)
  T_2_matrix <- get_performance_matrix(training_load, bounds_T_2, by_T_2)
  num_tau_1 <- (bounds_T_1[[2]] - bounds_T_1[[1]])/by_T_1 + 1
  num_tau_2 <- (bounds_T_2[[2]] - bounds_T_2[[1]])/by_T_2 + 1
  indexes_perf <- which(!is.na(performance))
  len_indexes_perf <- length(indexes_perf)

  ##!!!fix theta
  # Since a lot of the initializations are the same for each RLS,
  # this saves some time by just assigning the same objects to each
  # iteration of the algorithm
  RLS_initializations = list(
    matrix(c(.1,0), byrow = TRUE),
    delta*matrix(c(1,0,0,1), nrow=2),
    matrix(0,2,1),
    matrix(c(1,0,0,1), nrow=2),
    numeric(len_indexes_perf),
    0,
    1,
    as.list(rep(NA, len_indexes_perf))
  )

  #to store the results
  cost_array <- array(0, c(len_indexes_perf, num_tau_1, num_tau_2))
  theta_list <- list(rep(NA, num_tau_2*num_tau_1))

  # do the RLS algorithm for all combinations of tau_1 and tau_2
  for (i in 1:num_tau_1) {
    for (j in 1:num_tau_2) {
      tmp_ans = RLS_Algorithm(
        training_load,
        performance,
        indexes_perf,
        RLS_initializations,
        p_0,
        T_1_matrix[i,],
        -T_2_matrix[j,],
        alpha=alpha
      )
      theta_list[[i + num_tau_1*(j-1)]] = tmp_ans[[1]]
      cost_array[,i,j] = tmp_ans[[2]]
    }
  }

  # refines the best estimates to a usable output
  min_cost_vec <- numeric(len_indexes_perf)
  best_parameter_list <- as.list(rep(NA,len_indexes_perf))
  sm=0
  for (i in 1:len_indexes_perf) {
    tmp = which.min(cost_array[i,,])
    j = (tmp-1)%%num_tau_1 + 1
    k = (tmp-1)%/%num_tau_1 + 1
    if (i==1){
      min_cost_vec[[i]] = cost_array[i, j, k]
      sm = cost_array[i, j, k]
    }
    else{
      min_cost_vec[[i]] = cost_array[i, j, k]-alpha*cost_array[i-1, j, k]
      sm = sm + cost_array[i,j,k] - alpha*cost_array[i-1, j, k]
    }
    best_parameter_list[[i]] = as.list(theta_list[[j+num_tau_1*(k-1)]][[i]])
    best_parameter_list[[i]] = unlist(
      c(list(p_0),
        best_parameter_list[[i]],
        list(bounds_T_1[[1]] + (j-1)*by_T_1,
             bounds_T_2[[1]] + (k-1)*by_T_2
        )
      )
    )
  }
  list(best_parameter_list, min_cost_vec)
}


#' Helper function for RLS
#'
#' @inheritParams Influence
#' @param day A number. Tells the function what day to compute the performance
#'   on.
#' @inheritParams RLS_predicted_performance
#'
#' @return A number
#'
#' @examples
numeric_predictedPerformance <- function(
    params,
    day = length(training_load),
    training_load
){
  p_0 = params[[1]]
  k_1 = params[[2]]; k_2 = params[[3]]
  tau_1 = params[[4]]; tau_2 = params[[5]]
  T_1=0; T_2=0
  coef_1 = exp(-1/tau_1); coef_2 = exp(-1/tau_2)
  for (t in 1:day) {
    T_1 = coef_1*T_1 + training_load[[t]]
    T_2 = coef_2*T_2 + training_load[[t]]
  }
  p_0 + k_1*T_1 - k_2*T_2
}



#' The augmented RLS predition of performance.
#'
#' @param training_load A numerical vector. In the i'th spot contains the
#'   Training impulse on day i.
#' @param performance
#' @param p_0 The performance on day 0.
#' @param alpha A positve number in (0,1]. This is the "forgetting factor" in the
#'   RLS algorithm.
#' @param delta A non-negative number. Used to initialize a positive
#'   definite matrix in the RLS algorithm
#' @param bounds_T_1 Specifies the bounds in R for T_1
#'   (See `get_performance_matrix`)
#' @param by_T_1 Specifies the increment for R in T_1
#'   (See `get_performance_matrix`)
#' @param bounds_T_2 Specifies the bounds in R for T_2
#'   (See `get_performance_matrix`)
#' @param by_T_2 Specifies the increment for R in T_1
#'   (See `get_performance_matrix`)
#'
#' @return A list containing the predicted performance from the RLS algorithm,
#'   The list containing the parameters used to calculate the performance, and
#'   and the marginal SSE.
#' @export
#'
#' @examples fill in later
RLS_predicted_performance <- function(
    training_load,
    performance,
    p_0,
    alpha = 1,
    delta = 1,
    bounds_T_1,
    by_T_1,
    bounds_T_2,
    by_T_2
){
  #call the helper function
  tmp <- estimate_RLS_parameters(
    training_load,
    performance,
    p_0,
    alpha,
    delta,
    bounds_T_1,
    by_T_1,
    bounds_T_2,
    by_T_2
  )
  params_RLS <- tmp[[1]]
  SSE_list <- tmp[[2]]

  # To make the output list of parameters for each of the
  # days when we have parameters
  params_list = as.list(rep(NA,length(performance)))
  j = 1
  for(i in which(!is.na(performance))) {
    params_list[[i]] = params_RLS[[j]]
    j = j + 1
  }
  ####!!! make first parameter dependent on context !!!####
  params_list[[1]] = c(p_0,.1,0,50,0)

  #Fill in NA where there is no performance data.
  for (i in 2:length(params_list)) {
    if (is.na(params_list[[i]][[1]])==TRUE) {
      params_list[[i]] = params_list[[i-1]]
    }
  }

  # Going from parameters to predicted performance.
  day=as.list(1:length(training_load))
  out_list = purrr::map2(
    params_list,
    day,
    numeric_predictedPerformance,
    training_load
  )

  list(out_list, params_RLS, SSE_list)
}


# predicted_performance <- RLS_predicted_performance(Training.1,
#                                             Performance.1,
#                                  262,
#                                  alpha = .95,
#                                  delta = 1000,
#                                  bounds_T_1 = c(1,50),
#                                  by_T_1 = 1,
#                                  bounds_T_2 = c(1,50),
#                                  by_T_2 = 1
#                                  )

# predicted=unlist(predicted_performance)
#
# df.1=data.frame(day.1,predicted,Performance.1)
#
#
# plot=ggplot(df.1,aes(x=day.1, sy=predicted))+
#   geom_line(aes(y=predicted, colour="black"), size=1)+
#   geom_point(aes(y = Performance.1, color="red"), shape = 1)+
#   scale_color_manual("", values = c("black", "red"),
#                      labels = c("Predicted Performance",
#                                 "Actual Performance")
#   )
#
# plot
