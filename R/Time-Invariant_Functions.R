#Time invariant functions

#model functions
#clean up parameters to function

#' predictedPerformance
#'
#' Banister IR Time Invariant Predicted Performance
#'
#' @inheritParams Influence
#' @inheritParams RLS_predicted_performance
#' @param day A number that specifies the length of the return vector.
#'   Defaults to the length of `training_load`.
#'
#' @returns A numerical vector the length of `day`.
#' @export
#'
predictedPerformance=function(params,training_load,day=length(training_load)){
  p_0=params[[1]]; k_1=params[[2]]; k_2=params[[3]]; tau_1=params[[4]]; tau_2=params[[5]]
  out=c(rep(0,day))
  T_1=0; T_2=0
  coef_1=exp(-1/tau_1)
  coef_2=exp(-1/tau_2)
  for(t in 1:day){
    T_1=coef_1*(T_1)+training_load[[t]]
    T_2=coef_2*(T_2)+training_load[[t]]
    out[[t]]=p_0+k_1*T_1-k_2*T_2
  }
  return(out)
}


#' SSE
#'
#' The SSE of a set of a model.
#'
#' @inheritParams predictedPerformance
#' @param performance A Numerical vector that is the same
#'   length as `training_load`. This contains the data of the
#'   actual performance of the athlete.
#'
#' @return A number
#' @export
#'
#' @examples
SSE <- function(
    params,
    training_load,
    performance,
    day = length(training_load)
){
  Pred <- predictedPerformance(params,training_load,day)
  #Performance=Performance[1:day-1]
  error <- performance[!is.na(performance)] - Pred[which(!is.na(performance))]
  error <- error[!is.na(error)]
  SSE <- sum(error^2)
  SSE
}

#' Optimize Parameters
#'
#' Return an optimization of the time-invariant parameters, using R's built-in
#' optimization function.
#'
#' @inheritParams predictedPerformance
#'
#' @return A number
#'
optim_par <- function(
    params,
    training_load,
    Performance,
    day=length(Performance)
){
  x <- optim(
    par = v,
    fn = SSE,
    training_load = training_load,
    Performance = Performance,
    day=day
  )
  x$par
}
