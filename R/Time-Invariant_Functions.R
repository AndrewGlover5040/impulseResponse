#Time invariant functions

#model functions
#clean up parameters to function

#' predictedPerformance
#'
#' Banister IR Time Invariant Predicted Performance
#'
#' @param params A length-five numerical vector giving the value of (in order)
#'   p_0, k_1, k_2, tau_1, tau_2.
#' @param training_load A numerical vector containing the training load.
#'   Can be any length.
#' @param day A number that specifies the length of the return vector.
#'   Defaults to the length of `training_load`.
#' @returns A numerical vector the length of `day`.
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


SSE <- function(params,training_load,Performance,day=length(training_load)){
  Pred=predictedPerformance(params,training_load,day)
  #Performance=Performance[1:day-1]
  error=Performance[!is.na(Performance)] - Pred[which(!is.na(Performance))]
  error=error[!is.na(error)]
  SSE=sum(error^2)
  return(SSE)
}

optim_par <- function(v,training_load,Performance,day=length(Performance)){
  x=optim(par = v,
          fn = SSE, training_load = training_load,
          Performance = Performance,
          day=day)
  return(x$par)
}

