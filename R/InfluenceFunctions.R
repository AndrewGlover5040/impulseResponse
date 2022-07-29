#for testing the influence function in clark&skiba fig. 8b.

library(devtools)
####!!! check this is outputing the right value !!!####
#' Influcence Function
#'
#' @param params A numerical vector of length five. It encodes
#'   [c(p_0, k_1, k_2, tau_1, tau_2)].
#' @param startDay The day that training begins at
#' @param t_p The day that the performance is supposed to happen
#'
#' @return A vector of length t_p-startday
#' @export
#'
#' @examples lk;j
Influence=function(params, startDay, t_p = 0){
  k_1=params[2]; k_2=params[3]; tau_1=params[4]; tau_2=params[5]
  n=t_p-startDay-1
  out=c(rep(0,n))
  for(i in 1:n){
    out[i]=k_1*exp(-(n-i+1)/tau_1)-k_2*exp(-(n-i+1)/tau_2)
  }
  return(out)
}

#' Title
#'
#' @param params A numerical vector of length five. It encodes
#'   [c(p_0, k_1, k_2, tau_1, tau_2)].
#'
#' @return A number
#' @export
#'
#' @examples
get_t_n=function(params){
  k_1=params[2]; k_2=params[3]; tau_1=params[4]; tau_2=params[5]
  t_g=(tau_1*tau_2)/(tau_1-tau_2)*log(k_2/k_1)
}

#' get_t_g
#'
#' @param params A numerical vector of length five. It encodes
#'   [c(p_0, k_1, k_2, tau_1, tau_2)].
#'
#' @return
#' @export
#'
#' @examples
get_t_g=function(params){
  k_1=params[2]; k_2=params[3]; tau_1=params[4]; tau_2=params[5]
  t_g=(tau_1*tau_2)/(tau_1-tau_2)*log((k_2*tau_1)/(k_1*tau_2))
  return(t_g)
}





