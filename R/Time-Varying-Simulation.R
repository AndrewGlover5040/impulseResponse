## ----libraries, purl = TRUE------------------------------------------------------------------------------------------------------------------------------
library(ggplot2) # for graphing
library(patchwork) # to add graphs together
library(tibble) # tibbles




## ----Generating the parameters matrix, purl=TRUE---------------------------------------------------------------------------------------------------------
params_mat <- function(k_1, tau_1, k_2, tau_2, change_days=NULL, days) {
  if (length(k_1) == length(tau_1) && 
      length(k_1) == length(tau_2) &&
      length(k_1) == length(k_2) &&
      length(k_1) == length(change_days)+1
  ) {}
    else {stop("check length of parameters")} 
  out_matrix <- matrix(0, nrow = days, ncol = 4)
  colnames(out_matrix) <- c("k_1", "tau_1", "k_2", "tau_2")
  bound_1 <- 1; bound_2 <- days
  j <- 0 # counter for index of k_1, tau_1, etc
  for (elem in c(change_days, days)) {
    j <- j + 1
    bound_2 <- elem 
    for (i in bound_1:bound_2) {
      out_matrix[i, ] <- c(k_1[[j]], tau_1[[j]], k_2[[j]], tau_2[[j]])
    }
    bound_1 <- elem
  }
  return(out_matrix)
}


## ----Varying Performance, purl = TRUE--------------------------------------------------------------------------------------------------------------------
#' Title
#'
#' @param params_mat a n by 4 matrix, where n is the number of days
#' @param training_load an n dimensional vector 
#'
#' @return an n? dimensional vector with the performance
#' @export
#'
#' @examples
mat_to_perf <- function(p_0, params_mat, training_load) {
  days <- nrow(params_mat)
  perf_out <- c(rep(NA, days))
  T_1 <- 0; T_2 <- 0
  for (i in 1:days) {
    T_1 <- exp(-1/params_mat[i, "tau_1"])*T_1 + training_load[[i]]
    T_2 <- exp(-1/params_mat[i, "tau_2"])*T_2 + training_load[[i]]
    perf_out[[i]] <- p_0 + params_mat[i, "k_1"]*T_1 - params_mat[i, "k_2"]*T_2 
  }
  return(perf_out)
}


## ----purl = TRUE-----------------------------------------------------------------------------------------------------------------------------------------
perf_tv <- function(p_0,
                      k_1,
                      tau_1,
                      k_2,
                      tau_2,
                      change_days = NULL,
                      days,
                      training_stim,
                      lim = TRUE) {
  training_load <- c()
  limit <- 0
  if (training_stim[[1]] == "constant") {
    training_load <- c(rep(training_stim[[2]], days))
    k <- length(k_1)
    
    # See justification for this below
    limit <- p_0 + training_stim[[2]]*k_1[[k]]/(1-exp(-1/tau_1[[k]])) -
      training_stim[[2]]*k_2[[k]]/(1-exp(-1/tau_2[[k]])) 
  }
  tmp_matrix <- params_mat(k_1,
                           tau_1,
                           k_2,
                           tau_2,
                           change_days,
                           days)
  modeled_performance <- mat_to_perf(p_0, tmp_matrix, training_load)
  tmp_data <- tibble(
    "day" = c(0:days),
    "performance" = c(p_0, modeled_performance),
    "limit" = c(rep(limit, days + 1))
  )
  
  plot <- ggplot(tmp_data, aes(x = day)) +
    geom_point(aes(y = performance, color = "perf")) +
  if (lim == TRUE){
     plot <- plot + geom_line(aes(y = limit, color = "lim"))
  }
 
  
  #   scale_color_manual("Legend",
  #                      values = c("lim" = "#e31a1c", # this color comes from the theme "Paired"
  #                                 "perf" = "black"))
  out_list <- list()
  out_list$plot <- plot
  out_list$performance <- modeled_performance
  return(out_list)
}


