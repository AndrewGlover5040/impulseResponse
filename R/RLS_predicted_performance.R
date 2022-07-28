#RLS
#
library(ggplot2)
# dat <- load("DFfromPDF (1).Rda")
# Athlete.data <- df
#
# day.1=Athlete.data[[1]]
# params.1=c(256, 0.10, 0.10, 15, 11)
# Training.1=Athlete.data[[2]]
# Performance.1=Athlete.data[[3]]
# #Training.1
#
#
# library(readr)
# Winther_data <- read_csv("Winther.data.csv")
#
#
# Athlete.data=read.csv("clarkeSkibaData - Sheet1.csv")
#
# day.1=Athlete.data[[1]]
# params.1=c(262, 0.18, 0.23, 36, 21)
# Training.1=Athlete.data[[2]]
# Performance.1=Athlete.data[[3]]
f

#assumes training load has a zero day
get_performance_matrix <- function(training_load,
                                   bounds,
                                   by){
  #need to check that this is a positive integer
  num_tau=(bounds[[2]]-bounds[[1]])/by+1
  e_to_the_tau_vec=c(rep(0,num_tau))
  for (i in 1:num_tau){
    e_to_the_tau_vec[[i]]=exp(-1/(bounds[[1]]+(i-1)*by))
  }
  e_to_the_tau_vec=as.matrix(e_to_the_tau_vec)
  max_day=length(training_load)
  out_matrix=matrix(0,num_tau,max_day)
for (i in 2:max_day){
    out_matrix[,i]=e_to_the_tau_vec*(out_matrix[,i-1])+training_load[[i]]
  }

  out_matrix
}


##assume T_2 is posititve
RLS_Algorithm <- function(training_load,
                          performance,
                          indexes_perf,
                          initializations,
                          p_0,
                          T_1_vec,
                          T_2_vec,
                          alpha=1
){
  #initializations
  ##fix theta to be the right thing
  theta = initializations[[1]]
  #=matrix(c(.1,0),byrow = TRUE)

  P = initializations[[2]]
  #=delta*matrix(c(1,0,0,1),nrow=2)

  K=initializations[[3]]
  #=matrix(0,2,1)

  #make dim two identity matrix
  I=initializations[[4]]
  #=matrix(c(1,0,0,1),nrow=2)

  S_vector=initializations[[5]]
  #=numeric(len_indexes_perf)

  #sum
  S=initializations[[6]]
  #=0

  #counter
  j=initializations[[7]]
  #=1

  theta_list=initializations[[8]]
  #=as.list(rep(NA,len_indexes_perf))

  #actual recursion algorithm
  for (i in indexes_perf){
    x=matrix(c(T_1_vec[[i]],T_2_vec[[i]]),byrow=TRUE)

    #so we don't have to compute this twice
    denom=as.numeric(alpha+t(x)%*%P%*%x)

    #update K (with old P)
    K=P%*%x/denom

    #update P ##fix to make more efficent
    P=(I-P%*%x%*%t(x)/denom)%*%P/alpha

    #update theta
    error = performance[[i]]-as.numeric(p_0+t(x)%*%theta)
    theta = theta+K*(error)
    S = S*alpha+error^2
    S_vector[[j]] = S
    theta_list[[j]] = theta
    j=j+1
  }
  list(theta_list,S_vector)
}

#For testing RLS function
# RLS_Algorithm(Training.1, Performance.1,
#     indexes_perf, initializations,
#     256, T_matrix[36,], -T_matrix[36,],
#     alpha=.98)




#maybe make initializations list in the server environment
estimate_RLS_parameters <- function(training_load,
                                    performance,
                                    p_0,
                                    alpha=1,
                                    delta=1,
                                    bounds_T_1,
                                    by_T_1,
                                    bounds_T_2,
                                    by_T_2
){
  T_1_matrix = get_performance_matrix(training_load,bounds_T_1,by_T_1)
  T_2_matrix = get_performance_matrix(training_load,bounds_T_2,by_T_2)
  num_tau_1 = (bounds_T_1[[2]]-bounds_T_1[[1]])/by_T_1+1
  num_tau_2 = (bounds_T_2[[2]]-bounds_T_2[[1]])/by_T_2+1
  indexes_perf = which(!is.na(performance))
  len_indexes_perf = length(indexes_perf)

  ##fix theta
  RLS_initializations = list(matrix(c(.1,0),byrow = TRUE),
                             delta*matrix(c(1,0,0,1),nrow=2),
                             matrix(0,2,1),
                             matrix(c(1,0,0,1),nrow=2),
                             numeric(len_indexes_perf),
                             0,
                             1,
                             as.list(rep(NA,len_indexes_perf))
  )

  cost_array <- array(0,c(len_indexes_perf,num_tau_1,num_tau_2))
  theta_list <- list(rep(NA,num_tau_2*num_tau_1))

  #do the RLS algorithm for all combinations of tau_1 and tau_2
  for (i in 1:num_tau_1){
    for (j in 1:num_tau_2){
      tmp_ans = RLS_Algorithm(training_load,
                              performance,
                              indexes_perf,
                              RLS_initializations,
                              p_0,
                              T_1_matrix[i,],
                              -T_2_matrix[j,],
                              alpha=alpha
      )
      theta_list[[i+num_tau_1*(j-1)]] = tmp_ans[[1]]
      cost_array[,i,j] = tmp_ans[[2]]
    }
  }


  min_cost_vec <- numeric(len_indexes_perf)
  #min_cost_index_list <- list(rep(NA,len_indexes_perf))
  best_parameter_list <- as.list(rep(NA,len_indexes_perf))
  sm=0
  for (i in 1:len_indexes_perf){
    tmp = which.min(cost_array[i,,])
    j = (tmp-1)%%num_tau_1+1
    k = (tmp-1)%/%num_tau_1+1
    #min_cost_index_list[[i]]=c(j,k)
    if (i==1){
      min_cost_vec[[i]]=cost_array[i,j,k]
      sm=cost_array[i,j,k]
    }
    else{
      min_cost_vec[[i]]=cost_array[i,j,k]-alpha*cost_array[i-1,j,k]
      sm=sm+cost_array[i,j,k]-alpha*cost_array[i-1,j,k]
    }
    best_parameter_list[[i]] = as.list(theta_list[[j+num_tau_1*(k-1)]][[i]])
    best_parameter_list[[i]] = unlist(c(list(p_0), best_parameter_list[[i]],
                                        list(bounds_T_1[[1]]+(j-1)*by_T_1,
                                             bounds_T_2[[1]]+(k-1)*by_T_2
                                        )
    )
    )
  }

  list(params = best_parameter_list, SSE = min_cost_vec)
}

#
# tmp = estimate_RLS_parameters(Training.1,
#                       Performance.1,
#                       256,
#                       alpha = 1,
#                       delta = 1,
#                       bounds_T_1 = c(1,51),
#                       by_T_1 = 1,
#                       bounds_T_2 = c(1,51),
#                       by_T_2 = 1
#                       )
# tmp


numeric_predictedPerformance <- function(params,day=length(Training.Load),Training.Load){
  p_0=params[[1]]; k_1=params[[2]]; k_2=params[[3]]; tau_1=params[[4]]; tau_2=params[[5]]
  T_1=0; T_2=0
  coef_1=exp(-1/tau_1); coef_2=exp(-1/tau_2)
  for(t in 1:day){
    T_1=coef_1*T_1+Training.Load[[t]]
    T_2=coef_2*T_2+Training.Load[[t]]
  }
  return(p_0+k_1*T_1-k_2*T_2)
}



RLS_predicted_performance <- function(training_load,
                                      performance,
                                      p_0,
                                      alpha=1,
                                      delta=1,
                                      bounds_T_1,
                                      by_T_1,
                                      bounds_T_2,
                                      by_T_2){
  #get the parameters
  tmp <- estimate_RLS_parameters(training_load,
                                 performance,
                                 p_0,
                                 alpha,
                                 delta,
                                 bounds_T_1,
                                 by_T_1,
                                 bounds_T_2,
                                 by_T_2)
  params_RLS = tmp[[1]]
  SSE_list = tmp[[2]]
  params_list = as.list(rep(NA,length(performance)))
  j=1
  for(i in which(!is.na(performance))){
    params_list[[i]] = params_RLS[[j]]
    j=j+1
  }
  ####!!! make first parameter dependent on context !!!####
  params_list[[1]]=c(p_0,.1,0,50,0)

  for (i in 2:length(params_list)){
    if(is.na(params_list[[i]][[1]])==TRUE){
      params_list[[i]]=params_list[[i-1]]
    }
  }

  day=as.list(1:length(training_load))
  out_list=purrr::map2(params_list,day,numeric_predictedPerformance,training_load)

  return(list(out_list,params_RLS,SSE_list))
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
