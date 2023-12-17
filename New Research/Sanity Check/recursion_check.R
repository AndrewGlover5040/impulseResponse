
days_test <- 500
perf <- function(params,
                 training_load,
                 t) {
  p_0=params[[1]]; k_1=params[[2]]; k_2=params[[3]]; tau_1=params[[4]]; tau_2=params[[5]]
  T_1= 0; T_2 = 0
  for (i in 1:(t-1)){
    T_1 = T_1 + exp(-(t-i)/tau_1)*training_load[[i+1]]
    T_2 = T_2 + exp(-(t-i)/tau_2)*training_load[[i+1]]
  }
  return(p_0+k_1*T_1-k_2*T_2)
}

params <- c(500, 1, 2, 27, 10)
training_load <- c(rep(100, days_test))
perf(params,
     training_load,
     500)

perf_sim_1 <- c(500, perf_tv(p_0 = 500,
                             k_1 = 1,
                             tau_1 = 25,
                             k_2 = 2,
                             tau_2 = 10,
                             days = days_test,
                             training_stim = list("constant", 100))$performance
)
perf_sim_1

invariant_perf <- function(params, training_load, day=length(training_load)){
  p_0=params[[1]]; k_1=params[[2]]; k_2=params[[3]]; tau_1=params[[4]]; tau_2=params[[5]]
  out=c(rep(0,day))
  T_1=0; T_2=0
  coef_1=exp(-1/tau_1)
  coef_2=exp(-1/tau_2)
  for(t in 1:day){
    T_1=coef_1*(T_1+training_load[[t]])
    T_2=coef_2*(T_2+training_load[[t]])
    out[[t]]=p_0+k_1*T_1-k_2*T_2
  }
  return(out)
}

invariant_perf(params,
               training_load)


##############Clarke and Skiba Figure 6C
p0 <- 500
k1 <- 1
k2 <- 2
tau1 <- 27
tau2 <- 10
n <- 200
w <- c(rep(100,200))

##################################### Define Functions
NFn <- function(k2,w,day,tau2){
  temp <- rep(0,(day-1))
  for(i in 1:(day-1)){
    temp[i] <- exp(-(day-i)/tau2)*w[i]
  }
  return(k2*sum(temp))
}

PFn <- function(n,k1,tau1,w){
  temp <- rep(0,(day-1))
  for(i in 1:(day-1)){
    temp[i] <- exp(-(day-i)/tau1)*w[i]
  }
  return(k1*sum(temp))
}




########################################
#Create Negative and Positive training effects by calling functions
#Negative
NF_vec <- rep(0,(n-1))
for (i in 2:n){
  NF_vec[(i-1)] <- NFn(k2,w,i,tau2)
}

#Positive
PF_vec <- rep(0,(n-1))
for (i in 2:n){
  PF_vec[(i-1)] <- NFn(k1,w,i,tau1)
}

p0+PF_vec-NF_vec


plot(2:n,NF_vec,ylim=c(0,3000),type="l",lty=2,col="red")
lines(PF_vec,col="blue",lty=2)
lines(p0+PF_vec-NF_vec,type="l",col="green")
