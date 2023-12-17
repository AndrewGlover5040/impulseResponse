##############Clarke and Skiba Figure 6C
p0 <- 500
k1 <- 1
k2 <- 2
tau1 <- 27
tau2 <- 10
n <- 200
w <- c(rep(100,120),rep(30,7),rep(0,73))

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

plot(2:n,NF_vec,ylim=c(0,3000),type="l",lty=2,col="red")
lines(PF_vec,col="blue",lty=2)
lines(p0+PF_vec-NF_vec,type="l",col="green")
