library(readxl)
Research <- read_excel("C:/Users/amglo/Downloads/Research.xlsx")
View(Research)
library(tidyr)


data_1<- read_excel("C:/Users/amglo/Downloads/better_table.xlsx",sheet ="Sheet1")
training_load=data_1$"training impulse"

performance=data_2$actual_perf
data_2=read_excel("C:/Users/amglo/Downloads/better_table.xlsx",sheet ="Sheet2")


View(better_table)

# names(tableS1)[names(tableS1)=="FIGURE 8: FITTING THE IR MODEL FOR A CYCLIST & DERIVATION OF THE CORRESPONDING INFLUENCE CURVE"]="...1"
# drop_na(tableS1,"...1")
# Days=tableS1$...4

v=c(262,0.18,0.23,36,21)
v

w <- c(rep(100,120),rep(30,7),rep(0,73))

TE=function(day,k,tau,w,vec=T){
  sum=0
  if (vec==T){
    values=c()
    for(t in 2:(day)){
      for (s in 1:(t-1)){
        sum=sum+exp(-(t-s)/tau)*w[s]
      }
      values=append(values,k*sum)
      sum=0
    }
    return(values)
  }
  else{
    for (s in 1:(day-1)){
      sum=sum+exp(-(day-s)/tau)*w[s]
    }
    return(k*sum)
  }
  # if(vec==T){
  #   return(values)
  # }else{
  #   return(k*sum)
  # }
}


P=function(day,p_0,k1,k2,w_1,w_2,dat,vec=T){
  output=c(rep(p_0,day-1))
  out=p_0
  if(vec==T){
    output=output+TE(day,k1,w_1,dat)-TE(day,k2,w_2,dat)
    return(output)
  }
  else{
    out=out+TE(day,k1,w_1,dat,F)-TE(day,k2,w_2,dat,F)
  }
}


out <- p_0+TE(200,1,27,w)-TE(200,2,10,w)
out
# P=P(200,500,1,2,27,10,w)
#
# PTE
# NTE
# P
#
# plot(2:200,NTE,ylim=c(0,3000),type="l",lty=2,col="red")
# lines(PTE,col="blue",lty=2)
# lines(P,type="l",col="green")

SSE=function(v,Training.Load,Performance,day=166){
  ##v=(p_0,k1,k2,w_1,w_2)
  p_0=v[1]; k1=v[2]; k2=v[3]; w_1=v[4]; w_2=v[5]

  Predicted.Performance=P(day,p_0,k1,k2,w_1,w_2,Training.Load,vec=T)
  Performance=Performance[1:day-1]
  errors=c(rep(0,day-1))
  for(i in 1:(day-1)){
    if(Performance[i]!=0){
    errors[i]= Performance[i] - Predicted.Performance[i]
    }
  }
  SSE=0
  for (i in 1:length(errors)){
    SSE <- SSE+errors[i]^2
  }
  return(SSE)
}

val=SSE(v,training_load,performance)
val
rlang::last_error()

optim_par=function(v,Training.Load,Performance,day=length(Performance)){
  x=optim(par = c(256, 0.10, 0.10, 15, 11),
          fn = SSE, Training.Load = Training.Load,
          Performance = Performance,
          day=day)
  return(x)
}

