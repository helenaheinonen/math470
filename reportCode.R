library(extRemes)
library(ismev)
library(copula)

load("MaxAnnuelsAll.RData")
load("MaxMensuels.RData")

# fit GEV to yearly block maxima
aixYearly.data <- data.frame(year = 1976:2015, rainfall=maxann.data[18,])
out_aix_yearly <- gev.fit(aixYearly.data$rainfall) 

# fit GEV to monthly block maxima
aixMonthly.data<-as.vector(t(maxmonth.data[18,,]))
out_aix_monthly<-gev.fit(aixMonthly.data)

#diagnostic plots for each model
#observe deviation for monthly model
#par(mfrow=c(1,2))
gev.diag(out_aix_yearly)
gev.diag(out_aix_monthly)

#profile likelihood for 5 and 20 year return levels
par(mfrow=c(1,2))
gev.prof(out_aix_yearly,m=5,xlow=50,xup=150,nint=1000)
title("5 year return level based on yearly data")
gev.prof(out_aix_yearly,m=20,xlow=65,xup=200,nint=1000)
title("20 year return level based on yearly data")


# largest 2 order statistics copula

load("MonthlyOrderStat.RData")
aixMonthlyOrderStat<-as.array(monthOrderStat.data[,18,,])

aixMonthlyOrderStat.sample <- cbind(
  Max = as.vector(aixMonthlyOrderStat[1,,]),
  SecondMax = as.vector(aixMonthlyOrderStat[2,,])
)

#fit GEV to the largest, second largest order stats
MonthlyOrderStat.model <- rlarg.fit(aixMonthlyOrderStat.sample)

tau<-function(x, mu, sigma, xi){ #tau of x with parameters mu, sigma, xi
  return ((1+xi*((x-mu)/sigma))^(-1/xi))
}

H <- function(x,y){ # joint CDF of largest and second largest
  return (exp(-tau(y, mu, sigma, xi)))(1+tau(y, mu, sigma, xi)-tau(x, mu, sigma, xi))
}

Fx <- function(x, mu, sigma, xi){ #marginal for max
  return (exp(-tau(x, mu, sigma, xi)))
}

Fy <- function(y, mu, sigma,xi){ #marginal for 2nd largest
  return ((exp(-tau(y,mu, sigma,xi)))*(1+tau(y,mu, sigma,xi)))
}

#model based copula sample 
monthlyModelCopula.sample<-matrix(nrow=480,ncol=2)

for (i in 1:480){
  #apply marginal CDF to each data point
  monthlyModelCopula.sample[i,1]=Fx(aixMonthlyOrderStat.sample[i,1], MonthlyOrderStat.model$mle[1], MonthlyOrderStat.model$mle[2],MonthlyOrderStat.model$mle[3])
  monthlyModelCopula.sample[i,2]=Fy(aixMonthlyOrderStat.sample[i,2], MonthlyOrderStat.model$mle[1],MonthlyOrderStat.model$mle[2],MonthlyOrderStat.model$mle[3])
}

plot(monthlyModelCopula.sample,col="red",pch=16,xlab=expression(hat(F)[X](x[i])), ylab=expression(hat(F)[Y](y[i])))
#title("Model Based Copula Sample from Monthly Block Maxima")

#empirical copula sample
monthlyEmpCopula.sample <- cbind( rank(aixMonthlyOrderStat.sample[,1])/481 , rank(aixMonthlyOrderStat.sample[,2])/481)

plot(monthlyEmpCopula.sample,col="blue",pch=16,xlab=expression(R[x[i]]/n+1), ylab=expression(R[y[i]]/n+1))
#title("Empirical Copula Sample from Monthly Block Maxima")


# theoretical copula sample
set.seed(123)
MonthlyTheory_copula.sample <- matrix(nrow=480, ncol=2)
for (i in 1:480){
  MonthlyTheory_copula.sample[i,1] = runif(1, min=0, max=1)
  w<-runif(1,min=0,max=1)
  MonthlyTheory_copula.sample[i,2] = ((MonthlyTheory_copula.sample[i,1])*w*(1-log((MonthlyTheory_copula.sample[i,1])*w)))
}
plot(MonthlyTheory_copula.sample,pch=16,col="#82f251",xlab="U",ylab="V") 
curve(x*(1-log(x)), col="black", from=0, to=1, add=TRUE,lwd=2)

#plot all 3 to compare them
plot(monthlyModelCopula.sample, col="red", pch=16,xlab="U",ylab="V") #model based copula sample
points(monthlyEmpCopula.sample, col="blue",pch=16) #empirical copula 
points(MonthlyTheory_copula.sample, col="#82f251",pch=16) #theoretical copula sample 
curve(x*(1-log(x)), from=0, to=1, add=TRUE)
legend("topleft", legend=c("Model Based", "Empirical", "Theoretical","x(1-lnx)"), 
       col=c("red", "blue", "#82f251", "black"), pch=c(16, 16,16,NA), lwd=c(NA,NA,NA, 2), bty="n")



# testing my copula model (ie. does my data follow this model?)
# done via Rosenblatt bootstrapping 
psi<-function(u,w){
  return(u*w*(1-log(u*w)))
}
#numerical inverse of psi on region [0, 1/u]
inverse_psi <- function(row){
  root <- uniroot(function(w) psi(row[1],w) - row[2], lower=0.00000000000000000000001, upper=1/row[1])$root
  return(root)
}

Cpartial <- function(row){
  if (row[2] >= row[1]*(1-log(row[1]))){
    return (1)
  } else{
    return (inverse_psi(row))
  }
}

rosenblatt<-function(data){
  transformed.data <- matrix(nrow=nrow(data), ncol=ncol(data))
  transformed.data[, 1] = data[,1] #first column doesnt change
  transformed.data[, 2] = apply(data, 1, Cpartial) # second column is Cpartial(u,v) of each row 
  return (transformed.data)
}
plot(rosenblatt(monthlyEmpCopula.sample),pch=16,xlab="R(U)", ylab="R(V)") #suspiciously not as uniform as I'd want ; still need to double check that my rosenbaltt is performing correctly but i did test psi inverse is correct


rosenMonthly <- rosenblatt(monthlyEmpCopula.sample) #this is R(U_hat, V_hat) for the monthly data
rosenAnnual <- rosenblatt(emp_copula) #R(U_hat, V_hat) for the annual data

#joint empirical CDF 
D <- function(s,t,data){
  return(mean(data[,1] <= s & data[,2]<= t))
}


testStat<-function(data){
  stat <- 0
  for (i in 1:nrow(data)){
    stat <- stat + (D(data[i,1],data[i,2],data) - data[i,1]*data[i,2])^2
  }
  return(stat)
}

#bootstrap for monthly 

testStat_C_emp <- testStat(rosenMonthly) # want to see if this value is "too large"
testStat_C_modelBased <- testStat(rosenblatt(monthlyModelCopula.sample))
bootstraps <- numeric(1000)

for (k in 1:1000){
  hyp_cop <- matrix(nrow=480, ncol=2)
  for (i in 1:480){
    hyp_cop[i,1] = runif(1,min=0,max=1)
    w<-runif(1,min=0,max=1)
    hyp_cop[i,2] = ((hyp_cop[i,1])*w*(1-log((hyp_cop[i,1])*w)))
  }
  rosen_hyp <- rosenblatt(hyp_cop)
  bootstraps[k] <- testStat(rosen_hyp)
}

p <- mean(bootstraps > testStat_C_emp)
# p value is 0???? 
