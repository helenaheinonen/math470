n <-2

data <- matrix(nrow=1000,ncol=2)
for(i in 1:1000){
  data[i,] <- sort(runif(n),decreasing=TRUE)[c(1,2)]
}


f<-function(u,w){
  return((w*sqrt(1-u) + 1 - sqrt(1-u))^2)
}
theory_cop <- matrix(nrow=1000, ncol=2)
for (i in 1:1000){
  theory_cop[i,2] = runif(1, min=0, max=1)
  w<-runif(1,min=0,max=1)
  theory_cop[i,1] = f(theory_cop[i,2],w)
}

#ii switched the two columns so that they would align i dont think this changes anything but just noting
# wait but thats weird that the runif follows a curve though? something is suspicious here

rfrechet <- function(n) {
  U <- runif(n)  # Generate n uniform(0,1) samples
  return(1 / (-log(U)))  # Apply inverse CDF transformation
}

# Example: Generate 1000 Fréchet samples
set.seed(123)  # For reproducibility

sample_frechet <- matrix(nrow=1000,ncol=2)
for (i in 1:1000){
  sample_frechet[i,]<-sort(rfrechet(2),decreasing=TRUE)[c(1,2)]
  }

empCop_ind_frechet<-cbind(rank(sample_frechet[,1]/1000), rank(sample_frechet[,2]/1000))

library(copula)
plot(empCop_ind_frechet, col="blue", pch=16,xlab="Maximum",ylab="Second Largest",title="Empirical copula for an iid sample of Standard Frechet") #empirical copula
#points(theory_cop,col="green",pch=16) # theoretical copula
curve(x*(1-log(x)),from=0,to=1,add=TRUE)#black curve is the limit curve 
legend("topleft", legend=c("x(1-lnx)"), 
       col=c("black"), lwd=c(2), bty="n") 

#curve( 1-(1-sqrt(x))^2, from=0,to=1,add=TRUE,col="red")

# they lie above the conjectured limit curve which is odd considering the real data points (aix) lie below?

###### try to introduce some dependence into the simulation to see if the data points are below the curve 
###### this might suggest that the reason for our data points being below the curve is due to some (time) dependence structure 
# i.e. that max and second max tend to be clustered together (maybe from the same weather pattern?)

# sampling via moving maximum process 

#X_i <- max( a/(a+1) * Z_{I-1}, 1/(a+1) Z_i)


a <- 0.5
n <- 50

# code for just one simulation

MM.sim <- function(n,a){
  Z <- -1/log(runif(n+1))
  pmax((a/(a+1))*Z[1:n],(1/(a+1))*Z[2:(n+1)])
}

X <-MM.sim(n,a) # a single iteration
plot(log(X),pch=20) #plotting on the log-scale (or take a root), otherwise nothing is visible

#adapt to create a larger time dependent block sample
timeDepSample <- matrix(nrow=1000, ncol=2)


for (i in 1:1000){
  block <- MM.sim(n,a) # creates one block of size n 
  sorted <- sort(block,decreasing=TRUE)
  timeDepSample[i,1] = sorted[1] # take its max and second max
  timeDepSample[i,2] = sorted[2]
}

plot(timeDepSample)
timeDepEmpCop<-cbind(rank(timeDepSample[,1])/1000, rank(timeDepSample[,2])/1000)
plot(timeDepEmpCop, col="blue", pch=16) #empirical copula
curve(x*(1-log(x)),from=0,to=1,add=TRUE) 

# so technically they do lie below the curve but they seem to follow some other curve? which is more straight line
# and if i make a=0 (so iid then they follow the curve super nicely!!)
# weird observation, if they are *almost* iid but not quite so like a=0.1 we get a weird empirical copula
# but it definitely starts to creep away (and below) from the curve

# i dont know why it pushes up towards the curve eventually creating a straight line 

## here's the plots all next to each other using ggplot

library(ggplot2)
library(gridExtra)
a_values <- c(0, 0.2, 0.5, 0.9)
plot_list <- list()
# Loop over different a values
for (a in a_values) {
  timeDepSample <- matrix(nrow = 1000, ncol = 2)
  
  # Generate dependent sample
  for (i in 1:1000) {
    block <- MM.sim(n, a) # creates one block of size n 
    sorted <- sort(block, decreasing = TRUE)
    timeDepSample[i, 1] <- sorted[1] # Take max
    timeDepSample[i, 2] <- sorted[2] # Take second max
  }
  
  # Compute empirical copula
  timeDepEmpCop <- cbind(rank(timeDepSample[, 1]) / 1000, rank(timeDepSample[, 2]) / 1000)
  
  # Convert to dataframe
  df <- data.frame(U1 = timeDepEmpCop[, 1], U2 = timeDepEmpCop[, 2])
  
  # Create ggplot
  p <- ggplot(df, aes(x = U1, y = U2)) +
    geom_point(color = "blue", alpha = 0.5) +
    ggtitle(paste("Empirical Copula (a =", a, ")")) +
    theme_minimal() +
    xlab("Ranked Max") + 
    ylab("Ranked 2nd Max")
  p<-p + stat_function(fun = function(x) x*(1-log(x)), color = "red", size = 1)
  # Store plot
  plot_list[[as.character(a)]] <- p
}

# Display all plots in a grid
grid.arrange(grobs = plot_list, ncol = 2)

