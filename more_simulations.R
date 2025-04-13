# trying some simulations to see if as n increases for iid unif we draw near the curve


n <- 1000  # this is what will go to infinity, we start with just 1 unif, then 2 unifs, then 3 unifs, etc. and rank them
num_layers <- 5  # Number of layers

# Initialize a 3D array (100 rows, 2 columns, n layers)
simulations <- array(0, dim = c(1000, 2, n))

set.seed(123) 

for (layer in 1:n) {
  for (i in 1:1000) {
    sim_values <- runif(layer) # so the nth layer we take n runifs to rank 
    
    sorted_values <- sort(sim_values, decreasing = TRUE)
    
    # Store the max in first column, second max in second column
    simulations[i, 1, layer] <- sorted_values[1]  # Maximum
    simulations[i, 2, layer] <- sorted_values[2]  # Second highest
  }
}

#so in theory, simulations[,,2] should be the same as we had done before with the 2 runifs
library(copula)

par(mar = c(5, 6, 4, 2)) 
plot(pobs(simulations[,,2]), col="blue", pch=16,xlab=expression("Largest"),ylab=expression("Second Largest"))
#curve( 1-(1-sqrt(x))^2, from=0,to=1,add=TRUE)
#curve(3*x^(2/3)-2*x, from=0,to=1,add=TRUE,col="blue")
curve(x*(1-log(x)),from=0,to=1,add=TRUE,col="red") # the limit curve
points(pobs(simulations[,,100]), col="#82f251", pch=16)
legend("topleft", legend=c("n=2", "n=100", "x(1-lnx)"), 
       col=c("blue", "#82f251", "red"), pch=c(16, 16,NA), lwd=c(NA,NA, 2), bty="n")
#points(pobs(simulations[,,10]), col="green", pch=16)
# we want to see as n increases if they draw near the curve
# I think they do!!!!


# actually take a theoretical copula sample
#demonstrating convergence behavior as n increases

n<-100
set.seed(123)
theory_copula.sample <- matrix(nrow=1000, ncol=2) # copula sample of size 1000
for (i in 1:1000){
  theory_copula.sample[i,1] = runif(1, min=0, max=1)
  u<-theory_copula.sample[i,1]
  w<-runif(1,min=0,max=1)
  theory_copula.sample[i,2] = n*w*u^((n-1)/n) - (n-1)*w^(n/(n-1))*u
}

plot(theory_copula.sample, col="#82f251",pch=16,xlab="U", ylab="V")
curve(n*x^((n-1) / n) - (n-1)*x, from=0,to=1, add=TRUE,col="black")
curve(x*(1-log(x)), from=0,to=1, add=TRUE,col="red",)
legend("topleft", legend=c("n=100", expression(phi(x)), "x(1-lnx)"), 
       col=c("#82f251", "black", "red"), pch=c(16, NA,NA), lwd=c(NA,2, 2), bty="n")

