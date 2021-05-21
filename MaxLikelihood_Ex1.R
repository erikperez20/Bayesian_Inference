library(mvtnorm)

theta = seq(0,1,0.01)# Create sequence

n=100 # data size
sumx=60 # sum of x values
EMV=sumx/n # Maximum likelihood value

# Parameters for the beta distribution
alpha = 10 
beta = 15

likelihood = (theta^sumx)*(1-theta)^(n-sumx) # likelihood dist.
priori = dbeta(theta, alpha,beta) # prior dist.
posteriori = dbeta(theta , sumx+alpha,n-sumx+beta) # posteriori dist.

# Plot de dist. posteriori
plot(theta,posteriori,type='l') 
lines(theta,posteriori,col='black',lty=1)
# Plot de la distribution a priori
lines(theta,priori,col="blue",lty=2)
# Plot the la funcion de verosimilitud
# plot(theta,likelihood*10^293,col="green", lty=3)
lines(theta,likelihood*10^30,col="dark red", lty=1)
abline(v=EMV,col="red")

# Create legend
legend("topleft", legend=c("Posteriori", "Priori","Likelihood"),
       col=c("black", "blue","dark red"), lty = 1:2, cex=0.8)

