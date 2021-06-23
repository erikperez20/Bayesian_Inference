### Data: Yi ~ N(mu,sd^2)
### mu ~ N(mu0, sd0^2) Apriori 1
### sd^2 ~ InvGamma(a,b) Apriori 2

library(invgamma)
library(rjags)

### Calculate the posterior distribution samples with the Gibbs Sampling Method and
### using the RJAGS library

## Simulate data:
n = 1000 # Side of data
mu = 3.5 # Mean of data
std = 1.5 # Standard deviation of data
data = rnorm(n = n, mean = mu, sd = std) # simulate data comming from a normal dist
meanY = mean(data)
hist(data,breaks = 15)

# JAGS requiere que todos los datos estén en una lista
data <- list(X=data,n=n)

model_string <- textConnection("model{
for(i in 1:n){
X[i] ~ dnorm(mu, tau) #tau = 1/sigma^2
}
# A Prioris
tau ~ dgamma(0.5, 1)
sigma2 <- 1/tau
mu ~ dnorm(0, 0.001)
}")

inits <- list(mu=0,tau=1)
model <- jags.model(model_string,data = data, inits=inits, n.chains=2,quiet=TRUE)
update(model, 10000, progress.bar="none")
params <- c("mu","sigma2")
samples <- coda.samples(model,
                        variable.names=params,
                        n.iter=30000, progress.bar="none")
summary(samples)
plot(samples)












