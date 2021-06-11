### Data: Yi ~ N(mu,sd^2)
### mu ~ N(mu0, sd0^2) Apriori 1
### sd^2 ~ InvGamma(a,b) Apriori 2

library(invgamma)

### Calculate the posterior distribution samples with the Gibbs Sampling Method

## Simulate data:
n = 300 # Side of data
mu = 3.5 # Mean of data
std = 1.5 # Standard deviation of data 
data = rnorm(n = n, mean = mu, sd = std) # simulate data comming from a normal dist
meanY = mean(data)
hist(data,breaks = 15)

# Priori distribution parameters for mu given a known sigma^2:
mu0 = 0 # mean
var0 = 1000 # variance

# Priori distribution parameters for sigma^2 given a known mu:
a = 0.5
b = 1

xv = seq(0,3,0.1)
sol = dinvgamma(xv, shape = a, rate = b)
plot(xv, sol, type='l', main = 'Prior Distribution of std^2',
     xlab = "x", ylab="Distribution")

## Gibbs Sampling Algorithm:

num_iter = 30000 # number of iterations

# initial sampling values
mu_samp = 0 
var_samp = 1

# Arrays to store the sample values
mu_samples = matrix(0,num_iter)
var_samples = matrix(0,num_iter)
mu_s.plot = matrix(0,num_iter/100)
var_s.plot = matrix(0,num_iter/100)

for(i in 1:num_iter){
  # Posteriori values for conditional probability of 
  # mu|sd^2,Y ~ N(mu_post, std_post^2)
  mu_post = (n*meanY*1/var_samp + mu0/var0)/(n/var_samp + 1/var0)
  var_post = 1/(n*1/var_samp + 1/var0)
  
  # We sample mu assuming a known variance:
  mu_samp = rnorm(1, mean = mu_post, sd = sqrt(var_post))
  
  # Posteriori values for conditional probability of 
  # sd^2|mu,Y ~ InvGamma(a_post, b_post)
  a_post = n/2 + a
  b_post = 1/2*sum((data-mu_samp)^2) + b
  
  # We sample var assuming a known mean:
  var_samp = rinvgamma(1, shape = a_post, rate = b_post)
  
  # Store the sampling values:
  mu_samples[i] = mu_samp #store the mu sampling
  var_samples[i] = var_samp #store the var sampling

  if (i%%100 == 0){ # Store 300 values
    mu_s.plot[i/100] = mu_samp #store the mu sampling
    var_s.plot[i/100] = var_samp #store the var sampling
  }
}

iters = seq(1,num_iter,100)
plot(iters, mu_s.plot, type='l', main = 'Sample Values of Mu',
     xlab = "Iteration Number", ylab="Mu")
plot(iters, var_s.plot, type='l', main = 'Sample Values of std^2',
     xlab = "Iteration Number", ylab="std^2")

par(mfrow=c(1,1))
plot(var_samples,mu_samples,xlab="Sigma^2",ylab="mu",main="Joint posterior distribution")
abline(mean(data),0)
abline(v=var(data))

