######## Linear Regression with Gibbs Sampling ########

# yi = B0 + B1*xi + epsi,  epsi~N(0,sd^2)
# Yi ~ N(B0 + B1*xi, sd^2)

library(invgamma)

# Simulate Data:
set.seed(55)

n <- 200 # Number of data points
xvector <- rnorm(n, mean = 10, sd = 4) # x points
xvector <- sort(xvector) # Sort the values of the vector

beta0 <- 10 # beta0: intersection
beta1 <- 2.421 # beta1: slope
sd2 <- 6 # Variance of noise distribution
eps <- rnorm(n, mean = 0, sd = sqrt(sd2)) # Noise added to data

yvector <- beta0 + beta1*xvector + eps

# Plot data
plot(xvector, yvector, main = 'Plot of Data',
     xlab = "x", ylab="Y")

## Linear Regression using the least squares method:

LSM <- lm(yvector ~ xvector)
summary(LSM)

ypred <- predict(LSM)
length(ypred)

# Plot data
plot(xvector, yvector, main = 'Plot of Data',
     xlab = "x", ylab="Y")
lines(xvector, ypred, col='red',lty=1)
# Create legend
legend("topleft", legend=c("Data","Least Squares Method"),
       col=c("black","red"), lty = c(1,1), cex=0.8)

## Linear Regression using Gibbs Sampling

# Prior distributions:
# beta0 ~ N(mu0, sd02)
mu0 <- 0
sd02 <- 1000

# beta1 ~ N(mu1, sd12)
mu1 <- 0
sd12 <- 1000

# sd2 ~ IG(a,b)
a <- 0.5
b <- 1

## Gibbs Algorithm:
num_iter <- 30000 # number of iterations

# initial sampling values
beta0_samp <- 0
beta1_samp <- 10
var_samp <- 15

# Arrays to store the sample values
sample_store <- matrix(0, num_iter, 3)
sample.plot <- matrix(0, num_iter/100, 3)

for(i in 1:num_iter){
  # Posteriori values for conditional probability of 
  # beta0|beta1,sd2,Y ~ N(mu_post0, std_post0^2)
  numerator0 <- sum(yvector - beta1_samp*xvector)/var_samp + mu0/sd02 # numerator of expression
  denom0 <- n/var_samp + 1/sd02 #denominator of expression
  mu_post0 <- numerator0/denom0 # mean of beta0
  var_post0 <- 1/denom0 # variance of beta0
  
  # We sample beta0 assuming a known beta1 and variance:
  beta0_samp <- rnorm(1, mean = mu_post0, sd = sqrt(var_post0))
  
  # Posteriori values for conditional probability of 
  # beta1|beta0,sd2,Y ~ N(mu_post1, std_post1^2)
  numerator1 <- sum((yvector - beta0_samp)*xvector)/var_samp + mu1/sd12 # numerator of expression
  denom1 <- sum(xvector^2)/var_samp + 1/sd12 #denominator of expression
  mu_post1 <- numerator1/denom1 # mean of beta0
  var_post1 <- 1/denom1 # variance of beta0
  
  # We sample beta1 assuming a known beta0 and variance:
  beta1_samp <- rnorm(1, mean = mu_post1, sd = sqrt(var_post1))
  
  # Posteriori values for conditional probability of 
  # sd2|beta0,beta1,Y ~ IG(a_post, b_post)
  a_post <- n/2 + a
  b_post <- sum((yvector - beta0_samp - beta1_samp*xvector)^2)/2 + b
  
  # We sample variance assuming a known beta0 and beta1:
  var_samp <- rinvgamma(1, shape = a_post, rate = b_post)
  
  # Store the sampling values:
  sample_store[i,] <- c(beta0_samp, beta1_samp, var_samp) 
  
  if (i%%100 == 0){ # Store 300 values
    #store the sampled parameters
    sample.plot[i/100,] <- c(beta0_samp, beta1_samp, var_samp) 
  }
}

iters = seq(1,num_iter,100)
plot(iters, sample.plot[1:300,1], type='l', main = 'Sample Values of Beta0',
     xlab = "Iteration Number", ylab="Beta0")
plot(iters, sample.plot[1:300,2], type='l', main = 'Sample Values of Beta1',
     xlab = "Iteration Number", ylab="Beta1")
plot(iters, sample.plot[1:300,3], type='l', main = 'Sample Values of Variance',
     xlab = "Iteration Number", ylab="Variance")

ypredG <- sample_store[num_iter,1] + sample_store[num_iter,2]*xvector

# Plot data
plot(xvector, yvector, main = 'Plot of Data',
     xlab = "x", ylab="Y")
lines(xvector, ypred, col='red',lty=1)
lines(xvector, ypredG, col='blue',lty=1)
# Create legend
legend("topleft", legend=c("Data","Least Squares Method","Gibbs Sampling"),
       col=c("black","red","blue"), lty = c(1,1,1), cex=0.8)

