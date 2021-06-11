######## Multilinear Regression with Gibbs Sampling ########

# yi = B0 + B1*x1i + B2*x2i + ... + Bp*xpi + epsi,  epsi~N(0,sd^2)
# Yi ~ N(B0 + B1*xi + B2*x2i + ... + Bp*xpi, sd^2)


library(invgamma)

###### Simulate Data ######
set.seed(540)

dimension <- 5 # Work with 5 x vectors
n <- 200 # number of data points for each x vector

x_matrix = matrix(0, n, dimension)

# Create the xvectors and the xmatrix:
for (i in 1:dimension){
  minval = runif(1, min = -10,max = 10) # random min value
  maxval = runif(1, min = 20, max = 50) # random max value
  xvec <- runif(n, min = minval, max = maxval) # sample n values from unif. dist.
  xvec <- sort(xvec) # Sort the values of the vector 
  x_matrix[1:n,i] = xvec # store values in matrix
}

## Simulate beta values
beta0 <- 4.543 # beta0: intersection
beta_vec <- rnorm(dimension, mean = 0, sd = 12) # beta values for each vector x: slopes

## Simulate noise
sd2 <- 20 # Variance of noise distribution
eps <- rnorm(n, mean = 0, sd = sqrt(sd2)) # Noise added to data

## Calculate the y_vector
y_vec <- beta0 + rowSums(t(t(x_matrix) * beta_vec)) + eps

# Plot data
plot(x_matrix[1:n,5], y_vec, main = 'Plot of Data',
     xlab = "x", ylab="Y")


###### Linear Regression using Gibbs Sampling ######

# Prior distributions:
# beta_i ~ N(mu = 0, sd = 1000)
mu_vec <- rep(0,dimension + 1)
sd_vec <- rep(1000,dimension + 1)

# sd2 ~ IG(a,b)
a <- 0.5
b <- 1

## Gibbs Algorithm parameters:
num_iter <- 50000 # number of iterations

# initial sampling values
beta_samp <- rep(0,dimension+1)
# beta_samp <- c(0,1,2,3,4,5)
var_samp <- 15

## Gibbs Algorithm:
GLR_model <- function(y_data, x_m_data, betas, variance_samp,  mean_vals, var_vals, 
                      a_param, b_param, num_iterations){
  
  ### Gibbs Linear Regression Model ###
  # - y_data, x_m_data: y data and x vectors to fit. x_m_data shape is (num_data_points, num_vectors)
  # - betas: beta initial sampling values. shape: (num_betas)
  # - variance_samp: initial sampling variance
  # - mean_vals, var_vals: mean and variance values of a normal prior distribution for the beta parameters
  # - a_param, b_param: a and b parameters of inverse gamma prior distribution for the variance
  # - num_iterations: number of iterations in the Gibbs algorithm
  
  dim_beta <- length(betas) # number of beta parameters
  n_size <- length(y_data) # size of data

  # Arrays to store the sample values
  sample.store <- matrix(0, num_iterations, dim_beta + 1)
  sample.plot <- matrix(0, num_iterations/100, dim_beta + 1)

  for(i in 1:num_iterations){
    # Posteriori values for conditional probability of
    # beta0|beta1,sd2,Y ~ N(mu_post0, std_post0^2)

    R0i = y_data - rowSums(t(t(x_m_data) * betas[2:(dim_beta)])) # calculate the R0 term
    numerator0 <- sum(R0i)/var_samp + mu_vec[1]/sd_vec[1] # numerator of expression
    denom0 <- n_size/var_samp + 1/sd_vec[1] #denominator of expression

    # We sample beta0 assuming all known beta_i and variance:
    beta0_samp <- rnorm(1, mean = numerator0/denom0, sd = 1/sqrt(denom0))
    betas[1] = beta0_samp

    # Posteriori values for conditional probability of
    # betak|b0,b1,...,b_{k-1},b_{k+1},...,bp,sd2,Y ~ N(mu_post, std_post^2)
    for (k in 2:dim_beta){
      dropB <- betas[-k][2:(dim_beta - 1)] # vector of betas without the kth value
      dropx <- x_m_data[1:n_size,-(k-1)] # matrix of x without the kth value
      Rk <- y_data - betas[1] - rowSums(t(t(dropx) * dropB)) # Calculate the Rk term
      numeratork <- sum(Rk * x_m_data[1:n_size,(k-1)])/var_samp + mean_vals[k]/var_vals[k] # numerator of expression
      denomk <- sum(x_m_data[1:n_size,(k-1)]^2)/var_samp + 1/var_vals[k] #denominator of expression

      # We sample beta1 assuming a known beta0 and variance:
      beta_k_samp <- rnorm(1, mean = numeratork/denomk, sd = 1/sqrt(denomk))
      betas[k] <- beta_k_samp
      }

    # Posteriori values for conditional probability of
    # var_samp|beta_vec,Y ~ IG(a_post, b_post)
    a_post <- n_size/2 + a_param
    R0i2 = y_data - betas[1] - rowSums(t(t(x_m_data) * betas[2:(dim_beta)]))# Calculate the Rk2 term
    b_post <- sum(R0i2^2)/2 + b_param

    # We sample variance assuming a known beta0 and beta1:
    var_samp <- rinvgamma(1, shape = a_post, rate = b_post)
    
    # Store the sampling values:
    sample.store[i,] <- c(betas,var_samp)
    if (i%%100 == 0){ # Store 300 values
      #store the sampled parameters
      sample.plot[(i/100),] <- c(betas,var_samp)
      }
  }
  newlist <- list(sample.plot, sample.store)
  return(newlist)}

results <- GLR_model(y_vec, x_matrix, beta_samp, var_samp, mu_vec, sd_vec, a, b, num_iter)

results.plot = results[[1]]
results.store = results[[2]]

iters = seq(1,num_iter,100)
par(mfrow=c(3,3))
plot(iters, results.plot[1:(num_iter/100),1], type='l', main = 'Sample Values of Beta0',
     xlab = "Iteration Number", ylab="Beta 0")
plot(iters, results.plot[1:(num_iter/100),2], type='l', main = 'Sample Values of Beta1',
     xlab = "Iteration Number", ylab="Beta 1")
plot(iters, results.plot[1:(num_iter/100),3], type='l', main = 'Sample Values of Beta 2',
     xlab = "Iteration Number", ylab="Beta 2")
plot(iters, results.plot[1:(num_iter/100),4], type='l', main = 'Sample Values of Beta 3',
     xlab = "Iteration Number", ylab="Beta 3")
plot(iters, results.plot[1:(num_iter/100),5], type='l', main = 'Sample Values of Beta 4',
     xlab = "Iteration Number", ylab="Beta 4")
plot(iters, results.plot[1:(num_iter/100),6], type='l', main = 'Sample Values of Beta 5',
     xlab = "Iteration Number", ylab="Beta 5")
plot(iters, results.plot[1:(num_iter/100),7], type='l', main = 'Sample Values of Variance',
     xlab = "Iteration Number", ylab="Variance")


print("#### Results ####")
print('1) Beta 0 ')
sprintf('Expected value: %f',beta0)
sprintf('Obtained value: %f',mean(results.store[30001:num_iter,1]))
print('2) Beta 1 ')
sprintf('Expected value: %f',beta_vec[1])
sprintf('Obtained value: %f',mean(results.store[30001:num_iter,2]))
print('3) Beta 2 ')
sprintf('Expected value: %f',beta_vec[2])
sprintf('Obtained value: %f',mean(results.store[30001:num_iter,3]))
print('4) Beta 3 ')
sprintf('Expected value: %f',beta_vec[3])
sprintf('Obtained value: %f',mean(results.store[30001:num_iter,4]))
print('5) Beta 4 ')
sprintf('Expected value: %f',beta_vec[4])
sprintf('Obtained value: %f',mean(results.store[30001:num_iter,5]))
print('6) Beta 5 ')
sprintf('Expected value: %f',beta_vec[5])
sprintf('Obtained value: %f',mean(results.store[30001:num_iter,6]))
print('6) Variance ')
sprintf('Expected value: %f',sd2)
sprintf('Obtained value: %f',mean(results.store[30001:num_iter,7]))


