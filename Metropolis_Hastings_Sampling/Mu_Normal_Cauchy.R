### Data: Yi ~ N(mu,1)
### mu ~ Cauchy(0, 1) A priori

### Calculate the posterior distribution samples with the Metropolis-Hasting
### Sampling Method

## Simulate data:

set.seed(22)
Y <- dnorm(200,1,1)

# Function to compute the posterior distribution:
# Likelihood: Y[i] ~ N(theta,1)
# Prior:      theta ~ Cauchy(0,1)
posteriori <- function(theta, y, param1 = 0, param2 = 1){
  priori <- dcauchy(theta, location = param1, scale = param2)
  likelihood <- dnorm(y, mean = theta, sd = 1)
  return(priori*likelihood)}

## Independent proposal distribution:
independent_proposal <- function(theta, y){
  dist_func = dnorm(theta, mean= y, sd = 1)
  return(dist_func)}
## Candidate generated from the independent proposal
candidate_independent <- function(y){
  candidate = rnorm(1, mean = y, sd = 1)
  return(candidate)}

## Random Walk proposal distribution:
random_walk_proposal <- function(theta,theta_mean,sj2){
  dist_func = dnorm(theta, mean = theta_mean, sd = sqrt(sj2))
  return(dist_func)}
## Candidate generated from the random walk proposal dist.
candidate_random_walk <- function(theta, sj2){
  candidate = rnorm(1, mean = theta, sd = sqrt(sj2))
  return(candidate)}


###### Metropolis-Hasting Algorithm: ######

num_iter = 30000 # number of iterations

# initial sampling value
theta <- 0.5
# tuning parameter for a random walk proposed distribution
s2 = 0.1^2 


# Metropolis-Hasting sampling function:

MCMC_function <- function(theta_init, ydata, post_dist, proposed_dist,
                          num.iter, tuning_param){
  theta_val = theta_init
  
  # Arrays to store the sample values
  th_samples = matrix(0,num.iter)
  th_s.plot = matrix(0,num.iter/100)
  
  # Start iteration
  for(iter in 1:num.iter){
    if(proposed_dist == "Independent"){
      
      # Calculate candidate
      can <- candidate_independent(ydata)
      
      # Calculate the R term
      num1 <- post_dist(can, ydata) # first term pi(th*|y)
      denom1 <- post_dist(theta_val, ydata) # second term pi(th^(t-1)|y)
      num2 <- independent_proposal(theta_val, ydata) # 3th term pi(th^(t-1)|th*)
      denom2 <- independent_proposal(can, ydata) #4th term pi(th*|th^(t-1))
      }
    
    else if(proposed_dist == "Random Walk"){
      
      # Calculate candidate
      can <- candidate_random_walk(theta_val, tuning_param)

      # Calculate the R term
      num1 <- post_dist(can, ydata) # first term pi(th*|y)
      denom1 <- post_dist(theta_val, ydata) # second term pi(th^(t-1)|y)
      num2 <- random_walk_proposal(theta_val, can, tuning_param) # 3th term pi(th^(t-1)|th*)
      denom2 <- random_walk_proposal(can, theta_val, tuning_param) #4th term pi(th*|th^(t-1))
      }
    
    # Acceptance rate factor
    AcceptRate <- num1*num2/(denom1*denom2)
    # R factor: min between the acceptance factor and 1
    R <- min(1, AcceptRate)
    
    # Condition of acceptance of the theta candidate
    if (R == 1){
      theta_val <- can}
    else if (R < 1){
      unif <- runif(1)
      if (unif < R){
        theta_val = can}
    }
    
    # Store values in array
    th_samples[iter] <- theta_val
    if (iter%%100 == 0){ # Store 300 values
      th_s.plot[iter/100] <- theta_val} #store the mu sampling 
  }
  return(list(th_samples,th_s.plot))
}

# Compute the posterior on a grid for reference
theta_grid <- seq(-3,3,length=100)
dense      <- rep(0,100)
for(i in 1:100){
  dense[i]<- posteriori(theta_grid[i],Y)
}

### Random Walk proposed distribution

samples = MCMC_function(theta, Y, posteriori, "Random Walk", num_iter, s2)

results.store = samples[[1]]
results.plot = samples[[2]]

iters = seq(1,num_iter,100)
plot(iters, results.plot, type='l', main = 'Sample Values of theta using Random Walk proposed distribution',
     xlab = "Iteration Number", ylab="Theta")

summary(results.store[10000:num_iter])

hist(results.store[10000:num_iter],breaks=200,main="Posterior of theta distribution",freq = FALSE)
lines(theta_grid, dense*5,lty=1)

### Independent proposed distribution 

samples2 = MCMC_function(theta, Y, posteriori, "Independent", num_iter, s2)

results.store2 = samples2[[1]]
results.plot2 = samples2[[2]]

plot(iters, results.plot2, type='l', main = 'Sample Values of theta using Independent proposed distribution',
     xlab = "Iteration Number", ylab="Theta")

summary(results.store2[10000:num_iter])

hist(results.store2[10000:num_iter],breaks=200,main="Posterior of theta distribution",freq = FALSE)
lines(theta_grid, dense*5,col='red',lty=1)









