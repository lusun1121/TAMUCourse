setwd("D:/Google Drive/2021_SPRING/STAT611")

# Generate a random sample from N(2,1) distribution
dat_norm <- rnorm(300, 2, 1) 

# Construct the normal likelihood function
lik <- function(theta, dat){
  prod(dnorm(dat, theta, 1)) # Using dnorm function
}

# Plotting
mu <- seq(1,3, by = 0.001) # possible values of the mu parameter
plot(mu, sapply(mu, function(theta) lik(theta, dat_norm)),type = 'l', ylab = 'Likelihood')
abline(v=2, col=2)


     