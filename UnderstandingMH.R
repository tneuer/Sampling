rm(list=ls())
try(dev.off(), silent = T)
library(MASS)
library(mvtnorm)
source("/home/thomas/Uni/Nebenfach/SelectedTopics/Programme/RMS.R")

mu = c(1,2)
sigma = matrix(c(1,0.9,0.9,1), byrow = T, nrow = 2)

proposal1 <- function(value){
  new1 <- runif(1, min = value[1]-0.75, max = value[1]+0.75)
  new2 <- runif(1, min = value[2]-1, max = value[2]+1)
  return(c(new1, new2))
}

proposal2 <- function(value){
  return(rnorm(2,mean = value, sd= c(0.6, 0.4)))
}

pseudo_fenvelope <- function(value){
  return(0.9 * dmvnorm(value, mean = mu, sigma = matrix(c(2,0,0,2), nrow = 2)))
}

pseudo_renvelope <- function(n) {
  return(rmvnorm(n, mean = mu, sigma = matrix(c(2,0,0,2), nrow = 2)))
}

pseudo_joint <- function(x) {
  return(dmvnorm(x, mean = mu, sigma = sigma))
}

autoregressive <- function(value){
  z = runif(2, min = -1, max = 1)
  y = mu - (value - mu) + z
  return(y)
}

run_metropolis_MCMC <- function(startvalue, iterations, proposalfunction){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  counter = 0
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    alpha = dmvnorm(proposal, mean = mu, sigma = sigma) / dmvnorm(chain[i, ], mean = mu, sigma = sigma)
    if (runif(1) < alpha){
      chain[i+1,] = proposal
      counter = counter + 1
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  acc = counter/iterations
  return(list(sample = chain, acc = acc))
}

initial = c(0, 0)
iter = 8000

x = seq(-5, 5, by = 0.05)
y = seq(-1 ,6, by = 0.05)
z = dmvnorm(expand.grid(x,y), mu, sigma)
z = matrix(z, nrow = length(x))

burnIn = iter/4

par(mfrow = c(3,2))

start = Sys.time()
chain1 = run_metropolis_MCMC(startvalue = initial, iterations = iter, proposalfunction = proposal1)
acc = chain1$acc
chain1 = chain1$sample
contour(x,y,z, xlim = c(-2,4), col = 2, main = "Random Walk uniform")
points(chain1[burnIn:iter, 1], chain1[burnIn:iter, 2])
acf(chain1[burnIn:iter, 1], main = paste("Time needed: ", round(Sys.time() - start, 4), ", Acceptance: ", acc))


start = Sys.time()
chain2 = run_metropolis_MCMC(startvalue = initial, iterations = iter, proposalfunction = proposal2)
acc = chain2$acc
chain2 = chain2$sample
contour(x,y,z, xlim = c(-2,4), col = 2, main = "Random Walk normal")
points(chain2[burnIn:iter, 1], chain2[burnIn:iter, 2])
acf(chain2[burnIn:iter, 1], main = paste("Time needed: ", round(Sys.time() - start, 4), ", Acceptance: ", acc))


# start = Sys.time()
# chain3 = RMS(pseudo_joint, pseudo_renvelope, pseudo_fenvelope, initial, iter)
# acc = chain3$acc
# chain3 = chain3$sample
# contour(x,y,z, xlim = c(-2,4), col = 2, main = "Pseudorejection")
# points(chain3[burnIn:iter, 1], chain3[burnIn:iter, 2])
# acf(chain3[burnIn:iter, 1], main = paste("Time needed: ", round(Sys.time() - start, 4), ", Acceptance: ", acc))


start = Sys.time()
chain4 = run_metropolis_MCMC(startvalue = initial, iterations = iter, proposalfunction = autoregressive)
acc = chain4$acc
chain4 = chain4$sample
contour(x,y,z, xlim = c(-2,4), col = 2, main = "Autoregressive")
points(chain4[burnIn:iter, 1], chain4[burnIn:iter, 2])
acf(chain4[burnIn:iter, 1], main = paste("Time: ", round(Sys.time() - start, 4), ", Acceptance: ", acc))



# par(mfrow=c(2,1))
# sample_R = mvrnorm(6000, mu, sigma)
# contour(x,y,z, xlim = c(-2,4), col = 2, main = "Cholesky")
# points(sample_R)
# 
# contour(x,y,z, xlim = c(-2,4), col = 2, main = "Autoregression")
# points(chain4[burnIn:iter, 1], chain4[burnIn:iter, 2])