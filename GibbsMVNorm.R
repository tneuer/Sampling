dev.off()
rm(list=ls())

library("MASS", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library("mvtnorm", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")

#Analytic plot of multivariate normal
mu = c(0,0)
rho = 0.9
sigma = matrix(data = c(1,rho,rho,1), nrow = 2, ncol = 2)

x = seq(-3, 3, by = 0.05)
y = seq(-3 ,3, by = 0.05)
z = dmvnorm(expand.grid(x,y), mu, sigma)
z = matrix(z, nrow = length(x))


contour(x,y,z)

#Gibbs sampling implementationa
samplesize = 10000
spoint = c(-0,0)
old = c(0,0) #####
lim = 3

theta1 = c(spoint[1])
theta2 = c(spoint[2])
dev.off()
par(mfrow=c(2,2))
for (i in 1:samplesize){
  barplot(hist(theta1, plot = F, breaks = 30)$density, col = "yellow") #####
  barplot(c(0)) #######
  contour(x,y,z, main = paste(i, "/", samplesize), xlim = c(-lim,lim), ylim = c(-lim,lim)) ########
  points(theta1, theta2, col = 2) ########
  old[1] = spoint[1] #######
  old[2] = spoint[2] ########
  spoint[1] = rnorm(n = 1, mean = rho*spoint[2], sd = 1-rho**2) 
  spoint[2] = rnorm(n = 1, mean = rho*spoint[1], sd = 1-rho**2) 
  lines(c(old[1],spoint[1]), c(old[2],spoint[2]),col=4) ##########
  barplot(hist(theta2, plot = F, breaks = 30)$density, horiz = T, col = "yellow") ########
  theta1 = c(theta1, spoint[1])
  theta2 = c(theta2, spoint[2])
  Sys.sleep(0.09) ########
}

dev.off()
samplesize = 10000
spoint = c(-10,3)
theta1 = c(spoint[1])
theta2 = c(spoint[2])
for (i in 1:samplesize){
  spoint[1] = rnorm(n = 1, mean = rho*spoint[2], sd = 1-rho**2) 
  spoint[2] = rnorm(n = 1, mean = rho*spoint[1], sd = 1-rho**2) 
  theta1 = c(theta1, spoint[1])
  theta2 = c(theta2, spoint[2])
}

abscisse = seq(-5, 5, length.out = 100)
hist(theta1, breaks = 100, probability = T)
mean(theta1)
sd(theta1)
lines(density(theta1), col = 2, lwd = 3)
lines(abscisse, dnorm(abscisse), col = 4, lty = 2, lwd = 3)
legend("topright", legend = c("Sample", "True"), col = c(2,4), lty = c(1,2), lwd = 3)

hist(theta2, breaks = 100, probability = T)
mean(theta2)
sd(theta2)
lines(density(theta2), col = 2, lwd = 3)
lines(abscisse, dnorm(abscisse), col = 4, lty = 2, lwd = 3)
legend("topright", legend = c("Sample", "True"), col = c(2,4), lty = c(1,2), lwd = 3)





#####################More General Gibbs Sampler
dev.off()
rm(list=ls())

GibbsSampler <- function(data, p1, p2, init, iterations = 10000, burn_in = iterations/5, thin = 1) {
  #data: sample of the distribution one wants to approximate
  #p1,p2: conditional distributions with first argument data and second argument the other paramter
  #init: list of initial values for the parameters
  #iterations: integer number of generated tuples
  #burn_in: integer number of discarded samples for the first n events
  #autocorr: integer number to output every n-th value to reduce autocorrelation of the generated sample
  
  theta1 = theta2 = seq(0, 0, length.out = iterations+1)
  theta1[1] = init[1]
  theta2[1] = init[2]
  return_val = seq(burn_in+1, iterations, by = thin)
  for (i in 2:iterations+1) {
    theta1[i] = p1(data, theta2[i-1])
    theta2[i] = p2(data, theta1[i])
  }
  return(list(theta1 = theta1[return_val], theta2 = theta2[return_val]))
}


set.seed(44566)
y = rnorm(30, mean = 4, sd = 4)

cond_mu <- function(data, sigma2_inv) {
  new_sigma2 = 1/(0.25+30*sigma2_inv)
  new_mu = (-3/4 + sum(data)*sigma2_inv)*new_sigma2
  rand = rnorm(1, mean = new_mu, sd = sqrt(new_sigma2))
  return(rand)
}

cond_sigma <- function(data, mu) {
  new_a = 1.6+15
  new_b = 0.4 + 0.5*sum((data-mu)^2)
  rand = rgamma(1, shape = new_a, rate = new_b)
  return(rand)
}

analyze <- function(theta, hist = T, trace = T, simulate = F, nsim = 200) {
  print(summary(theta))
  print(paste("Standard deviation: ", sd(theta), "Variance: ", var(theta)))
  if (hist == T) {
    hist(theta, breaks =20) 
  }
  if (simulate == T) {
    plot(1:nsim, theta[1:nsim], type = "n", main = paste("First ", nsim," samples simulated"))
    abline(h = mean(theta), col = 2)
    for (i in 1:nsim) {
      points(i, theta[i], type = "b")
      Sys.sleep(0.22)
    }
  }
  if (trace == T) {
    plot(1:length(theta), theta, xlab = "Samplenumber", ylab = "Samplevalue", type = "b", main = "Traceplot")
  }
}

PleaseWork <- GibbsSampler(data = y, p1 = cond_mu, p2 = cond_sigma, init = c(0,1), iterations = 10000, burn_in = 4000)
t1 = PleaseWork$theta1
t2 = PleaseWork$theta2
analyze(t1, simulate = T, nsim = 50)
analyze(t2)
analyze(1/t2)
analyze(sqrt(1/t2))
