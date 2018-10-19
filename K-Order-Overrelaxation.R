dev.off()
rm(list=ls())

Overrelaxation <- function(data, p1, p2, init, K = 20, iterations = 10000, burn_in = iterations/5, autocorr = 1) {
  theta1 = theta2 = seq(0, 0, length.out = iterations+1)
  theta1[1] = init[1]
  theta2[1] = init[2]
  return_val = seq(burn_in+1, iterations, by = autocorr)
  thetatemp = rep(0, l = K)
  for (i in 2:iterations+1) {
    for (j in 1:K) {thetatemp[j] = p1(data, theta2[i-1])}
    thetasort = sort(thetatemp)
    theta1[i] = thetasort[K+1-which(thetasort == thetatemp[1])]
    for (j in 1:K) {thetatemp[j] = p2(data, theta1[i-1])}
    thetasort = sort(thetatemp)
    theta2[i] = thetasort[K+1-which(thetasort == thetatemp[1])]
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

analyze <- function(theta, hist = T, trace = T, acf = T, simulate = F, nsim = 200) {
  print(summary(theta))
  print(paste("Standard deviation: ", sd(theta), "Variance: ", var(theta)))
  par(mfrow=c(3,1))
  if (hist) {
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
  if (trace) {
    plot(1:length(theta), theta, xlab = "Samplenumber", ylab = "Samplevalue", type = "b", main = "Traceplot")
  }
  if (acf) {
    acf(theta)
  }
}

PleaseWork <- Overrelaxation(data = y, p1 = cond_mu, p2 = cond_sigma, init = c(0,1), iterations = 10000, burn_in = 4000)
t1 = PleaseWork$theta1
t2 = PleaseWork$theta2
analyze(t1, simulate = F, nsim = 50)
analyze(t2)
analyze(1/t2)
analyze(sqrt(1/t2))