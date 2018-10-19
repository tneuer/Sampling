try(dev.off(), silent = T)
rm(list=ls())

source("/home/thomas/Uni/Algorithmen/Sampler/A-R-Method.R")

start = Sys.time()
RMS = function(joint, renvelope, fenvelope, initial = 0, n = 10000, sim = F) {
  sample = matrix(0, nrow = n, ncol = length(initial))
  prev = initial
  counter1 = counter2 = 0
  
  while (counter1 < n) {
    counter2 = counter2 + 1
    z = renvelope(1)
    u = runif(1)
    
    if (u < joint(z) / fenvelope(z)) {
      u = runif(1)
      counter1 = counter1 + 1
      
      if (joint(prev) < fenvelope(prev)) C1 = 1
      else C1 = 0
      
      if (joint(z) < fenvelope(z)) C2 = 1
      else C2 = 0
    
      if (C1 == 1) alpha = 1 
      else if (C1 == 0 && C2 == 1) alpha = fenvelope(prev)/ joint(prev)
      else alpha = joint(z) * fenvelope(prev) / (joint(prev) * fenvelope(z))
      
      
      if (u < alpha) {
        sample[counter1, ] = z
        prev = z
      }
      else {
        sample[counter1, ] = prev
      }
    }
  }
  acc = n/counter2
  return(list(sample = sample, acc = acc))
}


if(T) {
  x = seq(-50, 250, by = 0.1)
  
  f1 <- function(x) {
    0.25*(dnorm(x, mean = 12, sd = 20) + dnorm(x, mean = -10, sd = 2) + dnorm(x, mean = 200, sd = 5) + dgamma(x-150, shape = 3, rate = 0.2))
  }
  
  renvelope = function(x) {
    runif(x, min = -50, max = 250)
  }
  
  fenvelope = function(x) {
    3*dunif(x, min = -50, max = 250)
  }
  
  sample = RMS(f1, renvelope, fenvelope, n = 30000)
  s = sample$sample
  hist(s, breaks = 150, probability = T, col = "yellow")
  lines(x, f1(x), lwd = 2, col = 2, lty = 2)
  lines(x, fenvelope(x), lty = 3, col = 4, lwd = 3)
  sample$acc
  mean(s)
  (12 - 10 + 200 + 150+3/0.2) /4
  # acf(s)
  Sys.time() - start
}
