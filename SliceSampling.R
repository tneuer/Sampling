# try(dev.off(), silent = T)
# rm(list=ls())

start = Sys.time()
SliceSampler = function(joint, initial, n = 3, stepwidth= 1, sim = F, pause = 1, roi = c(-10,10)) {
  sample = rep(0, l = n + 1)
  sample[1] = initial
  feval = 0
  if (sim == T) {
    x = seq(roi[1], roi[2], by = 0.1)
    plot(x, joint(x), type = "l")
    for (i in 2:(n+1)) {
      y.unif = runif(1, 0, joint(sample[i-1]))
      k = 1
      points(sample[i-1], y.unif, col = 2)
      x.range = c(sample[i-1] - stepwidth, sample[i-1] + stepwidth)
      lines(x.range, c(y.unif,y.unif), col = 4)
      while( y.unif < joint(x.range[1]) ) {
        feval = feval + 1
        x.range[1] = x.range[1] - stepwidth
        Sys.sleep(pause)
        lines(x.range, c(y.unif, y.unif), col = 4, lty = 2)
      }
      while( y.unif < joint(x.range[2]) ) {
        feval = feval + 1
        x.range[2] = x.range[2] + stepwidth
        Sys.sleep(pause)
        lines(x.range, c(y.unif, y.unif), col = 4, lty = 2)
      }
      proposal = runif(1, x.range[1], x.range[2])
      Sys.sleep(pause)
      points(proposal, y.unif, col = 3, type = "o")
      while (joint(proposal) < y.unif) {
        feval = feval + 1
        if (proposal < sample[i-1]) {x.range[1] = proposal}
        else {x.range[2] = proposal}
        Sys.sleep(pause)
        points(proposal, y.unif, col = 7, type = "o")
        proposal = runif(1, x.range[1], x.range[2])
        Sys.sleep(pause)
        points(proposal, y.unif, col = 5, type = "o")
      }
      sample[i] = proposal
      Sys.sleep(pause)
      plot(x, joint(x), type = "l")
      points(c(sample[i],sample[i-1]), c(y.unif, y.unif), col = c(3, 2), type = "o")
      Sys.sleep(pause)
      lines(x.range, c(y.unif, y.unif), col = 8)
      rug(sample)
      Sys.sleep(pause)
      lines(c(proposal, proposal), c(0, 1), col = 8)
      Sys.sleep(pause)
    }
  }
  else {
    for (i in 2:(n+1)) {
      y.unif = runif(1, 0, joint(sample[i-1]))
      k = 1
      x.range = c(sample[i-1] - stepwidth, sample[i-1] + stepwidth)
      while( y.unif < joint(x.range[1]) ) {
        feval = feval + 1
        x.range[1] = x.range[1] - stepwidth
      }
      while( y.unif < joint(x.range[2]) ) {
        feval = feval + 1
        x.range[2] = x.range[2] + stepwidth
      }
      proposal = runif(1, x.range[1], x.range[2])
      while (joint(proposal) < y.unif) {
        feval = feval + 1
        if (proposal < sample[i-1]) {x.range[1] = proposal}
        else {x.range[2] = proposal}
        proposal = runif(1, x.range[1], x.range[2])
      }
      sample[i] = proposal
    }
    
  }
  return(list(sample = sample, feval = feval))
}

if (F) {
x = seq(-50, 250, by = 0.1)
fx = function(x) 0.25*(dnorm(x, mean = 12, sd = 20) + dnorm(x, mean = -10, sd = 3) + dnorm(x, mean = 130, sd = 5) + dgamma(x-80, shape = 3, rate = 0.2))

plot(x, fx(x), type = "l")

sample1 = SliceSampler(fx, initial = 5, n = 30000, sim = T, stepwidth = 30, pause = 0.7, roi = c(-50, 150))
s = sample1$sample
f = sample1$feval
hist(s, breaks = 200, prob = T, col = "yellow")
lines(x, fx(x), col = 2, lwd = 2, lty = 2)
mean(s)
(12 - 10 + 130 + 95) /4
f
#acf(s)
Sys.time() - start
}

