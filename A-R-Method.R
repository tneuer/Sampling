try(dev.off(), silent = T)
rm(list=ls())

RS = function(joint, renvelope, fenvelope, n = 10000, sim = F, roi = c(-10,10), pause = 1) {
  if(sim) {
    x = seq(roi[1],  roi[2], length.out = 30)
    jx = joint(x)
    fx = fenvelope(x)
    ymax = max(c(jx, fx))
    plot(x, jx, col = 4, ylim = c(0, ymax), type = "l")
    lines(x, fx, col = 2, lty = 2)
  }
  sample = c()
  counter = 0
  while (length(sample) < n) {
    counter = counter + 1
    z = renvelope(1)
    u = runif(1, min = 0, max = 1)
    comp = joint(z) / (fenvelope(z))
    if (u <= comp) {
      sample = c(sample, z)
    }
    if (counter > 5) pause = 0.001
    if(sim && length(sample) != 0) {
      hist(sample, breaks = 60, probability = T, col = "yellow", ylim = c(0, ymax), xlim = c(roi[1], roi[2]), main = paste("Rate: ",round(length(sample)/counter, digits = 3),", Samples: ", length(sample)))
      if (u*fenvelope(z) <= joint(z)) {col = "green"}
      else {col = 2}
      points(z, u*fenvelope(z), col = col, pch = 19)
      lines(x, jx, col = 4)
      lines(x, fx, col = 2, lty = 2)
      Sys.sleep(pause)
    }
  }
  acc = n/counter
  return(list(sample = sample, acc = acc))
}


if(T) {
x = seq(-5, 15, by = 0.5)
joint = function(x) dnorm(x, mean = 5, sd = 2)
renvelope = function(x) rnorm(x, 5, 2)
fenvelope = function(x) 1.3*dnorm(x, 5, 2)

plot(x, joint(x), type = "l", col = 4, lwd = 2, ylim = c(0,0.26))
lines(x, fenvelope(x), col = 2, lty = 2) 
legend("topright", c("Target", "sampling pdf"), col = c(4,2), lty = c(1,2))

sample = RS(joint, renvelope, fenvelope, sim = T, roi = c(-5, 15), pause = 2, n = 30000)
s = sample$sample
plot(s)
hist(s, breaks = 40, probability = T, col = "yellow", ylim = c(0, 0.26))
plot(x, dnorm(x, mean = 5, sd = 2), col = 4, lwd = 2)
lines(x, fenvelope(x), col = 2, lty = 2)
#acf(s)
sample$acc
}

