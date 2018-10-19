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


