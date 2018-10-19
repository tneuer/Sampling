lower = 6
upper = 10
A = 1
B = -10
C = 8

f = function(x) {
  A*exp(B*(x-C))
}

x = seq(lower, upper, by = 0.05)
norm1 = integrate(f, lower, upper)$value
norm2 = A/B * (exp(B*(upper-C)) - exp(B*(lower-C)))

pdf = function(x) {
  1/norm1 * f(x)
}

plot(x, pdf(x))

sampler = function(r) {
  log(B*r*norm1/A + exp(B*(lower-C)))/B + C
}

val = runif(10000)
sample = sampler(val)

hist(sample, probability = T)
lines(x, pdf(x))
max(sample)
min(sample)
