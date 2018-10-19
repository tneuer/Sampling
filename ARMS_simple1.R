dev.off()
rm(list=ls())

start = Sys.time()
ARS <- function(joint, initial, samplesize = 5000, sim = F, pause = 1, bounds = c(-Inf,Inf), roi = c(-100, 100)) {
  h = function(x) log(joint(x))
  
  if (sim) {
    x = seq(roi[1], roi[2], by = 0.3)
    hx = h(x)
    plot(x, hx, type = "l")
  }
  
  T_k = sort(initial)
  h_val = h(T_k)
  l_T = length(T_k)
  
  derivate = function(f, x0) {
    #Funtion to calculate the derivative of f in x0
    h = 10e-4
    (f(x0+h) - f(x0-h))/(2*h)
  }
  h_deriv_val = derivate(h, T_k)
  
  ins_to_sorted = function(n, y, index = 0) {
    #Algorithm to find the position of a new element n in an already sorted list y
    l = length(y)
    if (l == 1) {
      if (y < n) {
        return(index+1)
      }
      else return(index)
    }
    if (n > y[l%/%2]) {
      ins_to_sorted(n, y[(l%/%2+1):l], index = index + l%/%2)
    }
    else {
      ins_to_sorted(n, y[1:(l%/%2)], index = index)
    }
  }
  
  u_k = function(x) {
    #new upper hull, former lower hull
    i = ins_to_sorted(x, T_k)
    if (i == 0 || i == length(T_k)) {
      return(-Inf)
    }
    return(((T_k[i+1] - x) * h_val[i] + (x - T_k[i]) * h_val[i+1])/(T_k[i+1] - T_k[i]))
  }
  
  exp_u = function(x) {
    return(exp(u_k(x)))
  }
  
  if (sim) {
    #Drawing initial upper and lower hull
    for (xi in x) {
      points(xi, u_k(xi))
    }
  }
  
  calc_norms = function(pos) {
    #Calculates the normalizaton constant for the bounded exponential function
    exp(h_val[pos])/h_deriv_val[pos] * (exp(h_deriv_val[pos]*(T_k[pos+1]-T_k[pos]))-exp(h_deriv_val[pos]*(T_k[pos]-T_k[pos])))
  }
  
  #Normalizing constants for piecewise exponentials, their weights and the partial sum.
  norms = calc_norms(seq(l_T-1))
  weights = norms/sum(norms)
  partials = c(0,cumsum(weights))
  
  choose_sampler = function(part) {
    #Chooses according to the weight from which piecewise exponential it should be sampled. 
    #The greater the weight the more probable is a sample from this interval
    i = ins_to_sorted(part, partials)
    return(sampler_exp(i, runif(1)))
  }
  
  sampler_exp = function(i, r) {
    #Using the inverse cdf transformation of a bounded exponential to sample
    # -i: interval index of sample
    # -r: random uniform number
    # -output: independent random sample frm bounded exponential
    y = h_deriv_val[i]*r*norms[i]/exp(h_val[i]) + exp(h_deriv_val[i]*(T_k[i]-T_k[i]))
    s = T_k[i] + log(y)/h_deriv_val[i]
    return(s)
  }
  
  MH_step = function(sample, newx) {
    u = runif(1)
    curx = sample[length(sample)]
    jcurx = joint(curx)
    jnewx = joint(newx)
    alpha = jnewx * min(jcurx, exp_u(curx) )/ (jcurx * min(jnewx, exp_u(newx)))
    if (u < alpha) sample = c(sample, newx)
    else sample = c(sample, curx)
  }
  
  sample = c(0)
  history_impr = c()
  count = improvals = tmp = feval = 0
  
  #start sampling
  while (count < samplesize) {
    print(count)
    
    x_proposal = choose_sampler(runif(1))
    if (x_proposal < bounds[1]+1e-3 || x_proposal > bounds[2]-1e-3) {
      next()
    }
    w = runif(1)
    u_val = u_k(x_proposal)

    #Rejection test
    if (w <= exp(h(x_proposal) - u_val)) {
      sample = MH_step(sample, x_proposal)
      count = count + 1
      feval = feval + 3
    }
    #Include the rejected point into the set for the uper and lower hull
    else {
      pos = ins_to_sorted(x_proposal, T_k) +1
      T_k = append(T_k, x_proposal, after = pos - 1)
      l_T = length(T_k)
      new_h_val = h(x_proposal)
      h_val = append(h_val, new_h_val, after = pos-1)
      new_h_deriv_val = derivate(h, x_proposal)
      h_deriv_val = append( h_deriv_val, new_h_deriv_val, after = pos-1)
      feval = feval + 3
      if (pos == 1) {
        new_norm = calc_norms(pos)
        norms = append(norms, new_norm, after = pos-1)
        norms[pos+1] = calc_norms(pos+1)
      }
      else if (pos == length(T_k)) {
        new_norm = calc_norms(pos)
        norms = append(norms, new_norm, after = pos-1)
        norms[pos-1] = calc_norms(pos-1)
      }
      else {
        new_norm = calc_norms(pos)
        norms = append(norms, new_norm, after = pos-1)
        norms[pos-1] = calc_norms(pos-1)
        norms[pos+1] = calc_norms(pos+1)
      }
      weights = norms/sum(norms)
      partials = c(0,cumsum(weights))
      improvals = improvals + 1
      tmp = count - tmp
      if (sim) {
        plot(x, hx, type = "l", col = 2, lwd = 3, main = paste("improved ", improvals, " times; Last one: ", tmp))
        u_vals = seq(l = length(x))
        for (i in 1:length(x)) {
          u_vals[i] = u_k(x[i])
        }
        points(x_proposal, log(w) + u_val)
        lines(x, u_vals, col = 3, lwd = 0.9)
        print(paste("improved ", improvals, " times; Last one: ", tmp))
        Sys.sleep(pause)
      }
      history_impr = c(history_impr, tmp)
      tmp = count
    }
    
  }
  return(list(sample = sample, improvals = improvals, history = history_impr, feval =  feval))
  
}

f1 <- function(x) {
  0.25*(dnorm(x, mean = 12, sd = 20) + dnorm(x, mean = -10, sd = 3) + dnorm(x, mean = 100, sd = 5) + dgamma(x-50, shape = 3, rate = 0.04))
}

GO = ARS(joint = f1, initial = c(-12, -3, 10,14, 90, 110, 70, 55), sim = T,  pause = 0.5, samplesize = 20000, bounds = c(-Inf, Inf), roi = c(-50, 150))
s = GO$sample
h = GO$history
f = GO$feval
plot(h)
abline(h = mean(h))
acf(s)
hist(s, probability = T, breaks = 100, col = "yellow")

x = seq(-100, 250, by = 0.1)
lines(x, f1(x), lty = 2, col = 4, lwd = 2)

print(paste(GO$improvals, mean(s), sd(s), f))
Sys.time() - start
