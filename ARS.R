# try(dev.off(), silent = T)
# rm(list=ls())

start = Sys.time()
ARS <- function(joint, initial, samplesize = 40000, sim = F, pause = 0.05, bounds = c(-Inf,Inf), roi = c(-100, 100)) {
  h = function(x) log(joint(x))
  
  if (sim) {
    x = seq(roi[1], roi[2], by = 0.2)
    hx = h(x)
    plot(x, hx, type = "l")
  }
  
  T_k = sort(initial)
  h_val = h(T_k)
  l_T = length(T_k)
  
  derivate = function(f, x0) {
    d <- rep(0, length(x0))
    #Funtion to calculate the derivative of f in x0
    h = 10e-4
    for (i in 1:length(x0)) {
    d[i] <- (f(x0[i]+h) - f(x0[i]-h))/(2*h)
    }
    return(d)
  }
  h_deriv_val = derivate(h, T_k)
  
  if (length(which(h_deriv_val<0)) == 0 || length(which(h_deriv_val>0)) == 0 ) {
    #Check if the initial values fulfill certain conditions
    print("One value to the left as well as on the rigth side of the mode is needed. Now: ") 
    print(h_deriv_val)
    print("At least one should be negative and one positive!")
    stop()
  }
  
  intersection = function(pos) {
    #calculate intersections of the tangents
    (h_val[pos+1] - h_val[pos] - T_k[pos+1]*h_deriv_val[pos+1] + T_k[pos]*h_deriv_val[pos])/(h_deriv_val[pos] - h_deriv_val[pos+1])
  }
  
  ins_to_sorted = function(n, y, index = 0) {
    #Algorithm to find the position of a new element in an already sorted list
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
    #upper hull
    i = ins_to_sorted(x, intersections)
    return(h_val[i] + (x - T_k[i]) * h_deriv_val[i])
  }
  
  l_k = function(x) {
    #lower hull
    i = ins_to_sorted(x, T_k)
    if (i == 0 || i == length(T_k)) {
      return(-Inf)
    }
    return(((T_k[i+1] - x) * h_val[i] + (x - T_k[i]) * h_val[i+1])/(T_k[i+1] - T_k[i]))
  }
  
  exp_u = function(x) {
    return(exp(u_k(x)))
  }
  
  #Initial intersections
  intersections = intersection(1:(length(T_k)-1))
  intersections = c(-Inf, intersections, Inf)
  
  if (sim) {
    #Drawing initial upper and lower hull
    for (xi in x) {
      points(xi, u_k(xi))
      points(xi, l_k(xi))
    }
  }
  
  calc_norms = function(pos) {
    #Calculates the normalizaton constant for the bounded exponential function
    exp(h_val[pos])/h_deriv_val[pos] * (exp(h_deriv_val[pos]*(intersections[pos+1]-T_k[pos]))-exp(h_deriv_val[pos]*(intersections[pos]-T_k[pos])))
  }
  
  #Normalizing constants for piecewise exponentials, their weights and the partial sum.
  norms = calc_norms(seq(l_T))
  weights = norms/sum(norms)
  partials = c(0,cumsum(weights))
  
  choose_sampler = function(part) {
    #Chooses according to the weight from which piecewise exponential it should be sampled. 
    #The greater the weight the more probable is a samle from this piece
    i = ins_to_sorted(part, partials)
    return(sampler_exp(i, runif(1)))
  }
  
  sampler_exp = function(i, r) {
    #Using the inverse cdf transformation of a bounded exponential to sample
    y = h_deriv_val[i]*r*norms[i]/exp(h_val[i]) + exp(h_deriv_val[i]*(intersections[i]-T_k[i]))
    s = T_k[i] + log(y)/h_deriv_val[i]
    return(s)
  }
  
  sample = c()
  history_impr = c()
  count = improvals = tmp = feval = 0
  
  #start sampling
  while (length(sample) < samplesize) {
    print(count)
  
    x_proposal = choose_sampler(runif(1))
    if (x_proposal < bounds[1]+1e-3 || x_proposal > bounds[2]-1e-3) {
      next()
    }
    w = runif(1)
    u_val = u_k(x_proposal)
    
    #Squeezing test
    if (w <= exp(l_k(x_proposal) - u_val)) {
      sample = c(sample, x_proposal)
      count = count + 1
    }
    #Rejection test
    else if (w <= exp(h(x_proposal) - u_val)) {
      sample = c(sample, x_proposal)
      count = count + 1
      feval = feval + 1
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
      if (pos == 1) {
        new_intersect = intersection(pos)
        intersections = append(intersections, new_intersect, after = pos)
        new_norm = calc_norms(pos)
        norms = append(norms, new_norm, after = pos-1)
        norms[pos+1] = calc_norms(pos+1)
      }
      else if (pos == length(T_k)) {
        new_intersect = intersection(pos-1)
        intersections = append(intersections, new_intersect, after = pos-1)
        new_norm = calc_norms(pos)
        norms = append(norms, new_norm, after = pos-1)
        norms[pos-1] = calc_norms(pos-1)
      }
      else {
        new_intersect1 = intersection(pos-1)
        new_intersect2 = intersection(pos)
        intersections[pos] = new_intersect1
        intersections = append(intersections, new_intersect2, after = pos)
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
        plot(x, hx, type = "l", col = 2, lwd = 3, main = paste("improved ", improvals, " times; Last one: ", tmp, "; sample: ", length(sample), "/", samplesize))
        u_vals = seq(l = length(x))
        l_vals = seq(length.out = length(x))
        for (i in 1:length(x)) {
          u_vals[i] = u_k(x[i])
          l_vals[i] = l_k(x[i])
        }
        points(x_proposal, log(w) + u_val)
        lines(x, u_vals, col = 3, lwd = 0.9)
        lines(x, l_vals, col = 4, lwd = 0.9)
        #for (xi in x) {
          #points(xi, u_k(xi), col = 3, lwd = 0.3)
          #points(xi, l_k(xi), col = 4, lwd = 0.3)
        #}
        Sys.sleep(pause)
      }
      history_impr = c(history_impr, tmp)
      tmp = count
    }
    
  }
  return(list(sample = sample, improvals = improvals, history = history_impr, feval =  feval))
  
}

if(F) {
f1 <- function(x) {
  dweibull(x, shape = 10, scale = 20)
}

GO = ARS(joint = f1, initial = c(2, 30), sim = T,  pause = 0.5, samplesize = 40000, bounds = c(0, Inf), roi = c(0, 30))
s = GO$sample
h = GO$history
f = GO$feval
plot(h)
abline(h = mean(h))
acf(s)
hist(s, probability = T, breaks = 100, col = "yellow")

x = seq(-60, 70, by = 0.1)
lines(x, f1(x), lty = 2, col = 2, lwd = 2)

print(paste(GO$improvals, mean(s), sd(s), f))
Sys.time() - start
}

if(F) {
library("ars", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
start = Sys.time()

#Example 3: sample 20 values from a beta(1.3,2.7) distribution
f2<-function(x,a,b){(a-1)*log(x)+(b-1)*log(1-x)}
f2prima<-function(x,a,b){(a-1)/x-(b-1)/(1-x)}
mysample2<-ars(30000,f2,f2prima,x=c(0.3,0.6),m=2,lb=TRUE,xlb=0,ub=TRUE,xub=1,a=1.3,b=2.7)
hist(mysample2, probability = T, breaks = 20)
x = seq(0,1,by=0.01)
lines(x, dbeta(x, 1.3, 2.7))

Sys.time() - start
}
