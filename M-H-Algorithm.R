dev.off()
rm(list = ls())

start = Sys.time()
MH = function(joint, initial, samplesize = 10000, stepsize = seq(1,1,l=length(initial)), burn_in = 4000, thin = 1, sim = F, pause = 0.5, jumper = c(samplesize,1), roi = c(-20,20)) {
  # An Metropolis Hastings sampling algorithm. Its main advantage is the plot function to see what is happening.
  # Input:  joint -       joint distribution to be samped
  #         initial -     intial value for the Metropolis Algorithm. Also indicates the dimensionality of the pdf
  #         samplesize -  Desired samplesize, has to be increased if thin != 1.
  #         stepsize -    Used standard deviation for the Gaussian candidate generating function
  #         burn_in -     Number of samples you want to discard in the beginning until the algorithm hopefully converges
  #                       Does not affect samplesize
  #         thin -        Thinning parameter, integer value
  #         sim -         True if the sample process should be shown
  #         pause -       Only relevant if sim = T; How long in seconds beofre a new frame is shown
  #         jumper -      Only relevant if sim = T; Tuple of integers, first number indicates how many samples should be drawn,
  #                       second value declares how many values should be left out before next values are simulated.
  #                       Shows different stages of the sampling
  #         roi -         Only relevant if sim = T; In which region the joint density is relevant
  #         
  # Output: sample - vector of samples with length samplesize/thin
  #         acc - acceptance rate of the algorithm
  
  nr_param = length(initial) #dimensionality of the problem
  samples = matrix(0, nrow = samplesize/thin, ncol = nr_param) #rows are different variables, columns are the samples
  tmp = initial #saves last accepted sample for the next step
  user_pause = pause #Saved for simulation later
  changer = F #Variable to make jumper cahnge possible
  
  #candidate generating function is by default multivariate normal
  proposal = function(theta) {
    rnorm(n = nr_param, mean = theta, sd = stepsize)
  }
  
  #candidate generating function, independence chain
  # proposal = function(theta) {
  #   runif(1, -20,250)
  # }
  symmetric = T #Saves calculation time if symmetric candidate genrating function is used (like Gaussian)
  
  x = seq(roi[1],roi[2],length.out = 100) #needed for plotting the candidate generating function
  counter = counter2 = counter3 = 0 #counter: to calculate acceptance rate; counter2 and counter 3 handle the jumper possibility
  ask1 = T #Variables to correspond to users opinion if the jumper should go on
  ask2 = F
  plotlast = 1000 #to save computational time, for the traceplot only the last samples are shown
  par(mfrow=c(1,2))
  
  for (i in -burn_in:samplesize) {
    proposal_val = proposal(tmp) #generate proposal
    
    #calculation of step probability
    if (symmetric) {
      acc = joint(proposal_val)/joint(tmp)
    }
    else {
      acc = joint(proposal_val)/joint(tmp) * dnorm(tmp, mean = proposal_val, sd = stepsize)/dnorm(proposal_val, mean = tmp, sd = stepsize)
    }
    
    #Draw new random uniform to compare with acc
    randunif = runif(1)
    new = randunif < acc
    if (new) {
      tmp2 = tmp
      tmp = proposal_val #proposal gets promoted to new position for the next draw
      if (i > 0 && (i %% thin == 0)) { #Burn_in is over and thinning condition fullfilled
        counter = counter + 1
        samples[i / thin, ] = proposal_val
      }
    }
    else if (i > 0 && (i %% thin == 0)) {
      tmp2 = tmp
      samples[i / thin, ] = tmp #old value is taken again
    }
    
    #The actual algorithm is done, now comes the plotting part.
    
    if (sim && nr_param == 1 && i > 0 && i %% thin == 0) { #Only one dimensional distributions are plotted, no burn_in
      if ((counter2 == 0) && (counter3 == 0)) {
        if (ask1 && ask2) {
          answer = readline(prompt = "Speed up? [Y / N / Never / Finish / Change] ")
          if (answer == "Y" || answer == "y") {
            pause = 0.005
          }
          else if (answer == "N" || answer == "n") {
            pause = user_pause
          }
          else if (answer == "Never" || answer == "never") {
            ask1 = F
          }
          else if (answer == "Finish" || answer == "finish" || answer == "f") {
            sim = F
          }
          else if (answer == "Change" || answer == "changer" || answer == "c") {
            changer = T
          }
        }
        
        if (changer) {
          jumper[1] = as.integer(readline(prompt = "Steps drawn: "))
          jumper[2] = as.integer(readline(prompt = "Steps left out: "))
        }
        counter2 = jumper[2]
        counter3 = jumper[1]
        ask2 = T
        changer = F
      }
      if(counter3 != 0) {
        counter3 = counter3 - 1
      }
      else if (counter2 != 0) {
        counter2 = counter2 - 1
        next()
      }
      
      
      if (new) color = "green" #color for accepted value is green
      else color = 2 #color for denied value is red
      plot(1, 1, xlim = c(1,samplesize), ylim = c(x[1], x[length(x)]), type = "n", main = paste("n = ",i, "acc = ", round(counter/i, digits = 3)))
      lines(-dnorm(x, mean = tmp2, sd = stepsize)*samplesize-i+samplesize, x, type = "l") #density function
      lines(-dnorm(x, mean = 0, sd = 1)*samplesize-i+samplesize, x, lty = 2, col = 6) #density function
      lines(c(0,samplesize), c(proposal_val, proposal_val), col = color) #line through the current proposal
      lines(c(0,samplesize), c(tmp2, tmp2), col = 1) #line through the last accepted point
      lines(c(0, samplesize), c(0,0))
      #Traceplot
      if (i < plotlast+1) {
        points(seq(samplesize, samplesize-i+1), samples[1:i, ]) 
      }
      else {
        points(seq(samplesize, samplesize-i+1)[(i-plotlast):i], samples[(i-plotlast):i, ])
      }
      #Histogram of the samples
      barplot(hist(samples[1:(i/thin), ], plot = F, breaks = 20, prob = T)$counts, horiz = T, col = "yellow", main = paste("Prob = ", round(acc, 4), "Uniform = ", round(randunif, 4)))
      
      Sys.sleep(pause)
    }
  }
  dev.off()
  return(list(sample = samples, acc = counter/samplesize))
}

if (T) {
  
  f = function(x) {
    0.25*(dnorm(x, mean = 12, sd = 20) + dnorm(x, mean = -10, sd = 2) + dnorm(x, mean = 200, sd = 5) + dgamma(x-150, shape = 3, rate = 0.2))
  }
  roi = c(-20, 250)
  
  f = function(x) {
    dnorm(x, mean = 0, sd = 1)
  }
  roi = c(-5,5)
  
  sample = MH(joint = f, initial = 1, samplesize = 30000, stepsize = 1, sim = F, pause = 5, burn_in = 0, jumper = c(10,0), thin = 17, roi = roi)
  sample$acc
  s = sample$sample
  plot(s, type = "l", main = paste("Acceptance rate: ",sample$acc))
  hist(s, probability = T, col = "yellow", breaks = 130, main = "")
  roi = seq(roi[1], roi[2], 0.01)
  lines(roi, f(roi), col = 2, lwd = 2, lty = 2)
  acf(s, main = "Thin = 10")
  mean(s)
  # (12 - 10 + 200 + 150+3/0.2) /4
  
  Sys.time() - start
}

