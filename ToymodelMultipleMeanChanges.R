###libraries used
library(ernest)
library(ggplot2)
library(ggdist)
library(patchwork)

###multiple mean change points

#Data input
#set.seed(123)
n <- 500 
tau1 <- 40 
tau2 <- 80
tau3 <- 280
Y  <- c(rnorm(tau1, mean=0, sd=1),rnorm(tau2-tau1,mean=2,sd=1),
        rnorm(tau3-tau2,mean=3,sd=1),rnorm(n-tau3,mean=6)) #samples

#set.seed(Sys.time())
#Loglikelihood for k change points


loglik_Mk <- function(k){
  part <- numeric(k+1)
  function(theta){
    if (k==0){
      part[1] <- sum((Y[1:n]-theta[1])^2) 
    }
    else{
      part[1] <- sum((Y[1:theta[2]]-theta[1])^2)
    }
    if (k>0){
      part[k+1] <- sum((Y[(theta[2*k]+1):n]-theta[2*k+1])^2)
    }
    if (k>1){
      for(i in 2:k){
        part[i] <- sum((Y[(theta[2*i-2]+1):theta[2*i]]-theta[2*i-1])^2)
      }
    }
    -n/2*log(2*pi)-1/2*sum(part)
  }
}

#priors for k change point model
uniftr <- function(u) {
  floor((n-1)*u) + 1
}
tau_Mk <- create_prior(vectorized_fn=uniftr, names="tau", repair="unique_quiet") #discrete uniform
#tau_Mk <- create_uniform_prior(names= "tau", lower=1,upper=n-1) #uniform
mu_Mk <- create_uniform_prior(names="mu",
                              lower= min(Y), 
                              upper =max(Y), repair="unique_quiet")
prior_Mk <- function(k){ ### for each proposed change point we get 2 more parameters
  finalprior <- mu_Mk+tau_Mk
  if (k==1){
    finalprior <- finalprior + mu_Mk
  }
  if (k>=2){
    for (i in 2:k){
      finalprior <- finalprior + mu_Mk+tau_Mk 
    }
    finalprior <- finalprior + mu_Mk
  }
  finalprior
}

#test to see what percentage we get correct
iterations <- 100
correct <- 0
logevidence_Mk <- matrix(NA, nrow = 100, ncol = 6)
for(i in 0:5){
  Sampler_Mk <-ernest_sampler(loglik_Mk(i), prior_Mk(i), nlive=500)
  for(j in 1:iterations){
  set.seed(Sys.time())
  Result_Mk <- generate(Sampler_Mk)
  logevidence_Mk[j,i+1] <- Result_Mk$log_evidence
  }
}

for(i in 1:iterations){
  max <- which.max(logevidence_Mk[i,])
  if (max==4){
    correct <- correct + 1
  }
}
correct

###
###results: number of change points 3, models (M0,...,M5)
###N=10,n=500,Delta=1, correct=18/100
###for N=20 we have 32/100 correct, 28 and 24 for n=50 (2 runs), 54 correct 
###for n=100, 75 correct for n=150, 82/100 for N=200, 89/100 for N=250