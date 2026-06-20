###libraries used
library(ernest)
library(ggplot2)
library(ggdist)
library(patchwork)

###Normal Mixture with nested sampling, be aware that labels could be switched

#data
#set.seed(123)
n=1000
tau=333
p <- 1/2
q <- 1/3
Y <- numeric(n)
for (i in 1:tau){
  if (runif(1) <= p){
    Y[i] <- rnorm(1,mean=0,sd=0.2)
  }
  else{
    Y[i] <- rnorm(1,mean=3,sd=0.3)
  }
}
for (i in (tau+1):n){
  if (runif(1) <= q){
    Y[i] <- rnorm(1,mean=6,sd=0.2)
  }
  else{
    Y[i] <- rnorm(1,mean=15,sd=1)
  }
}
hist(Y, breaks = 100, freq = FALSE)

#set.seed(Sys.time())
#The loglikelihoods of models 0 and 1

loglik_M0 <- function(theta){
  #list of parameters
  mu1 <- theta[1]
  mu2 <- theta[2]
  sd1 <- theta[3]
  sd2 <- theta[4]
  p1 <- theta[5]
  if(sd1 <= 0 || sd2 <= 0 || p1 <= 0 || p1 >= 1) return(-Inf)
  
  #stabilize exponentials
  log_mix <- function(logp, log1mp, logf1, logf2) {
    a <- logp + logf1
    b <- log1mp + logf2
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  
  #loglikelihood
  logdens1_1 <- dnorm(Y[1:n], mu1, sd1, log=TRUE)
  logdens1_2 <- dnorm(Y[1:n], mu2, sd2, log=TRUE)
  part1 <- sum(log_mix(log(p1), log(1-p1), logdens1_1, logdens1_2))
  part1
}

loglik_M1 <- function(theta){
  #list of parameters
  mu1 <- theta[1]
  mu2 <- theta[2]
  mu3 <- theta[3]
  mu4 <- theta[4]
  sd1 <- theta[5]
  sd2 <- theta[6]
  sd3 <- theta[7]
  sd4 <- theta[8]
  p1 <- theta[9]
  p2 <- theta[10]
  tau <- theta[11] #replace with as.integer(theta[11]) when using uniform prior
  
  if(sd1 <= 0 || sd2 <= 0 || sd3 <= 0 || sd4 <= 0) return(-Inf)
  if(p1 <= 0 || p1 >= 1 || p2 <= 0 || p2 >= 1) return(-Inf)
  
  #stabilize exponentials
  log_mix <- function(logp, log1mp, logf1, logf2) {
    a <- logp + logf1
    b <- log1mp + logf2
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  
  #loglikelihood
  logdens1_1 <- dnorm(Y[1:tau], mu1, sd1, log=TRUE)
  logdens1_2 <- dnorm(Y[1:tau], mu2, sd2, log=TRUE)
  part1 <- sum(log_mix(log(p1), log(1-p1), logdens1_1, logdens1_2))
  
  logdens2_1 <- dnorm(Y[(tau+1):n], mu3, sd3, log=TRUE)
  logdens2_2 <- dnorm(Y[(tau+1):n], mu4, sd4, log=TRUE)
  part2 <- sum(log_mix(log(p2), log(1-p2), logdens2_1, logdens2_2))
  
  part1 + part2
}

#priors for model 0 and 1
sdm <- sqrt(max((mean(Y)-Y)^2))

uniftr <- function(u) {
  floor((n-1)*u) + 1
}
tau_M1 <- create_prior(vectorized_fn=uniftr, names="tau")#discrete uniform
#tau_M1 <- create_uniform_prior(names= "tau", lower=1,upper=n-1) #uniform

prior_M0 <- create_uniform_prior(names=c("mu1","mu2","sd1","sd2","p1"), 
                                 lower=c(min(Y),min(Y),0,0,0), 
                                 upper=c(max(Y),max(Y),sdm,sdm,1))
prior_M1 <- create_uniform_prior(names=c("mu1","mu2","mu3","mu4",
                                         "sd1","sd2","sd3","sd4",
                                         "p1","p2"), 
                                 lower=c(min(Y),min(Y),min(Y),min(Y),
                                         0,0,0,0,0,0), 
                                 upper=c(max(Y),max(Y),max(Y),max(Y),
                                         sdm,sdm,sdm,sdm,1,1))+tau_M1

#samplers
sampler_M0 <- ernest_sampler(loglik_M0, prior_M0, nlive=1000,
                             sampler= multi_ellipsoid(enlarge = 1.5))
sampler_M1 <- ernest_sampler(loglik_M1, prior_M1, nlive=1000)

#Nested sampling run

result_M0 <- generate(sampler_M0)

result_M1 <- generate(sampler_M1)

#Summary of results
summary(result_M0)
summary(result_M1)

logBF_10 <- result_M1$log_evidence-result_M0$log_evidence
logBF_10

#posterior distributions

visualize(result_M0,mu1,.which="density")
visualize(result_M0,mu2,.which="density")
visualize(result_M0,sd1,.which="density")
visualize(result_M0,sd2,.which="density")
visualize(result_M0,p1,.which="density")
visualize(result_M1,mu1,.which="density")
visualize(result_M1,mu2,.which="density")
visualize(result_M1,mu3,.which="density")
visualize(result_M1,mu4,.which="density")
visualize(result_M1,sd1,.which="density")
visualize(result_M1,sd2,.which="density")
visualize(result_M1,sd3,.which="density")
visualize(result_M1,sd4,.which="density")
visualize(result_M1,p1,.which="density")
visualize(result_M1,p2,.which="density")
visualize(result_M1,tau,.which="density")