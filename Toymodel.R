###
#Made by: Thomas Ammerlaan
#Institution: Leiden University

###libraries used
library(ernest)
library(ggplot2)
library(ggdist)
library(patchwork)

###mean change model
#samples
#set.seed(123)
n <- 100
tau <- 50
Y <- c(rnorm(tau),rnorm(n-tau, mean= 1))

hist(Y, breaks = 25, freq = FALSE)

#set.seed(Sys.time())
#loglikelihoods of model 0 and 1


loglik_M0 <- function(theta){
  mu <- theta[1]
  sum(dnorm(Y[1:n], mean = mu, sd = 1, log = TRUE))
}
loglik_M1 <- function(theta){
  mu1 <-theta[1]
  mu2 <-theta[2]
  tau <-theta[3] #replace with as.integer(theta[3]) when using uniform prior
  part1 <- sum(dnorm(Y[1:tau], mean = mu1, sd = 1, log = TRUE))
  part2 <- sum(dnorm(Y[(tau+1):n], mean = mu2, sd = 1, log = TRUE))
  part1+part2
}

#priors of model 0 and 1
prior_M0 <- create_normal_prior(names=c("mu1","dummy"))

norm_M1 <- create_normal_prior(names=c("mu1","mu2"))
#unif_M1 <- create_uniform_prior(names="tau", lower=1, upper=n-1)

uniftr <- function(u) {
  floor((n-1)*u) + 1
}
unif_M1 <- create_prior(vectorized_fn=uniftr, names="tau")#discrete uniform
#tau_M1 <- create_uniform_prior(names= "tau", lower=1,upper=n-1) #uniform
prior_M1 <- norm_M1+unif_M1

#theoretical values of evidences
Z_0 <- 1/(sqrt(2*pi)^n*sqrt(n+1))*exp(-1/2*(sum((Y-mean(Y))^2)+n*mean(Y)^2/(n+1)))
Z_1 <- 0
for (i in 1:(n-1)){
  part1 <- 1/(sqrt(i+1)*sqrt((n-i)+1))
  part2 <- exp(-1/2*(sum((Y[1:i]-mean(Y[1:i]))^2)+i*mean(Y[1:i])^2/(i+1)))
  part3 <- exp(-1/2*(sum((Y[(i+1):n]-mean(Y[(i+1):n]))^2)+(n-i)*mean(Y[(i+1):n])^2/((n-i)+1)))
  Z_1 <- Z_1+part1*part2*part3
}
Z_1 <- Z_1/((n-1)*sqrt(2*pi)^n)
BFtheoretical_10 <- Z_1/Z_0


sampler_M0 <- ernest_sampler(loglik_M0, prior_M0, nlive=1000)
sampler_M1 <- ernest_sampler(loglik_M1, prior_M1, nlive=1000)
result_M0 <- generate(sampler_M0)
result_M1 <- generate(sampler_M1)

#results

summary(result_M0)
summary(result_M1)

#Bayesfactor BF_10 (both theoretical and numerical)
BF_10 <- exp(result_M1$log_evidence-result_M0$log_evidence)
BF_10

BFtheoretical_10

#posterior distributions
visualize(result_M0,mu1, .which = "density")
visualize(result_M1,mu1, .which = "density")
visualize(result_M1,mu2, .which = "density")
visualize(result_M1,tau, .which = "density")

#nested sampling error distribution
Zhat0 <- numeric()
Zhat1 <- numeric()

for(i in 1:5000){
  set.seed(Sys.time())
  sampler_M0 <- ernest_sampler(loglik_M0, prior_M0, nlive=100)
  sampler_M1 <- ernest_sampler(loglik_M1, prior_M1, nlive=100)
  result_M0 <- generate(sampler_M0)
  result_M1 <- generate(sampler_M1)
  Zhat0[i] <- result_M0$log_evidence
  Zhat1[i] <- result_M1$log_evidence
}

H0 <- result_M0$information

H1 <- result_M1$information

epsilon0 <- seq(min(Zhat0-log(Z_0)),max(Zhat0-log(Z_0)),by=0.001)

epsilon1 <- seq(min(Zhat1-log(Z_1)),max(Zhat1-log(Z_1)),by=0.001)

hist(Zhat0[1:100]-log(Z_0),freq=FALSE, 
     main='log(Zhat0)-log(Z0) with gaussian overlay')

points(epsilon0 ,2*exp(-100*epsilon0^2/(2*(H0+epsilon0)))
       ,type='l', col='red')

hist(Zhat1[1:200]-log(Z_1),freq=FALSE, 
     main = 'log(Zhat1)-log(Z1) with gaussian overlay')

points(epsilon1 ,2*exp(-100*epsilon1^2/(2*(H1+epsilon1)))
       ,type='l', col='red')

statsZhat <- matrix(c(Zhat0,Zhat1,Zhat0-log(Z_0),Zhat1-log(Z_1))
                    ,nrow=5000,ncol=4)

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

###mean and variance change
#data input
#set.seed(123)
n <- 100
tau <- 50
Y  <- c(rnorm(tau, mean=0, sd=1),rnorm(n-tau,mean=2,sd=2))
hist(Y, breaks = 25, freq = FALSE)

#set.seed(Sys.time())
#loglikelihoods for model 0 and 1

loglik_M0 <- function(theta){
  mu <- theta[1]
  sd <- theta[2]
  -n/2*(log(2*pi)+2*log(sd)) - 0.5*sum(((Y-mu)/theta[2])^2)
}

loglik_M1 <- function(theta){ 
  mu1 <- theta[1]
  mu2 <- theta[2]
  sd1 <- theta[3]
  sd2 <- theta[4]
  tau <- theta[5] #replace with as.integer(theta[5]) when using uniform prior
  part1 <- sum(((Y[1:tau]-mu1)/sd1)^2)
  part2 <- sum(((Y[(tau+1):n] - (mu2))/sd2)^2)
  -n/2*log(2*pi) - tau*log(sd1)-(n-tau)*log(sd2)-0.5*(part1 + part2)
}

#priors for model 0 and 1
sdm <- sqrt(max((mean(Y)-Y)^2))

uniftr <- function(u) {
  floor((n-1)*u) + 1
}
tau_M1 <- create_prior(vectorized_fn=uniftr, names="tau") #discrete uniform
#tau_M1 <- create_uniform_prior(names= "tau", lower=1,upper=n-1) #uniform
Prior_M0 <- create_uniform_prior(names=c("mu", "sd"), 
                                 lower= c(min(Y),0), 
                                 upper =c(max(Y),sdm))
Prior_M1 <- create_uniform_prior(names=c("mu1", "mu2","sd1","sd2"), 
                                 lower= c(min(Y),min(Y),0,0), 
                                 upper =c(max(Y),max(Y),sdm,sdm))+tau_M1

#nested sampling
sampler_M0 <- ernest_sampler(loglik_M0, Prior_M0, nlive=1000)
sampler_M1 <- ernest_sampler(loglik_M1, Prior_M1, nlive=1000)

#results
result_M0 <- generate(sampler_M0)
result_M1 <- generate(sampler_M1)

summary(result_M0)
summary(result_M1)

logBF_10 <- result_M1$log_evidence-result_M0$log_evidence
logBF_10

#posterior distributions

visualize(result_M0,mu, .which = "density")
visualize(result_M0,sd, .which = "density")
visualize(result_M1,mu1, .which = "density")
visualize(result_M1,mu2, .which = "density")
visualize(result_M1,sd1, .which = "density")
visualize(result_M1,sd2, .which = "density")
visualize(result_M1,tau, .which = "density")

sim_M0 <- calculate(Result_M0, ndraws = 1000)
plot(sim_M0, which = c("weight", "likelihood","evidence"))
sim_M1 <- calculate(Result_M1, ndraws = 1000)
plot(sim_M1, which = c("weight", "likelihood","evidence"))

visualize(result_M0,.which="trace")
visualize(result_M1,.which = "trace")
visualize(result_M0, .which = "density")
visualize(result_M1,.which = "density")

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

###Normal Mixture model with MCMC 
#set.seed(123)
n=100
tau=0
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
#set.seed(Sys.time())

#burn in period
burnin <- 10000

#Used for proposed samples
logit <- function(p) log(p / (1 - p))
invlogit <- function(z) 1 / (1 + exp(-z))

#We can use the loglikelihood functions from above

loglik_M0 <- function(theta){
  mu1 <- theta$mu1
  mu2 <- theta$mu2
  s1  <- theta$s1
  s2  <- theta$s2
  p1  <- theta$p1
  if(s1 <= 0 || s2 <= 0 || p1 <= 0 || p1 >= 1) return(-Inf)
  
  log_mix <- function(logp, log1mp, logf1, logf2) {
    a <- logp + logf1
    b <- log1mp + logf2
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  
  logdens1_1 <- dnorm(Y[1:n], mu1, s1, log=TRUE)
  logdens1_2 <- dnorm(Y[1:n], mu2, s2, log=TRUE)
  part1 <- sum(log_mix(log(p1), log(1-p1), logdens1_1, logdens1_2))
  part1
}

loglik_M1 <- function(theta){
  mu1 <- theta$mu1
  mu2 <- theta$mu2
  s1  <- theta$s1
  s2  <- theta$s2
  p1  <- theta$p1
  
  mu3 <- theta$mu3
  mu4 <- theta$mu4
  s3  <- theta$s3
  s4  <- theta$s4
  p2  <- theta$p2
  
  tau <- round(theta$tau)
  tau <- max(1, min(n - 1, tau))
  
  if(s1 <= 0 || s2 <= 0 || s3 <= 0 || s4 <= 0) return(-Inf)
  if(p1 <= 0 || p1 >= 1 || p2 <= 0 || p2 >= 1) return(-Inf)
  
  log_mix <- function(logp, log1mp, logf1, logf2) {
    a <- logp + logf1
    b <- log1mp + logf2
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  
  logdens1_1 <- dnorm(Y[1:tau], mu1, s1, log=TRUE)
  logdens1_2 <- dnorm(Y[1:tau], mu2, s2, log=TRUE)
  part1 <- sum(log_mix(log(p1), log(1-p1), logdens1_1, logdens1_2))
  
  logdens2_1 <- dnorm(Y[(tau+1):n], mu3, s3, log=TRUE)
  logdens2_2 <- dnorm(Y[(tau+1):n], mu4, s4, log=TRUE)
  part2 <- sum(log_mix(log(p2), log(1-p2), logdens2_1, logdens2_2))
  
  part1 + part2
}
#priors for model 0 and 1
logprior_M0 <- function(theta, x) {
  
  xmin <- min(x)
  xmax <- max(x)
  n    <- length(x)
  sdmax <- sqrt(max((mean(x)-x)^2))
  params <- c(theta$mu1, theta$mu2)
  
  # means
  if(any(params < xmin | params > xmax))
    return(-Inf)
  
  # standard deviations
  sds <- c(theta$s1, theta$s2)
  
  if(any(sds <= 0 | sds > sdmax))
    return(-Inf)
  
  # mixture weights
  ps <- c(theta$p1)
  
  if(any(ps <= 0 | ps >= 1))
    return(-Inf)
  
  # uniform prior
  return(0)
}
logprior_M1 <- function(theta, x) {
  
  xmin <- min(x)
  xmax <- max(x)
  n    <- length(x)
  sdmax <- sqrt(max((mean(x)-x)^2))
  
  params <- c(theta$mu1, theta$mu2,
              theta$mu3, theta$mu4)
  
  # means
  if(any(params < xmin | params > xmax))
    return(-Inf)
  
  # standard deviations
  sds <- c(theta$s1, theta$s2,
           theta$s3, theta$s4)
  
  if(any(sds <= 0 | sds > sdmax))
    return(-Inf)
  
  # mixture weights
  ps <- c(theta$p1, theta$p2)
  
  if(any(ps <= 0 | ps >= 1))
    return(-Inf)
  
  # changepoint
  if(theta$tau <= 1 || theta$tau >= (n - 1))
    return(-Inf)
  
  # uniform prior
  return(0)
}

# Posteriors of model 0 and 1
logpost_M0 <- function(x, theta) {
  loglik_M0(theta) + logprior_M0(theta, x)
}

logpost_M1 <- function(x, theta) {
  loglik_M1(theta) + logprior_M1(theta, x)
}

# for the BIC

llmax_M0 <- -10e300
llmax_M1 <- -10e300

#Metropolis Hasting sampling

mh_sampler_M0 <- function(x,
                          n_iter = 10000,
                          step_mu = 0.2,
                          step_sigma = 0.05,
                          step_p = 0.1) {
  
  n <- length(x)
  
  #initial values, based on available statistics
  
  theta <- list(
    mu1 = mean(x) - 1,
    mu2 = mean(x) + 1,
    s1  = sd(x),
    s2  = sd(x),
    p1  = 0.5
  )
  
  current_lp_M0 <- logpost_M0(x, theta)
  
  # Storage
  
  out <- matrix(NA, nrow = n_iter, ncol = 5)
  
  colnames(out) <- c(
    "mu1","mu2","s1","s2","p1"
  )
  
  accept <- 0
  
  # MCMC loop
  
  for(iter in 1:n_iter) {
    
    prop <- theta
    
    # Means
    
    prop$mu1 <- rnorm(1, theta$mu1, step_mu)
    prop$mu2 <- rnorm(1, theta$mu2, step_mu)
    
    # SDs (log scale, to make sure the proposal is positive)
    
    prop$s1 <- exp(rnorm(1, log(theta$s1), step_sigma))
    prop$s2 <- exp(rnorm(1, log(theta$s2), step_sigma))
    
    # Mixture weights (logit scale, so that the proposal is between 0 and 1)
    
    z1 <- rnorm(1, logit(theta$p1), step_p)
    
    prop$p1 <- invlogit(z1)
    
    # Acceptance step
    
    prop_lp_M0 <- logpost_M0(x, prop)
    
    #we keep the maximum likelihood
    if ((loglik_M0(prop) > llmax_M0) && (iter>burnin)){
      llmax_M0 <- loglik_M0(prop)
    }
    
    log_alpha <- prop_lp_M0 - current_lp_M0
    
    if(log(runif(1)) < log_alpha) {
      
      theta <- prop
      current_lp_M0 <- prop_lp_M0
      accept <- accept + 1
    }
    
    # Store sample
    
    out[iter, ] <- c(
      theta$mu1,
      theta$mu2,
      theta$s1,
      theta$s2,
      theta$p1
    )
  }
  
  cat("Acceptance rate =", accept / n_iter, "\n")
  
  return(c(as.data.frame(out),llmax_M0))
}

mh_sampler_M1 <- function(x,
                          n_iter = 10000,
                          step_mu = 0.2,
                          step_sigma = 0.05,
                          step_p = 0.1,
                          step_tau = 2) {
  
  n <- length(x)
  
  # Initial values
  
  theta <- list(
    mu1 = mean(x) - 1,
    mu2 = mean(x) + 1,
    s1  = sd(x),
    s2  = sd(x),
    p1  = 0.5,
    
    mu3 = mean(x) - 1,
    mu4 = mean(x) + 1,
    s3  = sd(x),
    s4  = sd(x),
    p2  = 0.5,
    
    tau = n / 2
  )
  
  current_lp_M1 <- logpost_M1(x, theta)
    
  # Storage
  
  out <- matrix(NA, nrow = n_iter, ncol = 11)
  
  colnames(out) <- c(
    "mu1","mu2","s1","s2","p1",
    "mu3","mu4","s3","s4","p2","tau"
  )
  
  accept <- 0
  
  # MCMC loop
  
  for(iter in 1:n_iter) {
    
    prop <- theta
    
    # Means
    
    prop$mu1 <- rnorm(1, theta$mu1, step_mu)
    prop$mu2 <- rnorm(1, theta$mu2, step_mu)
    
    prop$mu3 <- rnorm(1, theta$mu3, step_mu)
    prop$mu4 <- rnorm(1, theta$mu4, step_mu)
    
    # SDs (log scale, to make sure the proposal is positive)
    
    prop$s1 <- exp(rnorm(1, log(theta$s1), step_sigma))
    prop$s2 <- exp(rnorm(1, log(theta$s2), step_sigma))
    
    prop$s3 <- exp(rnorm(1, log(theta$s3), step_sigma))
    prop$s4 <- exp(rnorm(1, log(theta$s4), step_sigma))
    
    # Mixture weights (logit scale, so that the proposal is between 0 and 1)
    
    z1 <- rnorm(1, logit(theta$p1), step_p)
    z2 <- rnorm(1, logit(theta$p2), step_p)
    
    prop$p1 <- invlogit(z1)
    prop$p2 <- invlogit(z2)
    
    # Changepoint
    
    prop$tau <- floor(rnorm(1, theta$tau, step_tau))-1
    
    # Acceptance step
    
    prop_lp_M1 <- logpost_M1(x, prop)
    
    #we keep the maximum likelihood
    if ((loglik_M1(prop) > llmax_M1) && (iter>burnin)){
      llmax_M1 <- loglik_M1(prop)
    }
    
    log_alpha <- prop_lp_M1 - current_lp_M1
    
    if(log(runif(1)) < log_alpha) {
      
      theta <- prop
      current_lp_M1 <- prop_lp_M1
      accept <- accept + 1
    }
    
    # Store sample
    
    out[iter, ] <- c(
      theta$mu1,
      theta$mu2,
      theta$s1,
      theta$s2,
      theta$p1,
      theta$mu3,
      theta$mu4,
      theta$s3,
      theta$s4,
      theta$p2,
      theta$tau
    )
  }
  
  cat("Acceptance rate =", accept / n_iter, "\n")
  
  return(c(as.data.frame(out),llmax_M1))
}

#results iterations are chosen to compare it to the time nested sampling uses.
#~5 minutes for M0 and ~25 minutes for M1

samples_M0 <- mh_sampler_M0(Y, n_iter = 1000000)
llmax_M0 <- samples_M0[[6]]
samples_M1 <- mh_sampler_M1(Y, n_iter = 5000000)
llmax_M1 <- samples_M1[[12]]

samples_M0 <- as.data.frame(samples_M0[1:5])
samples_M1 <- as.data.frame(samples_M1[1:11])
burninsamples_M0 <- samples_M0[(burnin+1):1000000,]
burninsamples_M1 <- samples_M1[(burnin+1):5000000,]
bic0 <- -2 * llmax_M0 + 5 * log(n)
bic1 <- -2 * llmax_M1 + 11 * log(n)

logBF10_MCMC <- (bic0 - bic1)/2
logBF10_MCMC

#posteriors
hist(burninsamples_M0$mu1, probability =TRUE, main="Posterior samples of mu1", xlab="value")
hist(burninsamples_M0$mu2, probability =TRUE, main="Posterior samples of mu2", xlab="value")
hist(burninsamples_M0$s1, probability =TRUE, main="Posterior samples of s1", xlab="value")
hist(burninsamples_M0$s2, probability =TRUE, main="Posterior samples of s2", xlab="value")
hist(burninsamples_M0$p1, probability =TRUE, main="Posterior samples of p1", xlab="value")

hist(burninsamples_M1$mu1, probability =TRUE, main="Posterior samples of mu1", xlab="value")
hist(burninsamples_M1$mu2, probability =TRUE, main="Posterior samples of mu2", xlab="value")
hist(burninsamples_M1$mu3, probability =TRUE, main="Posterior samples of mu3", xlab="value")
hist(burninsamples_M1$mu4, probability =TRUE, main="Posterior samples of mu4", xlab="value")
hist(burninsamples_M1$s1, probability =TRUE, main="Posterior samples of s1", xlab="value")
hist(burninsamples_M1$s1, probability =TRUE, main="Posterior samples of s2", xlab="value")
hist(burninsamples_M1$s3, probability =TRUE, main="Posterior samples of s3", xlab="value")
hist(burninsamples_M1$s4, probability =TRUE, main="Posterior samples of s4", xlab="value")
hist(burninsamples_M1$p1, probability =TRUE, main="Posterior samples of p1", xlab="value")
hist(burninsamples_M1$p2, probability =TRUE, main="Posterior samples of p2", xlab="value")
hist(burninsamples_M1$tau, probability =TRUE, main="Posterior samples of tau", xlab="value")


