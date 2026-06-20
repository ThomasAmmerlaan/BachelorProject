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
    
    prop$tau <- rnorm(1, theta$tau, step_tau)
    
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

samples_M0 <- mh_sampler_M0(Y, n_iter = 100000)
llmax_M0 <- samples_M0[[6]]
samples_M1 <- mh_sampler_M1(Y, n_iter = 500000)
llmax_M1 <- samples_M1[[12]]

samples_M0 <- as.data.frame(samples_M0[1:5])
samples_M1 <- as.data.frame(samples_M1[1:11])
burninsamples_M0 <- samples_M0[(burnin+1):100000,]
burninsamples_M1 <- samples_M1[(burnin+1):500000,]
bic0 <- -2 * llmax_M0 + 5 * log(n)
bic1 <- -2 * llmax_M1 + 11 * log(n)

logBF10_MCMC <- (bic0 - bic1)/2
logBF10_MCMC