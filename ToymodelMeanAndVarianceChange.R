###libraries used
library(ernest)
library(ggplot2)
library(ggdist)
library(patchwork)

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