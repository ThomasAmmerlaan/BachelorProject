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

statsZhat <- matrix(c(Zhat0,Zhat1,Zhat0-log(Z_0),Zhat1-log(Z_1)),nrow=5000,ncol=4)

###Convergence

#livepoints
livepoint <- c(2,5,10,25,50,75,100,150,200,250,500,750,1000,1250,1500,1750,2000,
               2250,2500)

#statistics we want to keep
empirical_logevidence_M0 <-numeric()
empirical_logevidencevar_M0 <-numeric()
empirical_logevidenceerror_M0 <-numeric()
empirical_information_M0 <-numeric()

empirical_logevidence_M1 <-numeric()
empirical_logevidencevar_M1 <-numeric()
empirical_logevidenceerror_M1 <-numeric()
empirical_information_M1 <-numeric()

#nested sample runs
for (i in 1:length(livepoint)){
  logevidenceruni_M0 <- numeric()
  logerrorruni_M0 <- numeric()
  informationruni_M0 <- numeric()
  
  logevidenceruni_M1 <- numeric()
  logerrorruni_M1 <- numeric()
  informationruni_M1 <- numeric()
  for (j in 1:25){
    set.seed(Sys.time())
    sampler_M0 <- ernest_sampler(loglik_M0, prior_M0, nlive = livepoint[i])
    sampler_M1 <- ernest_sampler(loglik_M1, prior_M1, nlive = livepoint[i])
    results_M0 <- generate(sampler_M0)
    results_M1 <- generate(sampler_M1)
    
    logevidenceruni_M0 <- c(logevidenceruni_M0, results_M0$log_evidence)
    logerrorruni_M0 <- c(logerrorruni_M0, results_M0$log_evidence_err)
    informationruni_M0 <- c(informationruni_M0,results_M0$information) 
    
    logevidenceruni_M1 <- c(logevidenceruni_M1, results_M1$log_evidence)
    logerrorruni_M1 <- c(logerrorruni_M1, results_M1$log_evidence_err)
    informationruni_M1 <- c(informationruni_M1,results_M1$information) 
  }
  empirical_logevidence_M0 <- c(empirical_logevidence_M0,
                                mean(logevidenceruni_M0))
  empirical_logevidencevar_M0 <- c(empirical_logevidencevar_M0,
                                   var(logevidenceruni_M0))
  empirical_logevidenceerror_M0 <- c(empirical_logevidenceerror_M0, 
                                     mean(logerrorruni_M0))
  empirical_information_M0 <- c(empirical_information_M0, 
                                mean(informationruni_M0))
  
  empirical_logevidence_M1 <- c(empirical_logevidence_M1,
                                mean(logevidenceruni_M1))
  empirical_logevidencevar_M1 <- c(empirical_logevidencevar_M1,
                                   var(logevidenceruni_M1))
  empirical_logevidenceerror_M1 <- c(empirical_logevidenceerror_M1, 
                                     mean(logerrorruni_M1))
  empirical_information_M1 <- c(empirical_information_M1, 
                                mean(informationruni_M1))
}

#plots
plot(livepoint,empirical_logevidence_M0[1:length(livepoint)],
     xlab="number of livepoints",
     ylab="empirical logevidence",main="log(Z0)")
abline(h=log(Z_0),col="red")  
legend(x="bottomright", legend=c("theoretical value"),col=c("red"),lty=1)

plot(livepoint,empirical_logevidence_M1[1:length(livepoint)],
     xlab="number of livepoints",
     ylab="empirical logevidence",main="log(Z1)")
abline(h=log(Z_1),col="blue")  
legend(x="bottomright", legend=c("theoretical value"),col=c("blue"),lty=1)

plot(livepoint[1:length(livepoint)], 
     empirical_logevidencevar_M0[1:length(livepoint)],
     xlab="number of livepoints", 
     ylab="empirical variance of (log)evidence",main="M0")
points(livepoint[1:length(livepoint)],
       empirical_information_M0[1:length(livepoint)]/livepoint[1:length(livepoint)],
       type='l', col="red")
legend(x="topright", legend=c("emperical H"),col=c("red"),lty=1)

plot(livepoint[1:length(livepoint)], 
     empirical_logevidencevar_M1[1:length(livepoint)],
     xlab="number of livepoints", 
     ylab="empirical variance of (log)evidence",main="M1")
points(livepoint[1:length(livepoint)],
       empirical_information_M1[1:length(livepoint)]/livepoint[1:length(livepoint)],
       type='l', col="blue")
legend(x="topright", legend=c("emperical H"),col=c("blue"),lty=1)
