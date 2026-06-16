###libraries used
library(ernest)
library(ggplot2)
library(ggdist)
library(patchwork)
library(pracma)

#simple Gaussian 

logf <- function(x){
  x <- x[1]
  log(exp(-(x^2)))
}
  
#priors
uniform <- create_uniform_prior(names = c("mu","dummy"), 
                                lower = -1, upper = 1)
normal <- create_normal_prior(names=c("mu","dummy"))

#number of livepoints we run
livepoint <- c(2,5,10,25,50,100,250,500,750,1000,1250,1500,1750,2000,
               2500,3000,3500,4000,4500,5000,10000)
#The statistics we want to hold  
empirical_logevidenceunif <-numeric()
empirical_logevidence_sdunif <-numeric()
empirical_logevidenceerrorunif <-numeric()
empirical_informationunif <-numeric()
  
empirical_logevidencenormal <-numeric()
empirical_logevidence_sdnormal <-numeric()
empirical_logevidenceerrornormal <-numeric()
empirical_informationnormal <-numeric()

#nest sample runs for the livepoints 
for (i in 1:length(livepoint)){
  logevidenceruniunif <- numeric()
  logerrorruniunif <- numeric()
  informationruniunif <- numeric()
  
  logevidenceruninormal <- numeric()
  logerrorruninormal <- numeric()
  informationruninormal <- numeric()
  for (j in 1:25){
    set.seed(Sys.time())
    sampleruniform <- ernest_sampler(logf, uniform, nlive = livepoint[i])
    samplernormal <- ernest_sampler(logf, normal, nlive = livepoint[i])
    resultsuniform <- generate(sampleruniform)
    resultsnormal <- generate(samplernormal)
    
    logevidenceruniunif <- c(logevidenceruniunif, resultsuniform$log_evidence)
    logerrorruniunif <- c(logerrorruniunif, resultsuniform$log_evidence_err)
    informationruniunif <- c(informationruniunif,resultsuniform$information) 
    
    logevidenceruninormal <- c(logevidenceruninormal, resultsnormal$log_evidence)
    logerrorruninormal <- c(logerrorruninormal, resultsnormal$log_evidence_err)
    informationruninormal <- c(informationruninormal,resultsnormal$information) 
  }
  empirical_logevidenceunif <- c(empirical_logevidenceunif,
                                 mean(logevidenceruniunif))
  empirical_logevidence_sdunif <- c(empirical_logevidence_sdunif,
                                    var(logevidenceruniunif))
  empirical_logevidenceerrorunif <- c(empirical_logevidenceerrorunif, 
                                      mean(logerrorruniunif))
  empirical_informationunif <- c(empirical_informationunif, 
                                 mean(informationruniunif))
  
  empirical_logevidencenormal <- c(empirical_logevidencenormal,
                                   mean(logevidenceruninormal))
  empirical_logevidence_sdnormal <- c(empirical_logevidence_sdnormal,
                                      var(logevidenceruninormal))
  empirical_logevidenceerrornormal <- c(empirical_logevidenceerrornormal, 
                                        mean(logerrorruninormal))
  empirical_informationnormal <- c(empirical_informationnormal, 
                                   mean(informationruninormal))
}
plotl <-function(i){
  plot(livepoint[i:21], empirical_logevidence_sdunif[i:21]
       ,xlab="number of livepoints", 
       ylab="emperical variance",ylim=c(0,empirical_logevidence_sdnormal[i]))
  points(livepoint[i:21],empirical_informationunif[i:21]/livepoint[i:21]
         ,type='l')
  points(livepoint[i:21],empirical_informationnormal[i:21]/livepoint[i:21]
         ,type='l', col='red')
  points(livepoint[i:21], empirical_logevidence_sdnormal[i:21], col='red')
  points(seq(2, 10000, by = 0.1),(1/2*(log(3))-1/2+1/6)/seq(2, 10000, by = 0.1)
         ,type='l',lty=2, col='blue')
  points(seq(2, 10000, by = 0.1),
         (log(2)-log(sqrt(pi)*erf(1))-1/2+1/(exp(1)*sqrt(pi)*erf(1)))/seq(2, 10000, by = 0.1)
         ,type='l',lty=2,col='green')
}
plotl(1)
legend(x="topright", legend=c("uniform", "normal","theoretical normal","theoretical uniform"),
       col=c('black','red','blue','green'),lty=c(1,1,2,2))

plot(livepoint, exp(empirical_logevidencenormal), xlab = "number of livepoints", ylab="emperical evidence using normal prior")
abline(h=1/sqrt(3), col="red")
legend(x='bottomright', legend=c("theoretical value"),col="red", lty=1)
plot(livepoint,exp(empirical_logevidenceunif),xlab="number of livepoints", ylab="emperical evidence using uniform prior")
abline(h=1/2*sqrt(pi)*(erf(1)), col="blue")
legend(x='bottomright', legend=c("theoretical value"),col="blue", lty=1)

###For mean change 
n <- 100
tau <- 50
Y <- c(rnorm(tau),rnorm(n-tau, mean= 1))

#Theoretical evidence values
Z_2 <- 1/(sqrt(2*pi)^n*sqrt(n+1))*exp(-1/2*(sum((Y-mean(Y))^2)+n*mean(Y)^2/(n+1)))
Z_1 <- 0
for (i in 1:(n-1)){
  part1 <- 1/(sqrt(2*pi)^i*sqrt(i+1))
  part2 <- exp(-1/2*(sum((Y[1:i]-mean(Y[1:i]))^2)+i*mean(Y[1:i])^2/(i+1)))
  part3 <- 1/(sqrt(2*pi)^(n-i)*sqrt((n-i)+1))
  part4 <- exp(-1/2*(sum((Y[(i+1):n]-mean(Y[(i+1):n]))^2)+(n-i)*mean(Y[(i+1):n])^2/((n-i)+1)))
  Z_1 <- Z_1+part1*part2*part3*part4
}
Z_1 <- Z_1/(n-1)

#loglikelihoods
loglik_M0 <- function(theta){
  mu <- theta[1]
  -n/2*log(2*pi)-1/2*sum((Y[1:n]-mu)^2)
}
loglik_M1 <- function(theta){
  mu1 <-theta[1]
  mu2 <-theta[2]
  tau <-as.integer(theta[3])
  -n/2*log(2*pi)-1/2*sum((Y[1:tau]-mu1)^2)-1/2*sum((Y[(tau+1):n]-mu2)^2)
}

#priors
prior_M0 <- create_normal_prior(names=c("mu1","dummy"))

norm_M1 <- create_normal_prior(names=c("mu1","mu2"))
unif_M1 <- create_uniform_prior(names="tau", lower=1, upper=n-1)
prior_M1 <- norm_M1+unif_M1

#livepoints
livepoint <- c(2,5,10,25,50,100,250,500,750,1000,1250,1500,1750,2000,2500,3000)

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
  empirical_logevidence_M0 <- c(empirical_logevidence_M0,mean(logevidenceruni_M0))
  empirical_logevidencevar_M0 <- c(empirical_logevidencevar_M0,var(logevidenceruni_M0))
  empirical_logevidenceerror_M0 <- c(empirical_logevidenceerror_M0, mean(logerrorruni_M0))
  empirical_information_M0 <- c(empirical_information_M0, mean(informationruni_M0))
  
  empirical_logevidence_M1 <- c(empirical_logevidence_M1,mean(logevidenceruni_M1))
  empirical_logevidencevar_M1 <- c(empirical_logevidencevar_M1,var(logevidenceruni_M1))
  empirical_logevidenceerror_M1 <- c(empirical_logevidenceerror_M1, mean(logerrorruni_M1))
  empirical_information_M1 <- c(empirical_information_M1, mean(informationruni_M1))
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

plot(livepoint[6:16], empirical_logevidencevar_M0[6:length(livepoint)],
     xlab="number of livepoints", 
     ylab="empirical variance of (log)evidence",main="M0")
points(livepoint[6:16],
       empirical_information_M0[6:length(livepoint)]/livepoint[6:16],
       type='l', col="red")
legend(x="topright", legend=c("emperical H"),col=c("red"),lty=1)

plot(livepoint[6:16], empirical_logevidencevar_M1[6:length(livepoint)],
     xlab="number of livepoints", 
     ylab="empirical variance of (log)evidence",main="M1")
points(livepoint[6:16],
       empirical_information_M1[6:length(livepoint)]/livepoint[6:16],
       type='l', col="blue")
legend(x="topright", legend=c("emperical H"),col=c("blue"),lty=1)