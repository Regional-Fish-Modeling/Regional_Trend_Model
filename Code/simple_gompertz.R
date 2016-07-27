model{
  
  # Priors
  
  # Priors: random slopes (uncorrelated currently)
  # for(h in 1:nCovs){
  #   for(i in 1:nSites){
  #     b[h,i] ~ dnorm(mu.b[h], tau.b[h])
  #   }
  #   # hyperpriors
  #   mu.b[h] ~ dnorm(0, 0.001)
  #   tau.b[h] <- pow(sigma.b[h], -2)
  #   sigma.b[h] ~ dunif(0, 100)
  # }
  
  # fixed slops
  for(h in 1:nCovs) {
    b[h] ~ dnorm(0, 0.001)
  }
  
  lambda.0 ~ dunif(0, 500)
  r ~ dunif(0, 50)
  K ~ dunif(0, 100)
  sigma.nu ~ dunif(0, 10)
  tau.nu <- 1 / (sigma.nu * sigma.nu)
  iota ~ dunif(0, 100)
  
  for(t in 2:nYears) {
    nu[t-1] ~ dnorm(0, tau.nu)
  }
  
  # Abundance
  for(i in 1:nSites){
    N[i,1] ~ dpois(lambda.0)
    # log(lambda.0[i]) <- alpha.0 + site.0[i] + log(survey.length[i,1])
    for(t in 2:nYears){
      N[i,t] ~ dpois(lambda[i,t-1]) # t-1 is just for accounting
      log(lambda[i,t-1]) <- log(N[i,t-1] * exp(nu[t-1] + r * (1 - log(N[i,t-1] + 1) / log(K + 1))) + iota) + log(survey.length[i,t])
    }
  }
  
  # Detection
  for(i in 1:nSites){
    for(t in 1:nYears){
      y[i,t,1] ~ dbin(p[i,t], N[i,t])
      y[i,t,2] ~ dbin(p[i,t], N[i,t] - y[i,t,1])
      y[i,t,3] ~ dbin(p[i,t], N[i,t] - y[i,t,1] - y[i,t,2])
      
      p[i,t] <- 1/(1 + exp(-lp.lim[i,t]))
      lp.lim[i,t] <- min(999, max(-999, lp[i,t]))
      lp[i,t] <- p.mu + b1.p*prcp7day[i,t]
    }
  }
  
  # Priors: detection
  p.mean ~ dunif(0.1, 0.9)
  p.mu <- log(p.mean/(1-p.mean))
  b1.p ~ dnorm(0, 0.37)I(-3,3)
  # p.b2 ~ dnorm(0, 0.37)I(-3,3)
  
  # Derived parameters
  #for(i in 1:nSites) {
  for(t in 1:nYears) {
    N.region[t] <- sum(N[ ,t]) 
  }
  # }
  
}