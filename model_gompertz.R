model{
  
  # Priors
  r ~ dunif(0, 1)
  K ~ dunif(0, 10000)
  lambda1 ~ dunif(0, 10000)
  
  # Abundance
  for(i in 1:nSites){
    N[i,1] ~ dpois(lambda1)
    for(t in 2:nYears){
      N[i,t] ~ dpois(lambda[i,t-1]) # t-1 is just for accounting
      lambda[i,t-1] <- N[i,t-1] * exp(r * (1 - log(N[i,t-1] + 1) / log(K + 1)))
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
      lp[i,t] <- p.mu # + p.b1*prcp7day[i,t]
    }
  }
  
  # Priors: detection
  p.mean ~ dunif(0.1, 0.9)
  p.mu <- log(p.mean/(1-p.mean))
  p.b1 ~ dnorm(0, 0.37)I(-3,3)
  p.b2 ~ dnorm(0, 0.37)I(-3,3)
  
}