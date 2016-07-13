model{
  
  # Priors
  for(i in 1:nSites){
    site.r[i] ~ dnorm(0, tau.r)
  }
  tau.r <- pow(sigma.r, -2)
  sigma.r ~ dunif(0, 10)
  alpha.r ~ dunif(0, 10)
  
  # priors on initial abundance
  alpha.0 ~ dnorm(0, 0.001)
  for(i in 1:nSites) {
    site.0[i] ~ dnorm(0, tau.0)
  }
  tau.0 <- pow(sigma.0, -2)
  sigma.0 ~ dunif(0, 100)
  
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
    b[h] ~ dnorm(mu.b[h], tau.b[h])
    mu.b[h] ~ dnorm(0, 0.01)
    sigma.b[h] ~ dunif(0, 10)
    tau.b[h] <- pow(sigma.b[h], -2)
  }
  
  # Abundance
  for(i in 1:nSites){
    N[i,1] ~ dpois(lambda.0[i])
    log(lambda.0[i]) <- alpha.0 + site.0[i]
    for(t in 2:nYears){
      N[i,t] ~ dpois(lambda[i,t-1]) # t-1 is just for accounting
      lambda[i,t-1] <- N[i,t-1] * exp(r[i,t-1] * (1 - log(N[i,t-1] + 1) / log(K[i] + 1))) + rho[i,t-1] # + siteLength[i,j]
      }
  }
  
  # regession on intrinsic growth rate
  for(i in 1:nSites){
    for(t in 1:(nYears-1)){
  # r[i,t] <- alpha.r + site.r[i] + b[1,i]*fall.prcp[i,t] + b[2,i]*winter.prcp[i,t] + b[3,i]*spring.prcp[i,t] + b[4,i]*fall.tmean[i,t] + b[5,i]*winter.tmean[i,t] + b[6,i]*spring.tmean[i,t]
  r[i,t] <- alpha.r + site.r[i] + b[1]*fall.prcp[i,t] + b[2]*winter.prcp[i,t] + b[3]*spring.prcp[i,t] + b[4]*fall.tmean[i,t] + b[5]*winter.tmean[i,t] + b[6]*spring.tmean[i,t]
    }
  }
  
  # regression on carrying capacity (add drainage area)
  for(i in 1:nSites) {
    K[i] <- K.0 + K.site[i]
    K.site[i] ~ dnorm(0, tau.k) # prior on K.site
  }
  K.0 ~ dunif(0, 10000)
  # hyperpriors for K
  tau.k <- 1 / (sigma.k * sigma.k)
  sigma.k ~ dunif(0, 100)
  
  # add density independent factors
  for(i in 1:nSites){
    for(t in 1:(nYears-1)){
      log(rho[i,t]) <- ln.iota + eps.rho[i,t]
    }
  }
  
  # density-independent priors
  ln.iota ~ dnorm(0, 0.001)
  iota <- exp(ln.iota)
  for(i in 1:nSites){
    for(t in 1:(nYears-1)){
      eps.rho[i,t] ~ dnorm(0, tau.eps.rho)
    }
  }
  tau.eps.rho <- 1 / (sigma.eps.rho * sigma.eps.rho)
  sigma.eps.rho ~ dunif(0, 100)
  
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