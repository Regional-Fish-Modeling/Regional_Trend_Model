### Southern-range brook trout ADULT abundance model adapted from Kanno et al. (2015) Freshwater Biology
# NP Hitt <nhitt@usgs.gov> 1 November 2015

# Load working directory & libraries
#setwd("C:/Users/nhitt/Desktop/BKTSouth_Nmix")
#memory.limit(size=229355) # 7x the baseline memory limit (32765) 
library(reshape2)
# library(rjags)
library(jagsUI)
library(plyr)
library(ggplot2)
library(arm)
library(boot)
library(coda)
# library(R2WinBUGS)

# Fish count data
load("Data_FishCountAr.RData")
# Seasonal climate data
load("Data_SeasonalClimateStd.RData") # Standardized to site-specific mean=0 and sd=1
# Site covariate data
load("Data_SiteCovsStd.RData") # Standardized to dataset-level mean=0 and sd=1
# Detection covariate data
load("Data_DetectionCovsStd.RData") # Standardized to dataset-level mean=0 and sd=1

# Set up data structure
# nSites=326; nYears=34; nCovs=6; nPasses=3

# Adding 1 fish to abundances to avoid log(0) problem
ADUFish <- ADUFish + 1

# prcp7day=prcp7day.std, sampday=sampday.std

rownames(ADUFish) <- NULL
ADUFish
dim(ADUFish)

# Try on subset of data
ADUFish <- ADUFish[1:30,,]
dim(ADUFish)

# Check sites
# ADUFish[3,,]
# ADUFish[7,,]

nSites <- dim(ADUFish)[1]
nYears <- dim(ADUFish)[2]
nCovs <- 6

dat <- list(nSites=dim(ADUFish)[1], nYears=dim(ADUFish)[2], y=ADUFish, nCovs=nCovs, fall.prcp=FallPrcpStd, winter.prcp=WinterPrcpStd, spring.prcp=SpringPrcpStd, 
            fall.tmean=FallTmeanStd, winter.tmean=WinterTmeanStd, spring.tmean=SpringTmeanStd, 
            prcp7day=prcp7day.std, sampday=sampday.std, elev=elev.std)

# # Bundle data for Adult model
# dat <- list(nSites=nSites, nYears=nYears, nCovs=nCovs, y=ADUFish,
#             fall.prcp=FallPrcpStd, winter.prcp=WinterPrcpStd, spring.prcp=SpringPrcpStd, 
#             fall.tmean=FallTmeanStd, winter.tmean=WinterTmeanStd, spring.tmean=SpringTmeanStd, 
#             prcp7day=prcp7day.std, sampday=sampday.std, elev=elev.std)

# Set initial values
# alpha=array(runif(nSites,-5,5), dim=(nSites))
# p.b1=array(rnorm(1,0), 1),
# p.b2=array(rnorm(1,0), 1),

# Make decent starting values for N
N.init <- apply(dat$y, 1:2, sum)
Ni.max <- apply(N.init, 1, max, na.rm = TRUE)
Ni.max[which(Ni.max == -Inf)] <- max(Ni.max, na.rm = TRUE)
k <- which(is.na(N.init), arr.ind=TRUE)
N.init[k] <- Ni.max[k[,1]]
N.init <- (N.init + 1) * 2

inits <- function() list(N=N.init, #N=array(500, dim=c(nSites, nYears)),
                         p.mean=0.5)

parameters <- c("N", "K.0", "sigma.k", "alpha.r", "sigma.r", "b", "alpha.0", "sigma.0", "sigma.eps.rho", "iota", "p.mean", "N.region") #, "p")


# MCMC settings
ni <- 7000
nt <- 3
nb <- 5000
nc <- 3


# bugs.dir <- getwd() # "C:/Program Files/WinBUGS14/"

start.time = Sys.time()         # Set timer 
# Call BUGS from R 

out <- jags(data = dat, inits = inits, parameters.to.save = parameters, model.file = "model_gompertz.R", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time

# Summarize posteriors
print(out, dig = 3)

sum1 <- out$summary
# write.csv(sum1,'poisson.csv')

out$mean$K
out$mean$K

# Find which parmeters, if any, have Rhat > 1.1
which(out$summary[, c("Rhat")] > 1.1)

# Or see what max Rhat value is
max(out$summary[, c("Rhat")])

# str(out)
