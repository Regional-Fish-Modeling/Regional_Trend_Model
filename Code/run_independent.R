### Southern-range brook trout ADULT abundance model adapted from Kanno et al. (2015) Freshwater Biology
# NP Hitt <nhitt@usgs.gov> 1 November 2015

# Load working directory & libraries
#setwd("C:/Users/nhitt/Desktop/BKTSouth_Nmix")
#memory.limit(size=229355) # 7x the baseline memory limit (32765) 
# library(reshape2)
# library(rjags)
library(jagsUI)
# library(plyr)
# library(ggplot2)
# library(arm)
# library(boot)
library(coda)
library(readr)
library(dplyr)

#---------- set model -----------
testing <- TRUE


#---------- Load Data ----------

# Fish count data
load("Data/Data_FishCountAr.RData")
# Seasonal climate data
load("Data/Data_SeasonalClimateStd.RData") # Standardized to site-specific mean=0 and sd=1
# Site covariate data
load("Data/Data_SiteCovsStd.RData") # Standardized to dataset-level mean=0 and sd=1
# Detection covariate data
load("Data/Data_DetectionCovsStd.RData") # Standardized to dataset-level mean=0 and sd=1

SurveyLength <- read_csv("Data/BKT_SURVEY_LENGTH.csv")

# Set up data structure
# nSites=326; nYears=34; nCovs=6; nPasses=3

# Make survey length sites and years line up with fish data
SurveyLength <- data.frame(SiteID = dimnames(ADUFish)[[1]], stringsAsFactors = FALSE) %>%
  dplyr::left_join(SurveyLength) %>%
  dplyr::select(-SiteID)

# Convert survey length into matrix for JAGS
SurveyLength <- as.matrix(SurveyLength)

# Adding 1 fish to abundances to avoid log(0) problem
# ADUFish <- ADUFish + 1

# prcp7day=prcp7day.std, sampday=sampday.std

rownames(ADUFish) <- NULL
ADUFish <- ADUFish[ , 1:33, ]
dim(ADUFish)

# Try on subset of data
if(testing ==TRUE) {
  ADUFish <- ADUFish[1:30, 1:33, ]
}

# dim(ADUFish)

# Check sites
# ADUFish[3,,]
# ADUFish[7,,]

nSites <- dim(ADUFish)[1]
nYears <- dim(ADUFish)[2]
nCovs <- 6

dat <- list(nSites=dim(ADUFish)[1], nYears=dim(ADUFish)[2], y=ADUFish, nCovs=nCovs, fall.prcp=FallPrcpStd, winter.prcp=WinterPrcpStd, spring.prcp=SpringPrcpStd, 
            fall.tmean=FallTmeanStd, winter.tmean=WinterTmeanStd, spring.tmean=SpringTmeanStd, 
            prcp7day=prcp7day.std, sampday=sampday.std, elev=elev.std, survey.length = SurveyLength)

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
N.init <- apply(dat$y, 1:2, sum, na.rm = TRUE)
Ni.max <- apply(N.init, 1, max, na.rm = TRUE)
Ni.max[which(Ni.max == -Inf)] <- max(Ni.max, na.rm = TRUE)
k <- which(is.na(dat$y[ , , 1]), arr.ind=TRUE)
N.init[k] <- Ni.max[k[,1]]
N.init <- (N.init + 1) * 2
N.init <- ceiling(N.init) # / SurveyLength[1:nrow(N.init), ]

# simple initial N
# N.init <- array(2 * max(dat$y, na.rm = TRUE), dim=c(nSites, nYears)

inits <- function() list(N = N.init,
                         p.mean = runif(1, 0.4, 0.8))

  parameters <- c("N", "b", "alpha.0", "sigma.0", "p.mean", "N.region", "sigma.b", "mu.b", "b1.p") #, "p")

# MCMC settings
ni <- 1000
nt <- 1
nb <- 1
nc <- 3

start.time = Sys.time() # Set timer

# Call BUGS from R 
out <- jags(data = dat, inits = inits, parameters.to.save = parameters, model.file = "Code/model_independent.R", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time

# update to increase length of chains if need more - use for determining number needed for convergence
# out <- update(out, n.iter = 20000)

# Traceplots for parameters least likely to mix well or converge
  jagsUI::traceplot(out, parameters = c("alpha.0", "sigma.0", "sigma.b"))

# Whisker plots
whiskerplot(out, parameters = c("alpha.0", "alpha.r", "sigma.0", "sigma.k", "sigma.r", "sigma.b", "sigma.eps.rho"))

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
