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
library(tidyr)

#---------- set model -----------
model <- "gompertz"
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

#---------- Check covariate relationships ----------
str(FallPrcpStd)

ADU_catch <- as.data.frame(apply(ADUFish, 1:2, sum, na.rm = TRUE), stringsAsFactors = FALSE)
ADU_catch$Site <- dimnames(ADUFish)[1][[1]]
k <- which(is.na(ADUFish[ , , 1]), arr.ind=TRUE) # which have pass 1 = NA
ADU_catch[k] <- NA 
ADU_long <- tidyr::gather(ADU_catch, Year, Catch, -Site)

FallPrcpStd2 <- as.data.frame(FallPrcpStd, stringsAsFactors = FALSE)
colnames(FallPrcpStd) <- 1982:2015
FallPrcpStd2$Site <- dimnames(FallPrcpStd)[[1]]
FallPrcpStd_long <- tidyr::gather(FallPrcpStd2, Year, FallPrcpStd, -Site)

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * 0.5)
}

Pairs <- dplyr::left_join(ADU_long, FallPrcpStd_long)
pairs(Pairs[, 3:ncol(Pairs)], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, cex.labels = 2, font.labels = 2)

stocks <- data_frame(
  time = as.Date('2009-01-01') + 0:9,
  X = rnorm(10, 0, 1),
  Y = rnorm(10, 0, 2),
  Z = rnorm(10, 0, 4)
)

gather(stocks, stock, price, -time)
stocks %>% gather(stock, price, -time)
#---------- Set up data structure ---------
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
ADUFish <- ADUFish[ , 1:33, ] # 2015 data seem off for the covariates
dim(ADUFish)

# Try on subset of data
if(testing ==TRUE) {
  ADUFish <- ADUFish[sample(1:nrow(ADUFish), size = 40, replace = FALSE), , ]
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
# N.init <- apply(dat$y, 1:2, sum, na.rm = TRUE)
# Ni.max <- apply(N.init, 1, max, na.rm = TRUE)
# Ni.max[which(Ni.max == -Inf)] <- max(Ni.max, na.rm = TRUE)
# k <- which(is.na(dat$y[ , , 1]), arr.ind=TRUE)
# N.init[k] <- Ni.max[k[,1]]
# N.init <- (N.init + 1) * 2
# N.init <- ceiling(N.init) # / SurveyLength[1:nrow(N.init), ]

fillN <- function(N, p, lambda) { # Need to replace NAs in initial N with values, adjust others by detection probability
  for (i in 1:nrow(N)) {
    notNA <- which(!is.na(N[i, ]))
    size <- N[i, notNA] 
    size[size==0] <- 0.5
    addTo <- rbinom(length(notNA), size=size, prob=p)
    N[i, notNA] <- N[i, notNA] + addTo
    for (t in 1:ncol(N))
      if (is.na(N[i,t])) {
        if (t == ncol(N) || is.na(N[i, t+1])) {
          if (t==1) {
            N[i, t] <- rbinom(1, size=lambda, prob=p) #If no previous or subsequent value, use supplied default
          }
          else {
            N[i, t] <- N[i, t-1] #If no subsequent value, use previous
          }
        }
        else if (t==1)
          N[i, t] <- N[i, t+1] #If no previous value, use subsequent
        else
          N[i, t] <- ceiling((N[i, t-1] + N[i, t+1])/2) #Otherwise, use mean
      }
  }
  return(N)
}

# N.init <- fillN(N = dat$y[ , , 1], p = 0.3, lambda = median(dat$y[ , , 1], na.rm = TRUE))
# head(N.init)
# summary(N.init)

# simple initial N
N.init <- array(2 * max(dat$y, na.rm = TRUE), dim=c(nSites, nYears))
inits <- function() list(N = N.init,
                         lambda.0=runif(1, 10, 80),
                         r=runif(1, 0.1, 0.4), 
                         sigma.nu=runif(1), 
                         K=runif(1, 50, 100), 
                         iota=runif(1),
                         p.mean = runif(1, 0.5, 0.8))

if(model == "overdispersion") {
parameters <- c("N", "K.0", "sigma.k", "alpha.r", "sigma.r", "b", "alpha.0", "sigma.0", "sigma.eps.rho", "iota", "p.mean", "N.region", "b1.p", "K.0") #, "p")
} else {
  parameters <- c("N", "K.0", "sigma.k", "alpha.r", "sigma.r", "b", "alpha.0", "sigma.0", "p.mean", "N.region", "mu.b", "b1.p", "K.0", "sigma.sy") #, "p")
}

parameters <- c("N", "K", "p.mean", "N.region", "b1.p", "r", "lambda.0", "sigma.nu", "iota") #, "p")

# MCMC settings
ni <- 3000
nt <- 1
nb <- 1
nc <- 3


# bugs.dir <- getwd() # "C:/Program Files/WinBUGS14/"
mod <- ifelse(model == "overdispersion", "Code/model_gompertz_od.R", "Code/model_gompertz.R")

start.time = Sys.time() # Set timer

# Call BUGS from R 
out <- jags(data = dat, inits = inits, parameters.to.save = parameters, model.file = "Code/simple_gompertz.R", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time

# update to increase length of chains if need more - use for determining number needed for convergence
# out <- update(out, n.iter = 20000)

# Traceplots for parameters least likely to mix well or converge
if(model == "overdispersion") {
  jagsUI::traceplot(out, parameters = c("alpha.0", "alpha.r", "sigma.0", "sigma.k", "sigma.r", "sigma.b", "K.0", "sigma.eps.rho"))
  } else {
    jagsUI::traceplot(out, parameters = c("alpha.0", "alpha.r", "sigma.0", "sigma.k", "sigma.r", "b", "K.0"))
  }


jagsUI::traceplot(out, parameters = c("K", "p.mean", "b1.p", "r", "lambda.0", "sigma.nu", "iota"))

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
