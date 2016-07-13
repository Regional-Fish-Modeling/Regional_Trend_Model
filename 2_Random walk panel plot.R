###############################################################
# Extract pop'n average coefficents
# alpha0 = fixed intercept
# gamma = year affect common across all sites
# eta = site-specific intercept
T <- nYears
N <- nSites
# 
# FixedSiteTrend <- array(NA,c(out$n.sim,T,N) )
# dim(FixedSiteTrend)
# 
# for(i in 1:out$n.sim){
# for(t in 1:T){
# 	for(n in 1:N){
# 		FixedSiteTrend[i,t,n] <- exp((out$sims.list$alpha0[i] + out$sims.list$eta[i,n]) + out$sims.list$gamma[i,t] )
# 	}
# }	
# }
# 
# dim(FixedSiteTrend)
# head(FixedSiteTrend)
# 
# FixedSiteTrend2 <- matrix(NA, nrow=T, ncol=N)
# for(t in 1:T){
# 	for(n in 1:N){
# 		FixedSiteTrend2[t,n] <- mean(FixedSiteTrend[,t,n])
# 	}
# }
# 
# head(FixedSiteTrend2)
# summary(FixedSiteTrend[,1,1])
# 
# FixedSiteTrendLCI <- matrix(NA, nrow=T, ncol=N)
# for(t in 1:T){
# 	for(n in 1:N){
# 		FixedSiteTrendLCI[t,n] <- quantile(FixedSiteTrend[,t,n],0.025)
# 	}
# }	
# 
# head(FixedSiteTrendLCI)
# 
# FixedSiteTrendUCI <- matrix(NA, nrow=T, ncol=N)
# for(t in 1:T){
# 	for(n in 1:N){
# 		FixedSiteTrendUCI[t,n] <- quantile(FixedSiteTrend[,t,n],0.975)
# 	}
# }	
# 
# head(FixedSiteTrendUCI)


###########################
# SITE-SPECIFIC TRENDS
###########################
# u = site intercept
# xi = site trend

# Site=specific trend
# SiteTrend <- array(NA, c(out$n.sims, T,N) )
# for(t in 1:T){
# 	for(n in 1:N){
# 		SiteTrend[,t,n] <- exp(out$sims.list$xi[,n,t] + out$sims.list$alpha[,n])
# 	}
# }	

SiteTrend <- array(NA, c(dim(out$sims.list$N)[1], T, N) )
for(t in 1:T){
	for(n in 1:N){
		SiteTrend[,t,n] <- out$sims.list$N[,n,t]
	}
}
# SiteTrend <- out$sims.list$N


SiteTrend2 <- matrix(NA, nrow=T, ncol=N)
for(t in 1:T){
	for(n in 1:N){
		SiteTrend2[t,n] <- median(SiteTrend[,t,n])
	}
}	

head(SiteTrend2)
dim(SiteTrend2)

# Site=specific trend
SiteTrendLCI <- matrix(NA, nrow=T, ncol=N)
for(t in 1:T){
	for(n in 1:N){
		SiteTrendLCI[t,n] <- quantile(SiteTrend[,t,n],0.025)
	}
}	


head(SiteTrendLCI)

SiteTrendUCI <- matrix(NA, nrow=T, ncol=N)
for(t in 1:T){
	for(n in 1:N){
		SiteTrendUCI[t,n] <- quantile(SiteTrend[,t,n],0.975)
	}
}	


# Obtain posterior means for z's
# z.est <- numeric()
# for(n in 1:N){
# 	z.est[n] <- mean(out$sims.list$z[,n])
# }

abund <- matrix(NA, nrow=T, ncol=N)
for(t in 1:T){
  for(n in 1:N){
    abund[t,n] <- median(out$sims.list$N[,n,t])
  }
}	




# Plot subset of data

########### PLOT
res <- 6

sitePlot <- c(1:N)

name_figure <- 'AdultTrend.png'
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1
size.text = 1.0
x.label = 'Year'
y.label = 'Abundance'

plotYear <- 1982:2015


nf <- layout(matrix( c(1:(N)),nrow=5,ncol=2,byrow=T),  TRUE) 
layout.show(nf)
par(mar=c(0.5,0.5,0.5,0.5),oma=c(3,3.5,0,1),mai=c(0.1,0.1,0.1,0) )	

for(i in 1:N){

# min(SiteTrendLCI[,i])
  # max(SiteTrendUCI[,i])
  # max(dat$ratio[dat$state2==i])
  #+1e-06
  
plot(plotYear, abund[,1], axes=F, ylim=c(0,max(SiteTrendUCI[,i])), ylab='', xlab='', type='n')

i.for <- order( plotYear )
i.back <- order( plotYear , decreasing = TRUE )

x.polygon <- c( plotYear[i.for] , plotYear[i.back] )
y.polygon <- c( SiteTrendLCI[,i][i.for] , SiteTrendUCI[,i][i.back] )

polygon( x.polygon , y.polygon , col = "#A6CEE3" , border = NA )

# x.polygon2 <- c( plotYear[i.for] , plotYear[i.back] )
# y.polygon2 <- c( FixedSiteTrendLCI[,i][i.for] , FixedSiteTrendUCI[,i][i.back] )

# polygon( x.polygon2 , y.polygon2 , col = "gray" , border = NA )
# Add fitted lines
points(plotYear,SiteTrend2[,i], cex=0.8, pch=16,type='l',lty=1)
# points(plotYear,FixedSiteTrend2[,i], pch=22,type='l', lty=2) 
# Add estiamted N
points(plotYear, abund[,i], pch=16, cex=0.8)

# text(1991, 0.8, paste('z =',round(z.est[i],2)),font=1,col='black' )
text(1985,max(SiteTrendUCI[,i])-5, i)

axis(side=1, cex.axis=0.8, at=plotYear,  tck=-0.01, mgp=c(0,0.2,0) ) 
axis(side=2,cex.axis=0.8,font=1 ,tck=-0.01, mgp=c(0,0.15,0), las=1) 
mtext(x.label, line = 1, side = 1, cex = size.text, outer=T)
mtext(y.label, line = 1.5, side = 2, cex = size.text, outer=T)
box()
}

par(def.par)
dev.off()


#----------- Overall Regional Trend -------------

# convert to mean density?

N.region <- data.frame(Year = plotYear, Mean = out$mean$N.region, LCRI = out$q2.5$N.region, UCRI = out$q97.5$N.region, stringsAsFactors = FALSE)

ggplot(N.region, aes(Year, Mean)) + geom_line() + geom_point() + geom_ribbon(aes(ymin = LCRI, ymax = UCRI), alpha = 0.3) + ylab("Regional abundance") + theme_bw()

ggsave(filename = "Output/Figures/Regional_Abundance.png")

