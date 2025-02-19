### Streamflow sensitivity regimes of alpine catchments and their relationship with elevation, temperature, and glacier cover ###
# Code to make plots

#NA functions

mean.na <- function(x,...){if (all(is.na(x))){NA} else {mean(x, na.rm = TRUE)}}
sd.na <- function(x,...){if (all(is.na(x))){NA} else {sd(x, na.rm = TRUE)}}
sum.na <- function(x,...){if (all(is.na(x))){NA} else {sum(x, na.rm = TRUE)}}
min.na <- function(x,...){if (all(is.na(x))){NA} else {min(x, na.rm = TRUE)}}
max.na <- function(x,...){if (all(is.na(x))){NA} else {max(x, na.rm = TRUE)}}
which.max.na <- function(x,...){if (all(is.na(x))){NA} else {which.max(x)}}
which.min.na <- function(x,...){if (all(is.na(x))){NA} else {which.min(x)}}

#sign (+1, 0, -1)
sig.fun <- function(x){sign(x[2]-x[1])}

#difference
dif.fun <- function(x){(x[2]-x[1])}

#Mann-Kendall
mann.kendall <- function(t,x)
{
  #omit NA
  x <- x[which(is.finite(x))]
  t <- t[which(is.finite(x))]
  
  #sum of sign combinations
  Q <- sum(combn(x, m = 2, FUN = sig.fun, simplify = T))
  
  #variance, p-value
  n <- length(x)
  va <- 1/18 * n *(n - 1) * (2 *n + 5)
  p <- 1 - pnorm(abs(Q), mean = 0, sd = sqrt(va))
  
  #difference combations: time and value
  dt.all  <- combn(t, m = 2, FUN = dif.fun, simplify = T)
  dx.all  <- combn(x, m = 2, FUN = dif.fun, simplify = T)
  
  #slope and intercept
  slope <- median(dx.all / dt.all)
  intercept <- median(x) - slope * median(t)
  
  #output
  return(list(S = Q/sqrt(va), p_one_side = p,
              slope = slope, intercept=intercept, median=median(x), no=length(x)))
}

#moving average; circular
ma.filtc <- function(x, window)
{
  m <- floor(window/2)
  y <- filter(x, rep(1/window, window), sides = 2, method = "convolution", circular = TRUE)
  return(y[1:length(x)])
}


library(relaimpo)
library(car)

# plotting discharge discharge, temperature and precipitation annual variability for 25 Swiss catchments

out_T <- read.table(file="D:/R/weekly_streamflow_sensitivities/Weekly_Sensitivity_to_T.txt")
out_P <- read.table(file="D:/R/weekly_streamflow_sensitivities/Weekly_Sensitivity_to_P.txt")
out_Q <- read.table(file="D:/R/weekly_streamflow_sensitivities/Weekly_Sensitivity_to_Q.txt")
out_isotherm <- read.table(file="D:/R/weekly_streamflow_sensitivities/0deg_isotherm_regime.txt")

#path <- "/Users/mweiler/Markus/Research Projects/KHR/Abflussdaten/"

elim <- c(min(out_T$elev_ezg),max(out_T$elev_ezg))

etit <- "Elevation(m)"

height <- as.numeric(out_T$elev_ezg)
glac <- as.numeric(out_T$Glacier)

elev <- seq(500,3000,100)
glaciers <- seq(0,60,2.5)

# Important Parameter

nr <- length(out_T$Q_Nr)

# linear interpolation

ws <- 52
in1 <- 13
out_T.p <- out_T[,in1:(in1+ws-1)]
out_T.coeff <- out_T[,(in1+ws):(in1+2*ws-1)]
out_T.var <- out_T[,(in1+2*ws):(in1+3*ws-1)]
out_P.p <- out_P[,in1:(in1+ws-1)]
out_P.coeff <- out_P[,(in1+ws):(in1+2*ws-1)]
out_P.var <- out_P[,(in1+2*ws):(in1+3*ws-1)]
out_Q.p <- out_Q[,in1:(in1+ws-1)]
out_Q.coeff <- out_Q[,(in1+ws):(in1+2*ws-1)]
out_Q.var <- out_Q[,(in1+2*ws):(in1+3*ws-1)]
#a0 <- out_Q[,(in1+3*ws):(in1+4*ws-1)]

tcoef <- matrix(NA, nrow=nr, ncol=52)
pcoef <- matrix(NA, nrow=nr, ncol=52)
qcoef <- matrix(NA, nrow=nr, ncol=52)


############ Compare sensitivity coefficients of catchments (Figure 1) #############

library("RColorBrewer")
#col <- colorRampPalette(brewer.pal(n=9, name='YlGnBu'))(30)
col <- terrain.colors(30, alpha = 1)
palette(col)

xx <- seq(1,52)*7-3
xtick = c(0,31,59,90,120,151,181,212,243,273,304,334,365)
xtlab = c("J","F","M","A","M","J","J","A","S","O","N","D",NA)
zero <- rep(0, length(xx))

pdf(file="D:/R/weekly_streamflow_sensitivities/Figure 1 new.pdf",width=8, height=9)

#par( mfrow = c(3,3),mar=c(c(3, 4, 1, 1) + 0.1), cex=0.9)
#par(xaxs = "i", yaxs = "i")
par(mfrow = c(3,3),
    mar = c(2.5, 3, 1.5, 0.5),  # Increased top margin slightly
    oma = c(0.5, 0.5, 0.5, 0),  # Added small outer margin at top
    mgp = c(2, 0.7, 0),         # Axis label positioning
    cex = 0.9,
    xaxs = "i", 
    yaxs = "i")

movavg <-1
# absolute sensitivity
for (i in 1:nr) {tcoef[i,] <- ma.filtc(as.numeric(out_T.coeff[i,]),movavg)}  
for (i in 1:nr) {pcoef[i,] <- ma.filtc(as.numeric(out_P.coeff[i,]),movavg)}
for (i in 1:nr) {qcoef[i,] <- ma.filtc(as.numeric(out_Q.coeff[i,]),movavg)}

ylim1 <- c(-0.8,2.6)
plot(xx,tcoef[1,], typ="l", ylim=ylim1, col=col[(height[1]-500)*.01], xlab="", ylab = expression(b[T]*" [(mm/d)/K]"), xaxt = "n", lwd = 2,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)
for (i in 1:nr) {points(xx,tcoef[i,], typ="l",col=col[ceiling((height[i]-500)*.01)], lwd = 2)}
legend(1,ylim1[2], seq(500,3000,length.out=11),col=col, fill=seq(1,25,length.out=11), horiz=F, title="H [m a.s.l.]", bty="n", xjust=-0.05, cex=0.6)

ylim1 <- c(-0.1,1.0)
plot(xx,pcoef[1,], typ="l", ylim=ylim1, col=col[(height[1]-500)*.01], xlab="", ylab = expression(b[P]*" [-]"), xaxt = "n", lwd = 2 ,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)
for (i in 1:nr) {points(xx,pcoef[i,], typ="l",col=col[ceiling((height[i]-500)*.01)], lwd = 2)}
#legend(183,ylim1[2], seq(500,3000,length.out=6),col=col, fill=seq(1,25,length.out=6), horiz=TRUE, title="Elevation [m a.s.l.])", bty="n", xjust=0.5, cex=0.8)

ylim1 <- c(-0.1,1.4)
plot(xx,qcoef[1,], typ="l", ylim=ylim1, col=col[(height[1]-500)*.01], xlab="", ylab = expression(b[Q]*" [-]"), xaxt = "n", lwd = 2,  xaxs="i" )
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)
for (i in 1:nr) {points(xx,qcoef[i,], typ="l",col=col[ceiling((height[i]-500)*.01)], lwd = 2)}
#legend(183,ylim1[2]+0.05, seq(500,3000,length.out=5),col=col, fill=seq(1,25,length.out=5), horiz=TRUE, title="Elevation [m a.s.l.]", bty="n", xjust=0.5, cex=0.5)

# relative sensitivity
for (i in 1:nr) {tcoef[i,] <- ma.filtc(as.numeric(out_T.coeff[i,]/out_Q.p[i,]*100),movavg)}   # relative sensitivity
for (i in 1:nr) {pcoef[i,] <- ma.filtc(as.numeric(out_P.coeff[i,]/out_Q.p[i,]*100),movavg)}
for (i in 1:nr) {qcoef[i,] <- ma.filtc(as.numeric(out_Q.coeff[i,]/out_Q.p[i,]*100),movavg)}

ylim1 <- c(-20,25)
plot(xx,tcoef[1,], typ="l", ylim=ylim1, col=col[(height[1]-500)*.01], xlab="", ylab = expression(S[T]*" [%/K]"), xaxt = "n", lwd = 2,  xaxs="i" , yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)
for (i in 1:nr) {points(xx,tcoef[i,], typ="l",col=col[ceiling((height[i]-500)*.01)], lwd = 2)}
#legend(183,ylim1[2], seq(500,3000,length.out=6),col=col, fill=seq(1,25,length.out=6), horiz=TRUE, title="Elevation [m a.s.l.]", bty="n", xjust=0.5, cex=0.8)

ylim1 <- c(-5,50)
plot(xx,pcoef[1,], typ="l", ylim=ylim1, col=col[(height[1]-500)*.01], xlab="", ylab = expression(S[P]*" [%/(mm/d)]"), xaxt = "n", lwd = 2 ,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)
for (i in 1:nr) {points(xx,pcoef[i,], typ="l",col=col[ceiling((height[i]-500)*.01)], lwd = 2)}
#print(max(pcoef[i,]))
#points(xx,pcoef[25,], typ="l",col=col[(height[25]-500)*.01], lwd = 2)
#legend(183,ylim1[2], seq(500,3000,length.out=6),col=col, fill=seq(1,25,length.out=6), horiz=TRUE, title="Elevation [m a.s.l.]", bty="n", xjust=0.5, cex=0.8)

ylim1 <- c(-10,200)
plot(xx,qcoef[1,], typ="l", ylim=ylim1, col=col[(height[1]-500)*.01], xlab="", ylab = expression(S[Q]*" [%/(mm/d)]"), xaxt = "n", lwd = 2 ,  xaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)
for (i in 1:nr) {points(xx,qcoef[i,], typ="l",col=col[ceiling((height[i]-500)*.01)], lwd = 2)}
# legend(183,ylim1[2], seq(500,3000,length.out=6),col=col, fill=seq(1,25,length.out=6), horiz=TRUE, title="Elevation [m a.s.l.]", bty="n", xjust=0.5, cex=0.8)
legend(1,ylim1[2], seq(500,3000,length.out=11),col=col, fill=seq(1,25,length.out=11), horiz=F, title="Elevation [m a.s.l.]", bty="n", xjust=-2.5, cex=0.7)

### explained variation
for (i in 1:nr) {tcoef[i,] <- ma.filtc(as.numeric(out_T.var[i,]),movavg)}  
for (i in 1:nr) {pcoef[i,] <- ma.filtc(as.numeric(out_P.var[i,]),movavg)}
for (i in 1:nr) {qcoef[i,] <- ma.filtc(as.numeric(out_Q.var[i,]),movavg)}

clim_coef = tcoef + pcoef + qcoef
print(quantile(clim_coef))

ylim1 <- c(0,1)
plot(xx,tcoef[1,], typ="l", ylim=ylim1, col=col[(height[1]-500)*.01], xlab="", ylab="Relative importance T [-]", xaxt = "n", lwd = 2,  xaxs="i" , yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)
for (i in 1:nr) {points(xx,tcoef[i,], typ="l",col=col[ceiling((height[i]-500)*.01)], lwd = 2)}
#legend(183,ylim1[2], seq(500,3000,length.out=6),col=col, fill=seq(1,25,length.out=6), horiz=TRUE, title="Elevation [m a.s.l.]", bty="n", xjust=0.5, cex=0.8)

ylim1 <- c(0,1)
plot(xx,pcoef[1,], typ="l", ylim=ylim1, col=col[(height[1]-500)*.01], xlab="", ylab="Relative importance P [-]", xaxt = "n", lwd = 2 ,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)
for (i in 1:nr) {points(xx,pcoef[i,], typ="l",col=col[ceiling((height[i]-500)*.01)], lwd = 2)}
#legend(183,ylim1[2], seq(500,3000,length.out=6),col=col, fill=seq(1,25,length.out=6), horiz=TRUE, title="Elevation [m a.s.l.]", bty="n", xjust=0.5, cex=0.8)

ylim1 <- c(0,1)
plot(xx,qcoef[1,], typ="l", ylim=ylim1, col=col[(height[1]-500)*.01], xlab="", ylab="Relative importance Q(t-1) [-]", xaxt = "n", lwd = 2,  xaxs="i" , yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)
for (i in 1:nr) {points(xx,qcoef[i,], typ="l",col=col[ceiling((height[i]-500)*.01)], lwd = 2)}
#legend(183,ylim1[2], seq(500,3000,length.out=6),col=col, fill=seq(1,25,length.out=6), horiz=TRUE, title="Elevation [m a.s.l.]", bty="n", xjust=0.5, cex=0.8)

dev.off() 

# model for relating glacier cover to elevation

#using only the data from the 25 catchments
mod <- lm(glac[which(glac > 0)] ~ height[which(glac > 0)])
glac_elev <- mod$coefficients[2]*elev+mod$coefficients[1]
glac_elev[which(glac_elev <= 0)] <- 0
plot(height,glac, xlim=c(500,3000))
lines(elev,glac_elev)

# using the data from all 150 km2 catchments in CH
# EZG <- read.csv(file="/Users/mweiler/Markus/Research Projects/KHR/Abflussdaten/EZG150_CH/EZG150_Glacier.csv", sep=";")
# mod <- lm(EZG$GC_2003[which(EZG$GC_2003 > 0)] ~ EZG$Hmean[which(EZG$GC_2003 > 0)])
# glac_elev <- mod$coefficients[2]*elev+mod$coefficients[1]
# glac_elev[which(glac_elev <= 0)] <- 0
# points(EZG$Hmean,EZG$GC_2003, col="red")
# lines(elev,glac_elev, col="red")


# apply linear models to look for elevation dependence

abs <- 1 #FALSE #TRUE			# if abs 1 than plot absolute sensitivity, if 2 then relative, if 3 than elasticity

# 1 moving average in time 
if(abs == 1){																	#absolute sensitivity
  for (i in 1:nr) {tcoef[i,] <- ma.filtc(as.numeric(out_T.coeff[i,]),movavg)}  
  for (i in 1:nr) {pcoef[i,] <- ma.filtc(as.numeric(out_P.coeff[i,]),movavg)}
}

if(abs == 2){	 
  for (i in 1:nr) {tcoef[i,] <- ma.filtc(as.numeric(out_T.coeff[i,]/out_Q.p[i,]*100),movavg)}   # relative sensitivity
  for (i in 1:nr) {pcoef[i,] <- ma.filtc(as.numeric(out_P.coeff[i,]/out_Q.p[i,]*100),movavg)}
}

if(abs == 3){	 
  for (i in 1:nr) {tcoef[i,] <- ma.filtc(as.numeric(out_T.coeff[i,]/out_Q.p[i,]),movavg)}   # elacticity
  for (i in 1:nr) {pcoef[i,] <- ma.filtc(as.numeric(out_P.coeff[i,]/out_Q.p[i,]*out_P.p[i,]),movavg)}
}


# loop to calculate LM depending on elevation and glacier cover for each week
# select interaction between glacier and elevation

inter = 0 # if 1, then with interaction


mod_T <- data.frame(r2=rep(NA,ws), p=rep(NA,ws), inter_c=rep(NA,ws), inter_p=rep(NA,ws), elev_c=rep(NA,ws), glac_c=rep(NA,ws), elev_p=rep(NA,ws), glac_p=rep(NA,ws), elev_r=rep(NA,ws), glac_r=rep(NA,ws),inter_r=rep(NA,ws)) 
mod_P <- data.frame(r2=rep(NA,ws), p=rep(NA,ws), inter_c=rep(NA,ws), inter_p=rep(NA,ws), elev_c=rep(NA,ws), glac_c=rep(NA,ws), elev_p=rep(NA,ws), glac_p=rep(NA,ws), elev_r=rep(NA,ws), glac_r=rep(NA,ws),inter_r=rep(NA,ws)) 

T_elev <- matrix(data=NA, nrow=ws, ncol=length(elev))					# partials for glac and elevation
T_glac <- matrix(data=NA, nrow=ws, ncol=length(glaciers))	
T_elev50 <- matrix(data=NA, nrow=ws, ncol=length(elev))	

P_elev <- matrix(data=NA, nrow=ws, ncol=length(elev))
P_glac <- matrix(data=NA, nrow=ws, ncol=length(glaciers))
P_elev50 <- matrix(data=NA, nrow=ws, ncol=length(elev))	

for (i in 1:ws) 
{
  pred <- as.numeric(tcoef[,i])
  height <- as.numeric(out_T$elev_ezg)
  glac <- as.numeric(out_T$Glacier)
  if(inter == 1){mod <- lm(pred ~ height * glac)}  else {mod <- lm(pred ~ height + glac)}
  T_elev[i,] <- predict(mod, data.frame(height=elev, glac=rep(0,length(elev))),se.fit = TRUE)$fit   # partial for glacier = 0
  T_glac[i,] <- predict(mod, data.frame(height=rep(2500,length(glaciers)), glac=glaciers),se.fit = TRUE)$fit   # partial for elevation = 2500
  T_elev50[i,] <- predict(mod, data.frame(height=elev, glac=glac_elev),se.fit = TRUE)$fit   # partial for glacier mit mittlerer Gletscher??nderung pro H??he
  
  mod_T$elev_c[i] <- mod$coefficients[2]
  mod_T$glac_c[i] <- mod$coefficients[3]
  res <- summary(mod)
  mod_T$elev_p[i] <- res$coefficients[2,4]
  mod_T$glac_p[i] <- res$coefficients[3,4]
  mod_T$p[i]  <- pf(res$fstatistic[1], res$fstatistic[2], res$fstatistic[3], lower.tail = FALSE)
  mod_T$r2[i] <- res$r.squared
  varimp <- calc.relimp(mod)
  mod_T$elev_r[i] <- 	varimp$lmg[1]
  mod_T$glac_r[i] <- 	varimp$lmg[2]
  if(inter == 1){	 
    mod_T$inter_c[i] <- mod$coefficients[4]
    mod_T$inter_p[i] <- res$coefficients[4,4]
    mod_T$inter_r[i] <- 	varimp$lmg[3]
  }
}

for (i in 1:ws) 
{
  pred <- as.numeric(pcoef[,i])
  height <- as.numeric(out_P$elev_ezg)
  if(inter == 1){mod <- lm(pred ~ height * glac)}  else {mod <- lm(pred ~ height + glac)}
  P_elev[i,] <- predict(mod, data.frame(height=elev, glac=rep(0,length(elev))),se.fit = TRUE)$fit   # partial for glacier = 0
  P_glac[i,] <- predict(mod, data.frame(height=rep(2500,length(glaciers)), glac=glaciers),se.fit = TRUE)$fit   # partial for elevation = 2500
  P_elev50[i,] <- predict(mod, data.frame(height=elev, glac=glac_elev),se.fit = TRUE)$fit  # partial for glacier mit mittlerer Gletscher??nderung pro H??he
  
  mod_P$elev_c[i] <- mod$coefficients[2]
  mod_P$glac_c[i] <- mod$coefficients[3]
  res <- summary(mod)
  mod_P$elev_p[i] <- res$coefficients[2,4]
  mod_P$glac_p[i] <- res$coefficients[3,4]
  mod_P$p[i]  <- pf(res$fstatistic[1], res$fstatistic[2], res$fstatistic[3], lower.tail = FALSE)
  mod_P$r2[i] <- res$r.squared
  varimp <- calc.relimp(mod)
  mod_P$elev_r[i] <- 	varimp$lmg[1]
  mod_P$glac_r[i] <- 	varimp$lmg[2]
  if(inter == 1){	 
    mod_P$inter_c[i] <- mod$coefficients[4]
    mod_P$inter_p[i] <- res$coefficients[4,4]
    mod_P$inter_r[i] <- 	varimp$lmg[3]
  }
}


############ Dependence of absolute sensitivites on elevation (Figure 2) #############

xtick = c(0,31,59,90,120,151,181,212,243,273,304,334,365)
xtlab = c("J","F","M","A","M","J","J","A","S","O","N","D",NA)
xx <- seq(1,52)*7-3
rb <- (hsv(h = c(seq(.02,.02,length.out=10), 0, seq(.6,.6,length.out=10)), v = c(seq(1,1,length.out=10), 1,seq(1,1,length.out=10)), s=c(seq(1,0.1,length.out=10), 0, seq(0.1,1,length.out=10))))  #blue white red
yg <- (hsv(h = seq(.175,.6,length.out=20), v = seq(1,1,length.out=20), s=seq(0.2,1,length.out=20)))  # yellow green blue
yl <- (hsv(h = c(seq(.175,0,length.out=10),seq(1,1-0.175,length.out=10)), v = seq(1,1,length.out=20), s=c(seq(0,1,length.out=10),seq(1,1,length.out=10))))  # yellow to lila
yg <- (hsv(h = c(seq(.16,.17,length.out=5),seq(.17,.75,length.out=15)), v = seq(1,1,length.out=20), s=c(seq(0.0,1,length.out=5),seq(1,1,length.out=15))))  # yellow green blue
bv	<- (hsv(h = c(seq(.6,.6,length.out=10), seq(.6,.8,length.out=10)), v = c(seq(1,1,length.out=10), seq(1,1,length.out=10)), s=c(seq(0.0,.8,length.out=10), seq(0.8,1,length.out=10)))) 					#white over blue to lila
wbb <- rev(rgb(seq(0,1,length.out=11), seq(0,1,length.out=11), seq(0,1,length.out=11)))			# white to black
wb <- rev(rgb(seq(0,1,length.out=11), seq(0,1,length.out=11), seq(1,1,length.out=11)))  # white to dark blue
wbt <- rev(rgb(seq(1,1,length.out=11), seq(1,1,length.out=11), seq(1,1,length.out=11), seq(0.1,0.9,length.out=11)))			# white to black but variable transparent

elim1 <- c(500,3210) #c(500,3710)
elim2 <- c(0,64)

abs <- 1

pdf(file="D:/R/weekly_streamflow_sensitivities/Figure 2 new.pdf",width=8, height=8)

#par( mfrow = c(2,2),mar=c(c(2, 4, 1, 1) + 0.1), cex=0.8)
#par(mfrow = c(2, 2), mar = c(2, 4, 1, 1) + 0.1, oma = c(0, 0, 0, 3), cex = 0.8)
layout_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4, ncol = 2, byrow = TRUE)
layout(layout_matrix, heights = c(3, 1.5, 3, 1.5))
#par(mar = c(1, 4, 1, 1) + 0.1, oma = c(0, 0, 0, 3), cex = 0.8)
par(oma = c(0, 0, 0, 3), cex = 0.8)

if(abs == 1){ maxt <- 1.5
maxp <- 1.0
title_legend_T <- expression(b[T]*" [(mm/d)/K]")
title_legend_P <- expression(b[P]*" [-]")}

if(abs == 2){ maxt <- 25
maxp <- 30
title_legend_T <- expression(S[T]*" [%/K]")
title_legend_P <- expression(b[P]*" [%/(mm/d)]")}

palette(rb)
img <- T_elev
temp <- which(img > maxt)
if(length(temp) > 0){img[temp] <- maxt}
par(mar = c(0, 4, 1, 1) + 0.1)
image(xx, elev, img, zlim=c(-maxt,maxt), ylim=elim1, xlim=c(0,365), col=rb, xlab="", ylab="Elevation [m a.s.l.]", xaxt = "n" )  
legend(1,elim1[2]*0.95, seq(-maxt,maxt,length.out=11),col=rb, fill=seq(1,21,length.out=11), horiz=FALSE, title=title_legend_T, bty="n", xjust=0, cex=0.7)
#sig <- matrix(rep((mod_T$elev_p < 0.05),length(elev)),nrow=52,ncol=length(elev))
#image(xx, elev, sig, col=c(rgb(.3,.3,.3,.3),rgb(1,1,1,0)), add=T, zlim=c(0,1))
sigmod <- which(mod_T$elev_p < 0.05)
points(xx[sigmod],zero[sigmod]+elim1[2]*0.98, pch=15, col="gray")
#filled.contour(xx, elev, img,levels=seq(-maxt,maxt,length.out=21),ylim=elim1, xlim=c(0,365),col=rb, xlab="", ylab="Elevation [m a.s.l.]", xaxt = "n")
abline(elim1[2]*0.96,0)
yy_iso <- out_isotherm$isotherm.1[1:52]
yy_iso[yy_iso>3050] = NA
yy_iso2 <- out_isotherm$isotherm.2[1:52]
lines(xx, yy_iso, col="grey50", lwd=2, lty=5)
lines(xx, yy_iso2, col="grey50", lwd=2, lty=2)

palette(wb)
img <- P_elev
temp <- which(img > maxp)
if(length(temp) > 0){img[temp] <- maxp}
par(mar = c(0, 4, 1, 1) + 0.1)
#image(xx, elev, img, zlim=c(0,maxp), ylim=elim1, xlim=c(0,365),col=wb, xlab="", ylab="Elevation [m a.s.l.]", xaxt = "n" )
image(xx, elev, img, zlim=c(0, maxp), ylim=elim1, xlim=c(0, 365), col=wb, xlab="", ylab="", xaxt="n", yaxt="n")
legend(1,elim1[2]*0.95, seq(0,maxp,length.out=6),col=wb, fill=seq(1,11,length.out=6), horiz=FALSE, title=title_legend_P, bty="n", xjust=0, cex=0.7)
#sig <- matrix(rep((mod_P$elev_p < 0.05),length(elev)),nrow=52,ncol=length(elev))
#image(xx, elev, sig, col=c(rgb(.3,.3,.3,.3),rgb(1,1,1,0)), add=T, zlim=c(0,1))
sigmod <- which(mod_P$elev_p < 0.05)
points(xx[sigmod],zero[sigmod]+elim1[2]*0.98, pch=15, col="gray")
abline(elim1[2]*0.96,0)
lines(xx, yy_iso, col="grey50", lwd=2, lty=5)
lines(xx, yy_iso2, col="grey50", lwd=2, lty=2)

# Define custom y-tick positions and labels
y_ticks <- c(3339, 2944, 2549, 2154, 1759, 1364, 969, 574)  # Y-tick positions
#y_labels <- c(8.3, 5.8, 3.3, 0.8,-1.7,-4.2,-6.7)  # Y-axis labels #<- 1:7
y_labels <- c(-6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0)  # Y-axis labels #<- 1:7 
axis(side = 4, at = y_ticks, labels = y_labels, col.axis = "black", las = 0)  #
mtext("Temperature [°C]", side = 4, line = 3, col = "black", cex=0.8)  #

par(xaxs = "i", yaxs = "i")
par(mar = c(3, 4, 1, 1) + 0.1)
plot(xx,mod_T$elev_c*-1000/5, col="grey40", lwd=2, type="l",
     ylim=c(-0.2,0.05), xlim=c(0,365), xlab="", ylab=expression(b[TH]*" [*/K]"), xaxt = "n",
     panel.first = abline(0, 0, col = "grey50", lwd = 1, lty = 2))
# Add x-axis
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))

par(mar = c(3, 4, 1, 1) + 0.1)
plot(xx, mod_P$elev_c*-1000/5, col="grey40", lwd=2, type="l",
     ylim=c(-0.05,0.1), xlim=c(0,365), xlab="", ylab="", xaxt="n", yaxt="n",
     panel.first = abline(0, 0, col = "grey50", lwd = 1, lty = 2))
# Add x-axis
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
# Add y-axis on the right side
axis(side = 4)
# Add y-axis label on the right side
mtext(expression(b[PH]*" [*/K]"), side = 4, line = 3, cex = par("cex"))

palette(rb)
img <- T_glac
temp <- which(img > maxt)
if(length(temp) > 0){img[temp] <- maxt}
par(mar = c(0, 4, 1, 1) + 0.1)
image(xx, glaciers, img, zlim=c(-maxt,maxt), ylim=elim2, xlim=c(0,365), col=rb, xlab="", ylab="Glacier Cover [%]", xaxt = "n" )  
#sig <- matrix(rep((mod_T$glac_p < 0.05),length(glaciers)),nrow=52,ncol=length(glaciers))
#image(xx, glaciers, sig, col=c(rgb(.3,.3,.3,.3),rgb(1,1,1,0)), add=T, zlim=c(0,1))
sigmod <- which(mod_T$glac_p < 0.05)
points(xx[sigmod],zero[sigmod]+elim2[2]*0.98, pch=15, col="gray")
abline(elim2[2]*0.96,0)

palette(wb)
img <- P_glac
temp <- which(img > maxp)
if(length(temp) > 0){img[temp] <- maxp}
#image(xx, glaciers, img, zlim=c(0,maxp), ylim=elim2, xlim=c(0,365),col=wb, xlab="", ylab="Glacier Cover [%]", xaxt = "n" )
par(mar = c(0, 4, 1, 1) + 0.1)
image(xx, glaciers, img, zlim=c(0,maxp), ylim=elim2, xlim=c(0,365),col=wb, xlab="", ylab="", xaxt = "n" , yaxt="n")
axis(4)
mtext("Glacier Cover [%]", side = 4, line = 3, col = "black", cex=0.8)  #
#sig <- matrix(rep((mod_P$glac_p < 0.05),length(glaciers)),nrow=52,ncol=length(glaciers))
#image(xx, glaciers, sig, col=c(rgb(.3,.3,.3,.3),rgb(1,1,1,0)), add=T, zlim=c(0,1))
sigmod <- which(mod_P$glac_p < 0.05)
points(xx[sigmod],zero[sigmod]+elim2[2]*0.98, pch=15, col="gray")
abline(elim2[2]*0.96,0)

# Define custom y-tick positions and labels
#y_ticks <- c(55, 42, 29, 16, 3)  # Y-tick positions
#y_labels <- c(-6, -4, -2, 0, 2)  # Y-axis labels #<- 1:7
#axis(side = 4, at = y_ticks, labels = y_labels, col.axis = "black", las = 1)  #
#mtext("Temperature [°C]", side = 4, line = 3, col = "black", cex=0.8)  #

par(xaxs = "i", yaxs = "i")
par(mar = c(3, 4, 1, 1) + 0.1)
plot(xx,mod_T$glac_c, col="grey40", lwd=2, type="l",
     ylim=c(-0.02,0.04), xlim=c(0,365), xlab="", ylab=expression(b[TG]*" [*/%]"), xaxt = "n",
     panel.first = abline(0, 0, col = "grey50", lwd = 1, lty = 2))
# Add x-axis
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))

par(mar = c(3, 4, 1, 1) + 0.1)
plot(xx, mod_P$glac_c, col="grey40", lwd=2, type="l",
     ylim=c(-0.01,0.01), xlim=c(0,365), xlab="", ylab="", xaxt="n", yaxt="n",
     panel.first = abline(0, 0, col = "grey50", lwd = 1, lty = 2))
# Add x-axis
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
# Add y-axis on the right side
axis(side = 4)
# Add y-axis label on the right side
mtext(expression(b[PG]*" [*/%]"), side = 4, line = 3, cex = par("cex"))

dev.off() 


### REDO WITH RELATIVE SENSTIVITIES

# Apply linear models to look for elevation dependence

abs <- 2 #FALSE #TRUE			# if abs 1 than plot absolute sensitivity, if 2 then relative, if 3 then elasticity

# 1 moving average in time 
if(abs == 1){																	#absolute sensitivity
  for (i in 1:nr) {tcoef[i,] <- ma.filtc(as.numeric(out_T.coeff[i,]),movavg)}  
  for (i in 1:nr) {pcoef[i,] <- ma.filtc(as.numeric(out_P.coeff[i,]),movavg)}
}

if(abs == 2){	 
  for (i in 1:nr) {tcoef[i,] <- ma.filtc(as.numeric(out_T.coeff[i,]/out_Q.p[i,]*100),movavg)}   # relative sensitivity
  for (i in 1:nr) {pcoef[i,] <- ma.filtc(as.numeric(out_P.coeff[i,]/out_Q.p[i,]*100),movavg)}
}

if(abs == 3){	 
  for (i in 1:nr) {tcoef[i,] <- ma.filtc(as.numeric(out_T.coeff[i,]/out_Q.p[i,]),movavg)}   # elasticity
  for (i in 1:nr) {pcoef[i,] <- ma.filtc(as.numeric(out_P.coeff[i,]/out_Q.p[i,]*out_P.p[i,]),movavg)}
}

# loop to calculate LM depending on elevation and glacier cover for each week
# select interaction between glacier and elevation

inter = 0 # if 1, then with interaction

mod_T <- data.frame(r2=rep(NA,ws), p=rep(NA,ws), inter_c=rep(NA,ws), inter_p=rep(NA,ws), elev_c=rep(NA,ws), glac_c=rep(NA,ws), elev_p=rep(NA,ws), glac_p=rep(NA,ws), elev_r=rep(NA,ws), glac_r=rep(NA,ws),inter_r=rep(NA,ws)) 
mod_P <- data.frame(r2=rep(NA,ws), p=rep(NA,ws), inter_c=rep(NA,ws), inter_p=rep(NA,ws), elev_c=rep(NA,ws), glac_c=rep(NA,ws), elev_p=rep(NA,ws), glac_p=rep(NA,ws), elev_r=rep(NA,ws), glac_r=rep(NA,ws),inter_r=rep(NA,ws)) 

T_elev <- matrix(data=NA, nrow=ws, ncol=length(elev))					# partials for glac and elevation
T_glac <- matrix(data=NA, nrow=ws, ncol=length(glaciers))	
T_elev50 <- matrix(data=NA, nrow=ws, ncol=length(elev))	

P_elev <- matrix(data=NA, nrow=ws, ncol=length(elev))
P_glac <- matrix(data=NA, nrow=ws, ncol=length(glaciers))
P_elev50 <- matrix(data=NA, nrow=ws, ncol=length(elev))	

for (i in 1:ws) 
{
  pred <- as.numeric(tcoef[,i])
  height <- as.numeric(out_T$elev_ezg)
  glac <- as.numeric(out_T$Glacier)
  if(inter == 1){mod <- lm(pred ~ height * glac)}  else {mod <- lm(pred ~ height + glac)}
  T_elev[i,] <- predict(mod, data.frame(height=elev, glac=rep(0,length(elev))),se.fit = TRUE)$fit   # partial for glacier = 0
  T_glac[i,] <- predict(mod, data.frame(height=rep(2500,length(glaciers)), glac=glaciers),se.fit = TRUE)$fit   # partial for elevation = 2500
  T_elev50[i,] <- predict(mod, data.frame(height=elev, glac=glac_elev),se.fit = TRUE)$fit   # partial for glacier mit mittlerer Gletscher??nderung pro H??he
  
  mod_T$elev_c[i] <- mod$coefficients[2]
  mod_T$glac_c[i] <- mod$coefficients[3]
  res <- summary(mod)
  mod_T$elev_p[i] <- res$coefficients[2,4]
  mod_T$glac_p[i] <- res$coefficients[3,4]
  mod_T$p[i]  <- pf(res$fstatistic[1], res$fstatistic[2], res$fstatistic[3], lower.tail = FALSE)
  mod_T$r2[i] <- res$r.squared
  varimp <- calc.relimp(mod)
  mod_T$elev_r[i] <- 	varimp$lmg[1]
  mod_T$glac_r[i] <- 	varimp$lmg[2]
  if(inter == 1){	 
    mod_T$inter_c[i] <- mod$coefficients[4]
    mod_T$inter_p[i] <- res$coefficients[4,4]
    mod_T$inter_r[i] <- 	varimp$lmg[3]
  }
}

for (i in 1:ws) 
{
  pred <- as.numeric(pcoef[,i])
  height <- as.numeric(out_P$elev_ezg)
  if(inter == 1){mod <- lm(pred ~ height * glac)}  else {mod <- lm(pred ~ height + glac)}
  P_elev[i,] <- predict(mod, data.frame(height=elev, glac=rep(0,length(elev))),se.fit = TRUE)$fit   # partial for glacier = 0
  P_glac[i,] <- predict(mod, data.frame(height=rep(2500,length(glaciers)), glac=glaciers),se.fit = TRUE)$fit   # partial for elevation = 2500
  P_elev50[i,] <- predict(mod, data.frame(height=elev, glac=glac_elev),se.fit = TRUE)$fit  # partial for glacier mit mittlerer Gletscher??nderung pro H??he
  
  mod_P$elev_c[i] <- mod$coefficients[2]
  mod_P$glac_c[i] <- mod$coefficients[3]
  res <- summary(mod)
  mod_P$elev_p[i] <- res$coefficients[2,4]
  mod_P$glac_p[i] <- res$coefficients[3,4]
  mod_P$p[i]  <- pf(res$fstatistic[1], res$fstatistic[2], res$fstatistic[3], lower.tail = FALSE)
  mod_P$r2[i] <- res$r.squared
  varimp <- calc.relimp(mod)
  mod_P$elev_r[i] <- 	varimp$lmg[1]
  mod_P$glac_r[i] <- 	varimp$lmg[2]
  if(inter == 1){	 
    mod_P$inter_c[i] <- mod$coefficients[4]
    mod_P$inter_p[i] <- res$coefficients[4,4]
    mod_P$inter_r[i] <- 	varimp$lmg[3]
  }
}


############ Dependence of relative sensitivites on elevation (Figure 3) #############

xtick = c(0,31,59,90,120,151,181,212,243,273,304,334,365)
xtlab = c("J","F","M","A","M","J","J","A","S","O","N","D",NA)
xx <- seq(1,52)*7-3
rb <- (hsv(h = c(seq(.02,.02,length.out=10), 0, seq(.6,.6,length.out=10)), v = c(seq(1,1,length.out=10), 1,seq(1,1,length.out=10)), s=c(seq(1,0.1,length.out=10), 0, seq(0.1,1,length.out=10))))  #blue white red
yg <- (hsv(h = seq(.175,.6,length.out=20), v = seq(1,1,length.out=20), s=seq(0.2,1,length.out=20)))  # yellow green blue
yl <- (hsv(h = c(seq(.175,0,length.out=10),seq(1,1-0.175,length.out=10)), v = seq(1,1,length.out=20), s=c(seq(0,1,length.out=10),seq(1,1,length.out=10))))  # yellow to lila
yg <- (hsv(h = c(seq(.16,.17,length.out=5),seq(.17,.75,length.out=15)), v = seq(1,1,length.out=20), s=c(seq(0.0,1,length.out=5),seq(1,1,length.out=15))))  # yellow green blue
bv	<- (hsv(h = c(seq(.6,.6,length.out=10), seq(.6,.8,length.out=10)), v = c(seq(1,1,length.out=10), seq(1,1,length.out=10)), s=c(seq(0.0,.8,length.out=10), seq(0.8,1,length.out=10)))) 					#white over blue to lila
wbb <- rev(rgb(seq(0,1,length.out=11), seq(0,1,length.out=11), seq(0,1,length.out=11)))			# white to black
wb <- rev(rgb(seq(0,1,length.out=11), seq(0,1,length.out=11), seq(1,1,length.out=11)))  # white to dark blue
wbt <- rev(rgb(seq(1,1,length.out=11), seq(1,1,length.out=11), seq(1,1,length.out=11), seq(0.1,0.9,length.out=11)))			# white to black but variable transparent

elim1 <- c(500,3210)
elim2 <- c(0,64)

abs <- 2

pdf(file="D:/R/weekly_streamflow_sensitivities/Figure 3 new.pdf",width=8, height=8)

#par( mfrow = c(2,2),mar=c(c(2, 4, 1, 1) + 0.1), cex=0.8)
#par(mfrow = c(2, 2), mar = c(2, 4, 1, 1) + 0.1, oma = c(0, 0, 0, 3), cex = 0.8)
layout_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4, ncol = 2, byrow = TRUE)
layout(layout_matrix, heights = c(3, 1.35, 3, 1.35))
#par(mar = c(1, 4, 1, 1) + 0.1, oma = c(0, 0, 0, 3), cex = 0.8)
par(oma = c(0, 0, 0, 3), cex = 0.8)

if(abs == 1){ maxt <- 1.5
maxp <- 1.0
title_legend_T <- expression(b[T]*" [(mm/d)/K]")
title_legend_P <- expression(b[P]*" [-]")}

if(abs == 2){ maxt <- 25
maxp <- 30
title_legend_T <- expression(S[T]*" [%/K]")
title_legend_P <- expression(S[P]*" [%/(mm/d)]")}

palette(rb)
img <- T_elev
temp <- which(img > maxt)
if(length(temp) > 0){img[temp] <- maxt}
par(mar = c(0, 4, 1, 1) + 0.1)
image(xx, elev, img, zlim=c(-maxt,maxt), ylim=elim1, xlim=c(0,365), col=rb, xlab="", ylab="Elevation [m a.s.l.]", xaxt = "n" )  
legend(1,elim1[2]*0.95, seq(-maxt,maxt,length.out=11),col=rb, fill=seq(1,21,length.out=11), horiz=FALSE, title=title_legend_T, bty="n", xjust=0, cex=0.7)
#sig <- matrix(rep((mod_T$elev_p < 0.05),length(elev)),nrow=52,ncol=length(elev))
#image(xx, elev, sig, col=c(rgb(.3,.3,.3,.3),rgb(1,1,1,0)), add=T, zlim=c(0,1))
sigmod <- which(mod_T$elev_p < 0.05)
points(xx[sigmod],zero[sigmod]+elim1[2]*0.98, pch=15, col="gray")
#filled.contour(xx, elev, img,levels=seq(-maxt,maxt,length.out=21),ylim=elim1, xlim=c(0,365),col=rb, xlab="", ylab="Elevation [m a.s.l.]", xaxt = "n")
abline(elim1[2]*0.96,0)
lines(xx, yy_iso, col="grey50", lwd=2, lty=5)
lines(xx, yy_iso2, col="grey50", lwd=2, lty=2)

palette(wb)
img <- P_elev
temp <- which(img > maxp)
if(length(temp) > 0){img[temp] <- maxp}
#image(xx, elev, img, zlim=c(0,maxp), ylim=elim1, xlim=c(0,365),col=wb, xlab="", ylab="Elevation [m a.s.l.]", xaxt = "n" )
par(mar = c(0, 4, 1, 1) + 0.1)
image(xx, elev, img, zlim=c(0, maxp), ylim=elim1, xlim=c(0, 365), col=wb, xlab="", ylab="", xaxt="n", yaxt="n")
legend(1,elim1[2]*0.95, seq(0,maxp,length.out=6),col=wb, fill=seq(1,11,length.out=6), horiz=FALSE, title=title_legend_P, bty="n", xjust=0, cex=0.7)
#sig <- matrix(rep((mod_P$elev_p < 0.05),length(elev)),nrow=52,ncol=length(elev))
#image(xx, elev, sig, col=c(rgb(.3,.3,.3,.3),rgb(1,1,1,0)), add=T, zlim=c(0,1))
sigmod <- which(mod_P$elev_p < 0.05)
points(xx[sigmod],zero[sigmod]+elim1[2]*0.98, pch=15, col="gray")
abline(elim1[2]*0.96,0)
lines(xx, yy_iso, col="grey50", lwd=2, lty=5)
lines(xx, yy_iso2, col="grey50", lwd=2, lty=2)

# Define custom y-tick positions and labels
y_ticks <- c(3339, 2944, 2549, 2154, 1759, 1364, 969, 574)  # Y-tick positions
#y_labels <- c(8.3, 5.8, 3.3, 0.8,-1.7,-4.2,-6.7)  # Y-axis labels #<- 1:7
y_labels <- c(-6, -4, -2, 0, 2, 4, 6, 8)  # Y-axis labels #<- 1:7
axis(side = 4, at = y_ticks, labels = y_labels, col.axis = "black", las = 0)  #
mtext("Temperature [°C]", side = 4, line = 3, col = "black", cex=0.8)  #

par(mar = c(3, 4, 1, 1) + 0.1)
par(xaxs = "i", yaxs = "i")
plot(xx,mod_T$elev_c*-1000/5, col="grey40", lwd=2, type="l",
     ylim=c(-3,2), xlim=c(0,365), xlab="", ylab=expression(S[HT]*" [*/K]"), xaxt = "n",
     panel.first = abline(0, 0, col = "grey50", lwd = 1, lty = 2))
# Add x-axis
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))

par(mar = c(3, 4, 1, 1) + 0.1)
plot(xx, mod_P$elev_c*-1000/5, col="grey40", lwd=2, type="l",
     ylim=c(-0,4), xlim=c(0,365), xlab="", ylab="", xaxt="n", yaxt="n",
     panel.first = abline(0, 0, col = "grey50", lwd = 1, lty = 2))
# Add x-axis
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
# Add y-axis on the right side
axis(side = 4)
# Add y-axis label on the right side
mtext(expression(S[HP]*" [*/K]"), side = 4, line = 3, cex = par("cex"))

palette(rb)
img <- T_glac
temp <- which(img > maxt)
if(length(temp) > 0){img[temp] <- maxt}
par(mar = c(0, 4, 1, 1) + 0.1)
image(xx, glaciers, img, zlim=c(-maxt,maxt), ylim=elim2, xlim=c(0,365), col=rb, xlab="", ylab="Glacier Cover [%]", xaxt = "n" ) 
#sig <- matrix(rep((mod_T$glac_p < 0.05),length(glaciers)),nrow=52,ncol=length(glaciers))
#image(xx, glaciers, sig, col=c(rgb(.3,.3,.3,.3),rgb(1,1,1,0)), add=T, zlim=c(0,1))
sigmod <- which(mod_T$glac_p < 0.05)
points(xx[sigmod],zero[sigmod]+elim2[2]*0.98, pch=15, col="gray")
abline(elim2[2]*0.96,0)

palette(wb)
img <- P_glac
temp <- which(img > maxp)
if(length(temp) > 0){img[temp] <- maxp}
#image(xx, glaciers, img, zlim=c(0,maxp), ylim=elim2, xlim=c(0,365),col=wb, xlab="", ylab="Glacier Cover [%]", xaxt = "n" )
par(mar = c(0, 4, 1, 1) + 0.1)
image(xx, glaciers, img, zlim=c(0,maxp), ylim=elim2, xlim=c(0,365),col=wb, xlab="", ylab="", xaxt = "n" , yaxt="n")
axis(4)
mtext("Glacier Cover [%]", side = 4, line = 3, col = "black", cex=0.8)  #
#sig <- matrix(rep((mod_P$glac_p < 0.05),length(glaciers)),nrow=52,ncol=length(glaciers))
#image(xx, glaciers, sig, col=c(rgb(.3,.3,.3,.3),rgb(1,1,1,0)), add=T, zlim=c(0,1))
sigmod <- which(mod_P$glac_p < 0.05)
points(xx[sigmod],zero[sigmod]+elim2[2]*0.98, pch=15, col="gray")
abline(elim2[2]*0.96,0)

par(xaxs = "i", yaxs = "i")
par(mar = c(3, 4, 1, 1) + 0.1)
plot(xx,mod_T$glac_c, col="grey40", lwd=2, type="l",
     ylim=c(-0.2,0.3), xlim=c(0,365), xlab="", ylab=expression(S[GT]*" [*/%]"), xaxt = "n",
     panel.first = abline(0, 0, col = "grey50", lwd = 1, lty = 2)) 
# Add x-axis
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))

par(mar = c(3, 4, 1, 1) + 0.1)
plot(xx, mod_P$glac_c, col="grey40", lwd=2, type="l",
     ylim=c(-0.1,0.3), xlim=c(0,365), xlab="", ylab="", xaxt="n", yaxt="n",
     panel.first = abline(0, 0, col = "grey50", lwd = 1, lty = 2))
# Add x-axis
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
# Add y-axis on the right side
axis(side = 4)
# Add y-axis label on the right side
mtext(expression(S[GP]*" [*/%]"), side = 4, line = 3, cex = par("cex"))

dev.off() 





