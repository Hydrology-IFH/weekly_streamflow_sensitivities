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


############ Supporting information #############s

#col <- (hsv(h = seq(.175,.6,length.out=nr), v = seq(1,1,length.out=nr), s=seq(0.2,1,length.out=nr))) 
library("RColorBrewer")
#col <- colorRampPalette(brewer.pal(n=9, name='YlGnBu'))(30)
col <- terrain.colors(30, alpha = 1)
palette(col)

xx <- seq(1,52)*7-3
xtick = c(0,31,59,90,120,151,181,212,243,273,304,334,365)
xtlab = c("J","F","M","A","M","J","J","A","S","O","N","D",NA)
zero <- rep(0, length(xx))

# Regime curves (Figure S1)
pdf(file="D:/R/weekly_streamflow_sensitivities/Figure S1.pdf",width=6, height=4)

par(mfrow = c(1,1), 
    mar = c(3, 4, 1, 1) + 0.1, 
    cex = 0.9)

ylim1 <- c(0.1,100)
plot(xx,out_Q.p[1,], typ="l", ylim=ylim1, col=col[(height[1]-500)*.01], xlab="", ylab="Average streamflow [mm/d]", xaxt = "n", lwd = 2 ,  xaxs="i", yaxs="i", log = "y")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
#abline(0,0)
for (i in 1:nr) {points(xx,out_Q.p[i,], typ="l",col=col[ceiling((height[i]-500)*.01)], lwd = 2)}
legend(1,ylim1[2], seq(500,3000,length.out=11),col=col, fill=seq(1,25,length.out=11), horiz=F, title="H [m a.s.l.]", bty="n", xjust=-0.1, cex=0.6)

dev.off() 


# Total variation explained (Figure S2)

### explained variance
for (i in 1:nr) {tcoef[i,] <- ma.filtc(as.numeric(out_T.var[i,]),movavg)}  
for (i in 1:nr) {pcoef[i,] <- ma.filtc(as.numeric(out_P.var[i,]),movavg)}
for (i in 1:nr) {qcoef[i,] <- ma.filtc(as.numeric(out_Q.var[i,]),movavg)}

clim_coef = tcoef + pcoef + qcoef

pdf(file="D:/R/weekly_streamflow_sensitivities/Figure S2.pdf",width=6, height=4)

par(mfrow = c(1,1), 
    mar = c(3, 4, 1, 1) + 0.1, 
    cex = 0.9)

ylim1 <- c(0,1)
plot(xx,clim_coef[1,], typ="l", ylim=ylim1, col=col[(height[1]-500)*.01], xlab="", ylab="Total variation explained (R²) [-]", xaxt = "n", lwd = 2 ,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)
for (i in 1:nr) {points(xx,clim_coef[i,], typ="l",col=col[ceiling((height[i]-500)*.01)], lwd = 2)}
legend(x = "bottomleft", legend = seq(500, 3000, length.out = 11), col = col, fill = seq(1, 25, length.out = 11), horiz = FALSE, title = "H [m a.s.l.]", bty = "n", cex = 0.6)

dev.off() 


# example catchments (Figure S3)

pdf(file="D:/R/weekly_streamflow_sensitivities/Figure S3.pdf",width=8, height=10)

#par( mfrow = c(3,3),mar=c(c(3, 4, 1, 1) + 0.1), cex=0.9)
#par(xaxs = "i", yaxs = "i")
par(mfrow = c(4,3),
    mar = c(2.5, 3, 1.5, 0.5),  # Increased top margin slightly
    oma = c(0.5, 0.5, 0.5, 0),  # Added small outer margin at top
    mgp = c(2, 0.7, 0),         # Axis label positioning
    cex = 0.9,
    xaxs = "i", 
    yaxs = "i")

# selected catchments
n1 = 25
n2 = 10
n3 = 19

### streamflow
ylim1 <- c(-0.8,2.6)
plot(xx,out_Q.p[n1,], typ="l", col=col[(height[n1]-500)*.01], xlab="", ylab = expression("Q [mm/d]"), xaxt = "n", lwd = 2,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)

ylim1 <- c(-0.1,1.0)
plot(xx,out_Q.p[n2,], typ="l", col=col[(height[n2]-500)*.01], xlab="", ylab = expression("Q [mm/d]"), xaxt = "n", lwd = 2 ,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)

ylim1 <- c(-0.1,1.4)
plot(xx,out_Q.p[n3,], typ="l", col=col[(height[n3]-500)*.01], xlab="", ylab = expression("Q [mm/d]"), xaxt = "n", lwd = 2,  xaxs="i" )
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)

### temperature sensitivity
ylim1 <- c(-20,25)
plot(xx,out_T.coeff[n1,], typ="l", col=col[(height[n1]-500)*.01], xlab="", ylab = expression(b[T]*" [(mm/d)/K]"), xaxt = "n", lwd = 2,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)
#legend(183,ylim1[2], seq(500,3000,length.out=6),col=col, fill=seq(1,25,length.out=6), horiz=TRUE, title="Elevation [m a.s.l.]", bty="n", xjust=0.5, cex=0.8)

ylim1 <- c(-5,50)
plot(xx,out_T.coeff[n2,], typ="l", col=col[(height[n2]-500)*.01], xlab="", ylab = expression(b[T]*" [(mm/d)/K]"), xaxt = "n", lwd = 2 ,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)

ylim1 <- c(-5,50)
plot(xx,out_T.coeff[n3,], typ="l", col=col[(height[n3]-500)*.01], xlab="", ylab = expression(b[T]*" [(mm/d)/K]"), xaxt = "n", lwd = 2,  xaxs="i" )
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)

### precipitation sensitivity
ylim1 <- c(0,1)
plot(xx,out_P.coeff[n1,], typ="l", col=col[(height[n1]-500)*.01], xlab="", ylab = expression(b[P]*" [-]"), xaxt = "n", lwd = 2 ,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)

ylim1 <- c(0,1)
plot(xx,out_P.coeff[n2,], typ="l", col=col[(height[n2]-500)*.01], xlab="", ylab = expression(b[P]*" [-]"), xaxt = "n", lwd = 2 ,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)

ylim1 <- c(0,1)
plot(xx,out_P.coeff[n3,], typ="l", col=col[(height[n3]-500)*.01], xlab="", ylab = expression(b[P]*" [-]"), xaxt = "n", lwd = 2,  xaxs="i" )
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)

### storage sensitivity
ylim1 <- c(0,1)
plot(xx,out_Q.coeff[n1,], typ="l", col=col[(height[n1]-500)*.01], xlab="", ylab = expression(b[Q]*" [-]"), xaxt = "n", lwd = 2 ,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)

ylim1 <- c(0,1)
plot(xx,out_Q.coeff[n2,], typ="l", col=col[(height[n2]-500)*.01], xlab="", ylab = expression(b[Q]*" [-]"), xaxt = "n", lwd = 2 ,  xaxs="i", yaxs="i")
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)

ylim1 <- c(0,1)
plot(xx,out_Q.coeff[n3,], typ="l", col=col[(height[n3]-500)*.01], xlab="", ylab = expression(b[Q]*" [-]"), xaxt = "n", lwd = 2,  xaxs="i" )
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
abline(0,0)

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


# Apply linear models to look for elevation dependence

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


####################################
# make plots combining the infos for the regression with elevation and glacier coverage (Figure S4)

xtick = c(0,31,59,90,120,151,181,212,243,273,304,334,365)
xtlab = c("J","F","M","A","M","J","J","A","S","O","N","D",NA)
zero <- rep(0, length(xx))

pdf(file="D:/R/weekly_streamflow_sensitivities/Figure S4ab.pdf",width=8, height=4)

par( mfrow = c(1,2),mar=c(c(4, 4, 2, 1) + 0.1), cex=0.8, xaxs="i", yaxs="i")

plot(xx,mod_T$r2, ylim=c(0,1.15), typ="l", lwd=2, xlab="", ylab="Explained Variation [-]", xaxt = "n" )
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
polygon(c(xx,rev(xx)),c(mod_T$elev_r,zero), col="orange", border=NA)	
polygon(c(xx,rev(xx)),c((mod_T$elev_r+mod_T$glac_r),rev(mod_T$elev_r)), col="cyan", border=NA)	

sigmod <- which(mod_T$p < 0.05)
points(xx[sigmod],zero[sigmod]+1, pch=16, col="gray")
sigmod <- which(mod_T$elev_p < 0.05 & mod_T$elev_c > 0.0)
points(xx[sigmod],zero[sigmod]+1.04, pch="+", col="orange")
sigmod <- which(mod_T$elev_p < 0.05 & mod_T$elev_c < 0.0)
points(xx[sigmod],zero[sigmod]+1.04, pch="-", col="orange")
sigmod <- which(mod_T$glac_p < 0.05 & mod_T$glac_c > 0.0)
points(xx[sigmod],zero[sigmod]+1.08, pch="+", col="cyan")
sigmod <- which(mod_T$glac_p < 0.05 & mod_T$glac_c < 0.0)
points(xx[sigmod],zero[sigmod]+1.08, pch="-", col="cyan")

if(inter == 1){
  polygon(c(xx,rev(xx)),c(mod_T$r2,rev(mod_T$elev_r+mod_T$glac_r)), col="blue", border=NA)
  sigmod <- which(mod_T$inter_p < 0.05 & mod_T$inter_c > 0.0)
  points(xx[sigmod],zero[sigmod]+1.12, pch="+", col="blue")
  sigmod <- which(mod_T$inter_p < 0.05 & mod_T$inter_c < 0.0)
  points(xx[sigmod],zero[sigmod]+1.12, pch="-", col="blue")}

plot(xx,mod_P$r2, ylim=c(0,1.15), typ="l", lwd=2, xlab="", ylab="Explained Variation [-]", xaxt = "n" )
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
polygon(c(xx,rev(xx)),c(mod_P$elev_r,zero), col="orange", border=NA)	
polygon(c(xx,rev(xx)),c((mod_P$elev_r+mod_P$glac_r),rev(mod_P$elev_r)), col="cyan", border=NA)	

sigmod <- which(mod_P$p < 0.05)
points(xx[sigmod],zero[sigmod]+1, pch=16, col="gray")
sigmod <- which(mod_P$elev_p < 0.05 & mod_P$elev_c > 0.0)
points(xx[sigmod],zero[sigmod]+1.04, pch="+", col="orange")
sigmod <- which(mod_P$elev_p < 0.05 & mod_P$elev_c < 0.0)
points(xx[sigmod],zero[sigmod]+1.04, pch="-", col="orange")
sigmod <- which(mod_P$glac_p < 0.05 & mod_P$glac_c > 0.0)
points(xx[sigmod],zero[sigmod]+1.08, pch="+", col="cyan")
sigmod <- which(mod_P$glac_p < 0.05 & mod_P$glac_c < 0.0)
points(xx[sigmod],zero[sigmod]+1.08, pch="-", col="cyan")

if(inter == 1){
  polygon(c(xx,rev(xx)),c(mod_P$r2,rev(mod_P$elev_r+mod_P$glac_r)), col="blue", border=NA)
  sigmod <- which(mod_P$inter_p < 0.05 & mod_P$inter_c > 0.0)
  points(xx[sigmod],zero[sigmod]+1.12, pch="+", col="blue")
  sigmod <- which(mod_P$inter_p < 0.05 & mod_P$inter_c < 0.0)
  points(xx[sigmod],zero[sigmod]+1.12, pch="-", col="blue")}

dev.off()


### REDO WITH RELATIVE SENSTIVITIES

# Apply linear models to look for elevation dependence

abs <- 2 #FALSE #TRUE			# if abs 1 than plot absolute sensitivity, if 2 then relative, if 3 than elasticity

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


####################################
# make plots (Figure S2) combining the infos for the regression with elevation and glacier coverage (Figure S4)

xtick = c(0,31,59,90,120,151,181,212,243,273,304,334,365)
xtlab = c("J","F","M","A","M","J","J","A","S","O","N","D",NA)
zero <- rep(0, length(xx))

pdf(file="D:/R/weekly_streamflow_sensitivities/Figure S4cd.pdf",width=8, height=4)

par( mfrow = c(1,2),mar=c(c(4, 4, 2, 1) + 0.1), cex=0.8, xaxs="i", yaxs="i")

plot(xx,mod_T$r2, ylim=c(0,1.15), typ="l", lwd=2, xlab="", ylab="Explained Variation [-]", xaxt = "n" )
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
polygon(c(xx,rev(xx)),c(mod_T$elev_r,zero), col="orange", border=NA)	
polygon(c(xx,rev(xx)),c((mod_T$elev_r+mod_T$glac_r),rev(mod_T$elev_r)), col="cyan", border=NA)	

sigmod <- which(mod_T$p < 0.05)
points(xx[sigmod],zero[sigmod]+1, pch=16, col="gray")
sigmod <- which(mod_T$elev_p < 0.05 & mod_T$elev_c > 0.0)
points(xx[sigmod],zero[sigmod]+1.04, pch="+", col="orange")
sigmod <- which(mod_T$elev_p < 0.05 & mod_T$elev_c < 0.0)
points(xx[sigmod],zero[sigmod]+1.04, pch="-", col="orange")
sigmod <- which(mod_T$glac_p < 0.05 & mod_T$glac_c > 0.0)
points(xx[sigmod],zero[sigmod]+1.08, pch="+", col="cyan")
sigmod <- which(mod_T$glac_p < 0.05 & mod_T$glac_c < 0.0)
points(xx[sigmod],zero[sigmod]+1.08, pch="-", col="cyan")

if(inter == 1){
  polygon(c(xx,rev(xx)),c(mod_T$r2,rev(mod_T$elev_r+mod_T$glac_r)), col="blue", border=NA)
  sigmod <- which(mod_T$inter_p < 0.05 & mod_T$inter_c > 0.0)
  points(xx[sigmod],zero[sigmod]+1.12, pch="+", col="blue")
  sigmod <- which(mod_T$inter_p < 0.05 & mod_T$inter_c < 0.0)
  points(xx[sigmod],zero[sigmod]+1.12, pch="-", col="blue")}

plot(xx,mod_P$r2, ylim=c(0,1.15), typ="l", lwd=2, xlab="", ylab="Explained Variation [-]", xaxt = "n" )
axis(side = 1, at = xtick+13, labels = xtlab, lwd.ticks=0)
axis(side = 1, at = xtick, labels = rep("",13))
polygon(c(xx,rev(xx)),c(mod_P$elev_r,zero), col="orange", border=NA)	
polygon(c(xx,rev(xx)),c((mod_P$elev_r+mod_P$glac_r),rev(mod_P$elev_r)), col="cyan", border=NA)	

sigmod <- which(mod_P$p < 0.05)
points(xx[sigmod],zero[sigmod]+1, pch=16, col="gray")
sigmod <- which(mod_P$elev_p < 0.05 & mod_P$elev_c > 0.0)
points(xx[sigmod],zero[sigmod]+1.04, pch="+", col="orange")
sigmod <- which(mod_P$elev_p < 0.05 & mod_P$elev_c < 0.0)
points(xx[sigmod],zero[sigmod]+1.04, pch="-", col="orange")
sigmod <- which(mod_P$glac_p < 0.05 & mod_P$glac_c > 0.0)
points(xx[sigmod],zero[sigmod]+1.08, pch="+", col="cyan")
sigmod <- which(mod_P$glac_p < 0.05 & mod_P$glac_c < 0.0)
points(xx[sigmod],zero[sigmod]+1.08, pch="-", col="cyan")

if(inter == 1){
  polygon(c(xx,rev(xx)),c(mod_P$r2,rev(mod_P$elev_r+mod_P$glac_r)), col="blue", border=NA)
  sigmod <- which(mod_P$inter_p < 0.05 & mod_P$inter_c > 0.0)
  points(xx[sigmod],zero[sigmod]+1.12, pch="+", col="blue")
  sigmod <- which(mod_P$inter_p < 0.05 & mod_P$inter_c < 0.0)
  points(xx[sigmod],zero[sigmod]+1.12, pch="-", col="blue")}

dev.off()





