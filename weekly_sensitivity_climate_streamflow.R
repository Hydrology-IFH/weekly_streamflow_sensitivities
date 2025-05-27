#NA functions

mean.na <- function(x,...){if (all(is.na(x))){NA} else {mean(x, na.rm = TRUE)}}
sum.na <- function(x,...){if (all(is.na(x))){NA} else {sum(x, na.rm = TRUE)}}
min.na <- function(x,...){if (all(is.na(x))){NA} else {min(x, na.rm = TRUE)}}
max.na <- function(x,...){if (all(is.na(x))){NA} else {max(x, na.rm = TRUE)}}
which.max.na <- function(x,...){if (all(is.na(x))){NA} else {which.max(x)}}
which.min.na <- function(x,...){if (all(is.na(x))){NA} else {which.min(x)}}

quantile.na <- function(x, probs = seq(0, 1, 0.25)){
  if (all(is.na(x))){rep(NA, length(probs))} else {quantile(x, probs = probs, 
                                                            type = 1, na.rm = TRUE)}}

#concatenate
con <- function(...){paste(..., sep = "")}

#quantile function
probs <- c(0.0, seq(0.1, 0.9, by = 0.1), 1.00)
probs <- seq(0.1, 0.9, 0.4)
qfun <- function(x){quantile(x, probs = probs, type = 1, na.rm = TRUE)}



#moving average; repetition at begin and end
ma.filt <- function(x, window)
{
  m <- floor(window/2)
  y <- filter(c(rep(x[1], m), x, rep(x[length(x)], m)),
              rep(1/window, window), sides = 2, method = "convolution")
  return(y[(m+1):(m+length(x))])
}

#moving average; circular
ma.filtc <- function(x, window)
{
  m <- floor(window/2)
  y <- filter(x, rep(1/window, window), sides = 2, method = "convolution", circular = TRUE)
  return(y[1:length(x)])
}



#concatenate
con <- function(...){paste(..., sep = "")}

#omit Feb 29
omit.feb29 <- function(x){x <- x[x$date != "02-29",]
x$date <- as.Date(con("2011-", x$date))
x[order(x$date),]}

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
              slope = slope, intercept=intercept, median=median(x)))
}



#####################################
#Snow and degree day function
cal.SWE <- function(prec, temp, cfmax = 4.5, tt = 0 )
{       	
  #input
  #prec: precipitation [mm/d]
  #temp: temperature [°C]
  #cfmax: degree-day factor [mm/(°C d)]
  #tt: threshold temperature  [°C]
  
  NAvalues 	   	<- is.na(temp) 
  temp[NAvalues]		<- 5
  
  NAvalues1 	   	<- is.na(prec) 
  prec[NAvalues1]		<- 0
  
  #dataframe
  a <- 0
  a <- data.frame(prec = prec, temp = temp)
  #time series length
  n <- nrow(a)
  
  #time independent procedures
  #if temperature <= tt: snow accumulation
  ob.cold <- a$temp <= tt
  a$accu <- 0
  a$accu[ob.cold] <- a$prec[ob.cold]
  
  #if temperature > tt: snow melt
  #potential snow melt based on day-degree method
  a$psm <- 0
  a$psm[!ob.cold] <- cfmax * (a$temp[!ob.cold] - tt)
  
  #if temperature > tt: rain
  a$rain <- 0
  a$rain[!ob.cold] <- a$prec[!ob.cold]
  
  #object initialization
  a$sp <- 0     #snow pack
  a$out <- 0    #melt water
  
  
  #loop over time
  for (i in 2:n)
  {
    #i <- 5
    
    #snow accumumation
    a$sp[i] <- a$sp[i-1] + a$accu[i]
    
    #snow melt
    a$out[i] <- min(a$sp[i], a$psm[i])
    a$sp[i] <- a$sp[i] - a$out[i]
    
    #end of loop over time
  }
  
  a[(NAvalues | NAvalues1),]		<- NA				# fill vlaues with NA that had no input
  
  b <- 0
  b <- data.frame(cbind(a$out, a$sp))
  names(b) <- c("out", "swe")
  return(b)
  
  #output
  #out: output [mm/d]
  #swe: snow pack [mm]
}



#______________________________________________________________________________
# Main Programm



library(relaimpo)   ## zur Berechnung der Variable Importance
# Groemping, U. (2006). Relative Importance for Linear Regression in R: The Package relaimpo. Journal of Statistical Software, 17(1), 1–27.
# Groemping, U. (2007). Estimators of Relative Importance in Linear Regression Based on Variance Decomposition. The American Statistician, 61, 139-147. doi:10.1198/000313007X188252
library(car)
library(dplyr)
library(lubridate)




path_p <- "D:/R/weekly_streamflow_sensitivities/data/"	
path_t <- "D:/R/weekly_streamflow_sensitivities/data_temp/"

path_q <- "D:/R/weekly_streamflow_sensitivities/BAFUvergletscherte/"

path <- "D:/R/weekly_streamflow_sensitivities/Abflussdaten/"

#stations and pairs

files <- read.csv(paste(path,"Klima_Abfluss_Vergleich.csv", sep =""), header=TRUE)

nr <- length(files$Q_Nr)

files$start <- rep(NA,nr)
files$end <- rep(NA,nr)

pdf(file=paste(path,"Weekly Sensitivity.pdf", sep ="" ), width=16, height=75)

# figure (matrix plot)
rb <- hsv(h = c(seq(.02,.02,length.out=5),0, seq(.6,.6,length.out=5)), v = c(seq(1,1,length.out=5),1, seq(1,1,length.out=5)), s=c(seq(1,0,length.out=5),0, seq(0,1,length.out=5)))

#layout(matrix(c(seq(from=1,by=1,length.out=25)),5,5,byrow=TRUE))  # anordnung 5*5

layout(matrix(c((seq(from=1,by=5,length.out=25)),(seq(from=2,by=5,length.out=25)),(seq(from=3,by=5,length.out=25)),(seq(from=4,by=5,length.out=25)),(seq(from=5,by=5,length.out=25))),25,5,byrow=TRUE))  

par(mar=c(c(4, 4, 3, 1) + 0.1), cex=0.7)   #(3, 3, 4, 2)

# define data frame for final results for annual data

precip <- data.frame(year=seq(1925,2012))
temp <- data.frame(year=seq(1925,2012))
discharge <- data.frame(year=seq(1925,2012))


# define data table for final results of climate sensitivity
out_T <- files										# Temperature sensitivity
out_T$p <- matrix(NA, nrow=nr, ncol=52)				# mean predictor
out_T$coeff <- matrix(NA, nrow=nr, ncol=52) 		# coefficients
out_T$var_imp <- matrix(NA, nrow=nr, ncol=52) 		# Variable importance

out_P <- files										# Precipitation sensitivity
out_P$p <- matrix(NA, nrow=nr, ncol=52)				# mean predictor
out_P$coeff <- matrix(NA, nrow=nr, ncol=52) 		# coefficients
out_P$var_imp <- matrix(NA, nrow=nr, ncol=52) 		# Variable importance

out_Q <- files										# Discharge lag1 sensitivity
out_Q$p <- matrix(NA, nrow=nr, ncol=52)				# mean predictor
out_Q$coeff <- matrix(NA, nrow=nr, ncol=52) 		# coefficients
out_Q$var_imp <- matrix(NA, nrow=nr, ncol=52) 		# Variable importance

out_SWE <- files										# SWE
out_SWE$week <- matrix(NA, nrow=nr, ncol=52) 		# Variable importance
out_SWE$SWE <- matrix(NA, nrow=nr, ncol=52) 		# Variable importance

out_BFI <- files										# BFI
out_BFI$BFI <- matrix(NA, nrow=nr, ncol=1) 	

#  Main Loop

for (ii in 1:nr)   # nr
{
  
  
  # # Berechnung für ganzes Jahr (1-12) oder nur für Sommer (4-9)
  m_start  <-  1
  m_end   <- 12
  
  
  
  # read data
  
  a <- read.table(paste(path,"Data_Q_",files$Q_Nr[ii],"_C_",files$C_Nr[ii],".txt", sep ="" ))
  
  files$start[ii] <- min(a$y)
  files$end[ii] <- max(a$y)
  
  
  # calculate BFI
  library(UKFE)
  bfi_result <- BFI(Q =  a$Q, Plot = FALSE)
  out_BFI$BFI[ii,] <- bfi_result
  
  #annual time series: average or sum
  
  # select only summer
  wo1 <- which((as.numeric(a$m) >= m_start)  & (as.numeric(a$m) <= m_end)) 
  
  
  func <- function(x,...){
    method <- list(...)
    lapply(method, function(m)do.call(
      function(z)aggregate(z, list(year = a$y[wo1]), m,
                           na.rm = TRUE), list(x)))}
  
  
  annu <- func(a$Q[wo1], sum.na)
  annu1 <- func(a$swe[wo1], max.na)
  annu2 <- func(a$P_obs[wo1], sum.na)
  annu3 <- func(a$T_obs[wo1], mean.na)
  annu4 <- func(a$swe1[wo1], max.na)
  
  annu5 <- func(a$Q[wo1]>-99, sum.na)
  
  annu <- as.data.frame(annu)
  annu1 <- as.data.frame(annu1)
  annu2 <- as.data.frame(annu2)
  annu3 <- as.data.frame(annu3)
  annu4 <- as.data.frame(annu4)	
  annu5 <- as.data.frame(annu5)
  
  annu$x[annu5$x < ((m_end-m_start+1)*28)] <- NA		# make sure enough data per year available to calulate mean or sum
  annu1$x[annu5$x < ((m_end-m_start+1)*28)] <- NA
  annu2$x[annu5$x < ((m_end-m_start+1)*28)] <- NA
  annu3$x[annu5$x < ((m_end-m_start+1)*28)] <- NA
  annu4$x[annu5$x < ((m_end-m_start+1)*28)] <- NA
  
  annu <- merge(annu, annu1, by ="year")
  annu <- merge(annu, annu2, by ="year")
  names(annu)[2:4] <- c("Q", "swe","P")	
  annu <- merge(annu, annu3, by ="year")
  annu <- merge(annu, annu4, by ="year")
  names(annu)[5:6] <- c("T", "swe1")
  annu$year <- as.numeric(annu$year)
  head(annu)
  
  
  # print data into overall table
  
  precip <- merge(precip, annu[,c(1,4)], by="year", all.x=TRUE)
  names(precip)[ii+1] <- files$Q_Nr[ii]
  temp <- merge(temp, annu[,c(1,5)], by="year", all.x=TRUE)
  names(temp)[ii+1] <- files$Q_Nr[ii]
  discharge <- merge(discharge, annu[,c(1,2)], by="year", all.x=TRUE)
  names(discharge)[ii+1] <- files$Q_Nr[ii]
  
  
  # monatliche Mittel, summe oder maxima
  am<-0
  am <- aggregate(a$T_obs, list(ym = a$ym), mean.na)
  am$date <- as.Date(con(am$ym, "-1"))
  am$m <- format(am$date, "%m")
  
  a1 <- aggregate(a$P_obs, list(ym = a$ym), sum.na)
  am <- merge(am, a1, by ="ym")
  
  names(am) <- c("ym", "T", "date", "m","P")
  
  a1 <- aggregate(a$swe, list(ym = a$ym), max.na)
  am <- merge(am, a1, by ="ym")
  
  a1 <- aggregate(a$Q, list(ym = a$ym), sum.na)
  am <- merge(am, a1, by ="ym")
  
  names(am) <- c("ym", "T", "date", "m","P","swe","Q")
  
  a1 <- aggregate(a$inp, list(ym = a$ym), sum.na)
  am <- merge(am, a1, by ="ym")
  
  names(am) <- c("ym", "T", "date", "m","P","swe","Q", "inp")
  
  am$y <- format(am$date, "%Y")
  
  
  
  
  # weekly Mittel, summe oder maxima
  aw<-0
  aw <- aggregate(a$T_obs, list(yw = a$yw), mean.na)
  aw$y <- as.character(floor(aw$yw))
  aw$date <- as.Date(paste(aw$y,as.character(round(((aw$yw-floor(aw$yw))*52+1)*7)),sep="-"), "%Y-%j")
  aw$w <- round(((aw$yw-floor(aw$yw))*52+1))
  
  a1 <- aggregate(a$P_obs, list(yw = a$yw), mean.na)			# alles in mm/Tag (Mittelwert)
  aw <- merge(aw, a1, by ="yw")
  
  names(aw) <- c("yw", "T", "y", "date", "w","P")
  
  a1 <- aggregate(a$swe, list(yw = a$yw), max.na)
  aw <- merge(aw, a1, by ="yw")
  
  a1 <- aggregate(a$Q, list(yw = a$yw), mean.na)
  aw <- merge(aw, a1, by ="yw")
  
  names(aw) <- c("yw", "T", "y", "date", "w","P","swe","Q")
  
  a1 <- aggregate(a$inp, list(yw = a$yw), mean.na)
  aw <- merge(aw, a1, by ="yw")
  
  names(aw) <- c("yw", "T", "y", "date", "w","P","swe","Q", "inp")
  
  
  aw$m <- format(aw$date, "%m")
  
  #  second snow elevation
  a1 <- aggregate(a$inp1, list(yw = a$yw), sum.na)
  aw <- merge(aw, a1, by ="yw")
  
  names(aw) <- c("yw", "T", "y", "date", "w","P","swe","Q", "inp", "m", "inp1")
  
  aw$Q_1 <- aw$Q
  aw$Q_1[2:length(aw$Q)] <- aw$Q[2:length(aw$Q)-1]   # Q of time step t-1 for autocorrelation
  
  
  # calc average SWE per week
  aw <- aw %>%
    mutate(week = week(date))
    weekly_avg_swe <- aw %>%
    group_by(week) %>%
    summarise(avg_swe = mean(swe, na.rm = TRUE))  # na.rm = TRUE to handle missing values
  #plot(weekly_avg_swe$week, weekly_avg_swe$avg_swe, xlab="Date", ylab="SWE", col="blue", lwd=1, type="l")
  
  
  # analysis with weekly data
  
  ny <- length(annu$y)
  ws <- 52
  
  start <- 1900
  end <- 2014
  
  pred_w <- matrix(NA, nrow=ny, ncol=ws)   #  prediction
  res_w <- matrix(NA, nrow=ny, ncol=ws)    #   residue
  res_w_l <- matrix(NA, nrow=ny, ncol=ws)  # smooth residu
  w_lm <- data.frame(w=seq(1,ws), t_c=rep(NA,ws), i_c=rep(NA, ws), i1_c=rep(NA,ws), r2=rep(NA,ws), t_r=rep(0,ws), i_r=rep(0, ws), i1_r=rep(0,ws),  q=rep(NA,ws), t_l=rep(0,ws), i_l=rep(0, ws), i1_l=rep(0,ws),  t_u=rep(0,ws), i_u=rep(0, ws), i1_u=rep(0,ws), t=rep(NA,ws), p=rep(NA,ws), b0=rep(NA,ws), a0=rep(NA,ws), partial_cor_P=rep(NA,ws))
  w_mod <- list()
  #  c = coefficient, r = variable importance, l = lower Konfidenzbereich, u = upper confidence limit
  
  # Model and Plot for 2 predictors
  
  for (mw in 1:ws) {
    good <- which((as.numeric(aw$w) == mw) & (is.finite(aw$Q))& (is.finite(aw$T)) & (is.finite(aw$P)) & (is.finite(aw$Q_1)) & (aw$y >= start) & (aw$y < end))
    
    if (length(good) > 15) {
      Q <- aw$Q[good]
      T <- aw$T[good]
      inp <- aw$P[good]     # precip 
      #inp <- aw$inp[good]   # snow
      swe <- aw$swe[good]
      inp1 <- aw$Q_1[good]			# autocorrelation
      Y <- aw$y[good]
      
      w_lm$q[mw] <- mean.na(Q)
      w_lm$t[mw] <- mean.na(T)
      w_lm$p[mw] <- mean.na(inp)
      
      #  auch noch Abfluss der Woche davor dazu um die Autokorrelation zu berücksichtigen
      
      wo1 <- match(Y, annu$y)
      
      #mod <- loess(Q ~ T + inp + inp1, degree=1, span=.75)   
      mod <- lm(Q ~ T + inp + inp1)  
      
      w_lm$a0[mw] <- 	mod$coefficients[1]
      w_lm$t_c[mw] <- 	mod$coefficients[2]
      w_lm$i_c[mw] <- 	mod$coefficients[3]
      w_lm$i1_c[mw] <- 	mod$coefficients[4]
      
      # berechnung der Konfidenzbereiche des Coeffizinten:  confint(mod)
      conf <- confint(mod)
      w_lm$t_l[mw] <- 	conf[2,1]
      w_lm$i_l[mw] <- 	conf[3,1]
      w_lm$i1_l[mw] <- 	conf[4,1]
      w_lm$t_u[mw] <- 	conf[2,2]
      w_lm$i_u[mw] <- 	conf[3,2]
      w_lm$i1_u[mw] <- 	conf[4,2]
      
      
      varimp <- calc.relimp(mod)
      w_lm$t_r[mw] <- 	varimp$lmg[1]
      w_lm$i_r[mw] <- 	varimp$lmg[2]
      w_lm$i1_r[mw] <- 	varimp$lmg[3]		
      
      w_mod[[mw]] <- mod
      
      w_lm$r2[mw] <- cor(mod$fitted,Q)^2	
      
      
      pred_w[wo1,mw] <- mod$fitted
      res_w[wo1,mw] <- mod$residuals
      mod1 <- loess((predict(mod)-Q) ~ Y, degree=2, span=.6)
      res_w_l[wo1,mw] <- predict(mod1)
      
      # partial corr test
      library(ppcor)
      data <- data.frame(Q = Q, inp = inp, T = T, inp1 = inp1)
      result <- pcor(data)
      w_lm$partial_cor_P[mw] <- result$estimate["Q", "inp"]
      
      
    }
  }
  
  # partial correlations
  #dev.new()
  #plot(w_lm$w, w_lm$i_r, type = "l", col = "blue", xlab = "w", ylab = "Values", ylim = range(c(w_lm$i_r, w_lm$partial_cor_P)), lwd = 2)
  #lines(w_lm$w, w_lm$partial_cor_P, type = "l", col = "orange", lwd = 2)
  #legend("bottomleft", legend = c("i_r", "partial_cor_P"), col = c("blue", "orange"), lty = 1, lwd = 2)
  
  
  #  write weekly sensitivity data for each station
  
  out_T$p[ii,] <- w_lm$t
  out_T$coeff[ii,] <- w_lm$t_c
  out_T$var_imp[ii,] <- w_lm$t_r
  
  out_P$p[ii,] <- w_lm$p
  out_P$coeff[ii,] <- w_lm$i_c
  out_P$var_imp[ii,] <- w_lm$i_r
  
  out_Q$p[ii,] <- w_lm$q
  out_Q$coeff[ii,] <- w_lm$i1_c
  out_Q$var_imp[ii,] <- w_lm$i1_r
  #out_Q$a0[ii,] <- w_lm$a0
  
  out_SWE$week[ii,] <- weekly_avg_swe$week
  out_SWE$SWE[ii,] <- weekly_avg_swe$avg_swe
  
  
  
  #  end main Loop
}

dev.off()

write.table(out_T,file=paste(path,"Weekly_Sensitivity_to_T.txt", sep ="" ))
write.table(out_P,file=paste(path,"Weekly_Sensitivity_to_P.txt", sep ="" ))
write.table(out_Q,file=paste(path,"Weekly_Sensitivity_to_Q.txt", sep ="" ))
write.table(out_SWE,file=paste(path,"Weekly_SWE.txt", sep ="" ))



### 

# Common parameters
xx <- seq(1, 52) * 7 - 3
xtick <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
xtlab <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D", NA)
plot_settings <- list(
  type = "l", col = "black", lwd = 2,
  xlab = "", ylab = "Spearman rank correlation",
  ylim = c(-1, 1), xaxt = "n", xaxs = "i"
)

# Custom plotting function
create_plot <- function(correlations, filename) {
  png(filename,
      width = 6, height = 4,
      units = "in", res = 300) # High-resolution compact format
  
  plot(xx, correlations,
       type = plot_settings$type,
       pch = plot_settings$pch,
       col = plot_settings$col,
       lwd = plot_settings$lwd,
       xlab = plot_settings$xlab,
       ylab = plot_settings$ylab,
       ylim = plot_settings$ylim,
       xaxt = plot_settings$xaxt,
       xaxs = plot_settings$xaxs)
  
  axis(1, at = xtick + 13, labels = xtlab, lwd.ticks = 0)
  axis(1, at = xtick, labels = rep("", length(xtick)))
  abline(h = 0, lty = 2, col = "grey")
  
  dev.off()
}

# Loop over columns and calculate correlations for BFI vs coeff
correlations <- numeric(52)
for (week in seq_len(ncol(out_Q$coeff))) {
  correlations[week] <- cor(out_BFI$BFI, out_Q$coeff[, week], method = "spearman", use = "complete.obs")
}
create_plot(correlations, "D:/R/weekly_streamflow_sensitivities/correlation_plot_BFI_bQ.png")

# Loop over columns and calculate correlations for BFI vs var_imp
correlations <- numeric(52)
for (week in seq_len(ncol(out_Q$var_imp))) {
  correlations[week] <- cor(out_BFI$BFI, out_Q$var_imp[, week], method = "spearman", use = "complete.obs")
}
create_plot(correlations, "D:/R/weekly_streamflow_sensitivities/correlation_plot_BFI_varimp.png")

# Loop over columns and calculate correlations for SWE vs bQ
correlations <- numeric(52)
for (week in seq_len(ncol(out_Q$coeff))) {
  correlations[week] <- cor(out_SWE$SWE, out_Q$coeff[, week], method = "spearman", use = "complete.obs")
}
create_plot(correlations, "D:/R/weekly_streamflow_sensitivities/correlation_plot_SWE_bQ.png")

# Loop over columns and calculate correlations for SWE vs bT
correlations <- numeric(52)
for (week in seq_len(ncol(out_T$coeff))) {
  correlations[week] <- cor(out_SWE$SWE, out_T$coeff[, week], method = "spearman", use = "complete.obs")
}
create_plot(correlations, "D:/R/weekly_streamflow_sensitivities/correlation_plot_SWE_bT.png")


###

# Calculate Spearman correlations for each catchment
n_catchments <- ncol(out_Q$coeff)  # Number of catchments (columns)
correlations <- numeric(n_catchments)  # Initialize vector to store correlations

for (catchment in 1:n_catchments) {
  correlations[catchment] <- cor(out_Q$coeff[, catchment], out_SWE$SWE[, catchment], 
                                 method = "spearman", use = "complete.obs")
}

# Sort correlations by elevation
sorted_indices <- order(out_Q$elev_ezg)  # Sort indices by elevation
sorted_correlations <- correlations[sorted_indices]
sorted_elevation <- out_Q$elev_ezg[sorted_indices]

# Plot correlations sorted by elevation
png("D:/R/weekly_streamflow_sensitivities/correlation_sorted_by_elevation.png", 
    width = 4, height = 6, units = "in", res = 300)

plot(sorted_correlations, sorted_elevation, type = "p", col = "black", pch = 19,
     ylab = "Elevation (m)", xlab = "Spearman rank correlation",
     xlim = c(-1, 1))
abline(v = 0, lty = 2, col = "grey")

dev.off()

# Calculate correlation
avg_coeff <- rowMeans(out_Q$coeff)
cor_result <- cor.test(avg_coeff, out_BFI$BFI, method = "spearman")
print(paste("Spearman correlation coefficient:", round(cor_result$estimate, 3)))

# Partial correlation
avg_coeff <- rowMeans(out_Q$coeff)
max_SWE <- apply(out_SWE$SWE, 1, max) 
partial_cor <- pcor.test(avg_coeff, out_BFI$BFI, max_SWE, method = "spearman")
cat("Partial Spearman Correlation (controlling for SWE):\n")
cat("Correlation coefficient:", round(partial_cor$estimate, 3), "\n")
