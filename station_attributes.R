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

library(dplyr)
library(lubridate)

# calculate station averages

path_p <- "D:/R/weekly_streamflow_sensitivities/data/"	
path_t <- "D:/R/weekly_streamflow_sensitivities/data_temp/"
path <- "D:/R/weekly_streamflow_sensitivities/Abflussdaten/"
files <- read.csv(paste(path,"Klima_Abfluss_Vergleich.csv", sep =""), header=TRUE)

nr <- length(files$Q_Nr)

results <- files										# Temperature sensitivity
results$Tm <- array(NA, nr)				# mean T
results$Pm <- array(NA, nr)					# mean P
results$snowfrac <- array(NA, nr)					# snow fraction
results$Tm_weekly <- matrix(NA, nr, 53)				# mean weekly T
results$Tm_weekly_min <- matrix(NA, nr, 53)				# min weekly T

#  Main Loop

for (ii in 1:nr)   # nr
{
  
  loc = files$C_Nr[ii]   # Klima Daten file
  
  # read precip data
  FILE   	<- paste(path_p,"order_17129_",loc,"_rhs150d0_1_data.txt", sep ="")
  data <- read.table(FILE, header = FALSE, skip=3, colClasses = "character")
  date <- as.Date((data[,2]),"%Y%m%d")
  data <- data[,3]
  NAvalues 	   	<- (data == "-") 
  data[NAvalues]		<- NA
  data <- as.numeric(data)
  a1 <- data.frame(date = date,P_obs = data) 
  a1$year <- as.numeric(format(a1$date, "%Y"))
  a1 <- subset(a1, year >= 1900 & year <= 2014)
  a1 <- data.frame(date = date,P_obs = data) 
  results$Pm[ii] <- mean.na(a1$P_obs)*365
  
  # read temp data 
  FILE   	<- paste(path_t,"order_17130_",loc,"_ths200d0_1_data.txt", sep ="")
  data <- read.table(FILE, header = FALSE, skip=3, colClasses = "character")
  date <- as.Date((data[,2]),"%Y%m%d")
  data <- data[,3]
  NAvalues 	   	<- (data == "-") 
  data[NAvalues]		<- NA
  data <- as.numeric(data)
  a2 <- data.frame(date = date,T_obs = data) 
  a2$year <- as.numeric(format(a2$date, "%Y"))
  a2 <- subset(a2, year >= 1900 & year <= 2014)
  a2 <- data.frame(date = date,T_obs = data) 
  results$Tm[ii] <- mean.na(a2$T_obs)
  
  # weekly temperature
  # TODO calc weakly mean first...
  #a2 <- a2 %>%
  #  mutate(week = week(date))
  #Tm_weekly <- a2 %>%
  #  group_by(week) %>%
  #  summarise(Tm_weekly = min(T_obs, na.rm = TRUE))  # na.rm = TRUE to handle missing values
  #results$Tm_weekly[ii,] <- Tm_weekly$Tm_weekly
  
  ###
  a2 <- a2 %>%
    mutate(week = week(date), year = year(date))  # add year to distinguish weeks across years
  # calculate weekly averages
  Tm_weekly_year <- a2 %>%
    group_by(year, week) %>%  #group by both year and week to avoid overlap across years
    summarise(Tm_weekly_year = mean(T_obs, na.rm = TRUE), .groups = 'drop')
  # calculate something using the weekly averages, e.g. mean
  Tm_weekly_avg <- Tm_weekly_year %>%
    group_by(week) %>%
    summarise(Tm_weekly_year = mean(Tm_weekly_year, na.rm = TRUE))
  results$Tm_weekly[ii,] <- Tm_weekly_avg$Tm_weekly_year
  # or min
  Tm_weekly_min <- Tm_weekly_year %>%
    group_by(week) %>%
    summarise(Tm_weekly_year = min(Tm_weekly_year, na.rm = TRUE))
  results$Tm_weekly_min[ii,] <- Tm_weekly_min$Tm_weekly_year
  
  # snow frac
  P_snow <- array(0, length(a2$T_obs))
  P_snow[a2$T_obs < 0 & !is.na(a2$T_obs)] <- a1$P_obs[a2$T_obs < 0 & !is.na(a2$T_obs)]
  results$snowfrac[ii] <- sum.na(P_snow) / sum.na(a1$P_obs)
  
  #a <- merge(a1, a2, by ="date")
  
}

### elevation and temperature
pdf("D:/R/weekly_streamflow_sensitivities/temp_vs_elevation.pdf", width = 4, height = 4.5)  
# plot
plot(results$Tm, results$elev, 
     xlab = "Temperature [°C]", 
     ylab = "Elevation [m a.s.l.]", 
     pch = 19,  
     col = "black", 
     cex = 1.0)
model <- lm(results$elev ~ results$Tm)
abline(model, col = "black", lwd = 2)
# correlation
correlation <- cor(results$Tm, results$elev)
print(round(correlation,3))
# predict catchment info
model <- lm(results$Tm ~ results$elev)
results$temp_ezg <- coef(model)[1] + coef(model)[2] * results$elev_ezg
dev.off()

### elevation and snow fraction
pdf("D:/R/weekly_streamflow_sensitivities/snowfraction_vs_elevation.pdf", width = 4, height = 4.5)
# plot
plot(results$snowfrac, results$elev, 
     xlab = "Snow fraction [-]", 
     ylab = "Elevation [m a.s.l.]", 
     pch = 19,  
     col = "black", 
     cex = 1.0)
model <- lm(results$elev ~ results$snowfrac)
abline(model, col = "black", lwd = 2)
# correlation
correlation <- cor(results$snowfrac, results$elev)
print(round(correlation,3))
# predict catchment info
model <- lm(results$snowfrac ~ results$elev)
results$snowfrac_ezg <- coef(model)[1] + coef(model)[2] * results$elev_ezg
dev.off()

### elevation and precipitation
pdf("D:/R/weekly_streamflow_sensitivities/precipitation_vs_elevation.pdf", width = 4, height = 4.5) 
# plot
plot(results$Pm, results$elev, 
     xlab = "Precipitation [mm/y]", 
     ylab = "Elevation [m a.s.l.]", 
     pch = 19,  
     col = "black", 
     cex = 1.0)
model <- lm(results$elev ~ results$Pm)
abline(model, col = "black", lwd = 2)
# correlation
correlation <- cor(results$Pm, results$elev)
print(round(correlation,3))
# predict catchment info
model <- lm(results$Pm ~ results$elev)
results$precip_ezg <- coef(model)[1] + coef(model)[2] * results$elev_ezg
dev.off()

### glacier cover and temperature
pdf("D:/R/weekly_streamflow_sensitivities/temp_vs_glaciers.pdf", width = 4, height = 4.5) 
# plot
plot(results$temp_ezg[results$Glacier>0], results$Glacier[results$Glacier>0], 
     xlab = "Temperature [°C]", 
     ylab = "Glacier cover [%]", 
     pch = 19,  
     col = "black", 
     cex = 1.0)
model <- lm(results$Glacier[results$Glacier>0] ~ results$temp_ezg[results$Glacier>0])
abline(model, col = "black", lwd = 2)
# correlation
correlation <- cor(results$temp_ezg[results$Glacier>0], results$Glacier[results$Glacier>0])
print(round(correlation,3))
dev.off()

### glacier cover and temperature
pdf("D:/R/weekly_streamflow_sensitivities/snowfrac_vs_glaciers.pdf", width = 4, height = 4.5) 
# plot
plot(results$snowfrac_ezg[results$Glacier>0], results$Glacier[results$Glacier>0], 
     xlab = "Snow fraction [-]", 
     ylab = "Glacier cover [%]", 
     pch = 19,  
     col = "black", 
     cex = 1.0)
model <- lm(results$Glacier[results$Glacier>0] ~ results$snowfrac_ezg[results$Glacier>0])
abline(model, col = "black", lwd = 2)
# correlation
correlation <- cor(results$snowfrac_ezg[results$Glacier>0], results$Glacier[results$Glacier>0])
print(round(correlation,3))
dev.off()

### temperature for specified elevation bands
elev <- seq(500,3500,500)
model <- lm(results$Tm ~ results$elev)
temp_predict <- coef(model)[1] + coef(model)[2] * elev
print(round(temp_predict,1))

temperature <- seq(-6,8,2)
model <- lm(results$elev ~ results$Tm)
elev_predict <- coef(model)[1] + coef(model)[2] * temperature
print(round(elev_predict,0))

### temperature for specified glacier covers
glaciers <- seq(0,60,10)
model <- lm(results$temp_ezg[results$Glacier>0] ~ results$Glacier[results$Glacier>0])
temp_predict <- coef(model)[1] + coef(model)[2] * glaciers
print(round(temp_predict,1))

temperature <- seq(-6,2,1)
model <- lm(results$Glacier[results$Glacier>0] ~ results$temp_ezg[results$Glacier>0])
glacier_predict <- coef(model)[1] + coef(model)[2] * temperature
print(round(glacier_predict,0))



### elevation and temperature for each week
pdf("D:/R/weekly_streamflow_sensitivities/weekly_temp.pdf", width = 4, height = 4.5)  

weekly_isotherm <- matrix(NA, 53, 2)	
temperature <- 0 #c(0,1,2)

for (w in 1:53) 
{
  
  plot(results$Tm_weekly[,w], results$elev, 
       xlab = "Temperature [°C]", 
       ylab = "Elevation [m a.s.l.]", 
       pch = 19,  
       col = "black", 
       cex = 1.0)
  
  model <- lm(results$elev ~ results$Tm_weekly[,w])
  abline(model, col = "black", lwd = 2)
  
  elev_predict <- coef(model)[1] + coef(model)[2] * temperature
  print(round(elev_predict,0))
  weekly_isotherm[w,1] <- elev_predict
    
  # correlation
  correlation <- cor(results$Tm_weekly[,w], results$elev)
  print(round(correlation,3))
  
  # weekly minima
  model <- lm(results$elev ~ results$Tm_weekly_min[,w])
  abline(model, col = "black", lwd = 2)
  
  elev_predict <- coef(model)[1] + coef(model)[2] * temperature
  print(round(elev_predict,0))
  weekly_isotherm[w,2] <- elev_predict
  
}
dev.off()


out_isotherm <- data.frame(week = seq(1, 53, 1), isotherm = weekly_isotherm)
pdf("D:/R/weekly_streamflow_sensitivities/weekly_temp_regime.pdf", width = 4, height = 4.5)  
plot(out_isotherm$week, out_isotherm$isotherm.1, 
     xlab = "Week [-]", 
     ylab = "Elevation [m a.s.l.]", 
     pch = 19,  
     col = "black", 
     cex = 1.0)
dev.off()

write.table(out_isotherm,file=paste(path,"0deg_isotherm_regime.txt", sep ="" ))
