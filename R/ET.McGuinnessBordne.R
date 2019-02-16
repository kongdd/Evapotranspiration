ET.McGuinnessBordne <- function(data, constants, ts="daily", 
                                message = "yes",save.csv ="yes",...) 
{
  #class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax' and 'Tmin', or 'Temp'")
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 

  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)

  # estimating evapotranspiration using McGuinness-Bordne

  ET_MB.Daily <- R_a * (Ta + 5)/ (constants$lambda*68)  # McGuinness-Bordne daily evapotranspiration by McGuinness-Bordne  (1972) (mm.day^-1) (Oudin et al., 2005
  
  ET.Daily <- ET_MB.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "McGuinness-Bordne"
  ET_type <- "Potential ET"
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type)
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  if (message =="yes") {
    message(ET_formulation, " ", ET_type)
    
    message("Timestep: ", ts)
    message("Units: mm")
    message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
    if (NA %in% res_ts) {
      message(length(res_ts), " ET estimates obtained; ", length(which(is.na(res_ts))), " NA output entries due to missing data")
      message("Basic stats (NA excluded)")
      message("Mean: ",round(mean(res_ts,na.rm=T),digits=2))
      message("Max: ",round(max(res_ts,na.rm=T),digits=2))
      message("Min: ",round(min(res_ts,na.rm=T),digits=2))
    } else {
      message(length(res_ts), " ET estimates obtained")
      message("Basic stats")
      message("Mean: ",round(mean(res_ts),digits=2))
      message("Max: ",round(max(res_ts),digits=2))
      message("Min: ",round(min(res_ts),digits=2))
    }
    #class(results) <- funname
  }
  if (save.csv == "yes") {
    # write to csv file
    for (i in 1:length(results)) {
      namer <- names(results[i])
      write.table(as.character(namer), file="ET_McGuinnessBordne.csv", dec=".", 
                  quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
      write.table(data.frame(get(namer,results)), file="ET_McGuinnessBordne.csv", col.names=F, append= T, sep=',' )
    }
    invisible(results)
  } else {
    return(results)
  }
}