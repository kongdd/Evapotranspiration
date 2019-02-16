ET.ChapmanAustralian <- function(data, constants, ts="daily",PenPan=T, solar="sunshine hours", alpha=0.23,
                                 message = "yes",save.csv ="yes",...) 
{
  #class(data) <- funname
  
  # Check of specific data requirement
  if (PenPan == TRUE) { # Calculate Class-A pan evaporation using PenPan formula
    if (is.null(data$Tmax)|is.null(data$Tmin)) { 
      stop("Required data missing for 'Tmax' and 'Tmin', or 'Temp'")
    }
    if (is.null(data$va)|is.null(data$vs)) {
      if (is.null(data$RHmax)|is.null(data$RHmin)) {
        stop("Required data missing: need either 'va' and 'vs', or 'RHmax' and 'RHmin' (or 'RH')")
      }
    } 
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("Required data missing for 'uz' or 'u2'")
    }
    if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
      stop("Required data missing for 'Rs'")
    } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
      stop("Required data missing for 'n'")
    } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
      stop("Required data missing for 'Cd'")
    } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
      stop("Required data missing for 'Precip'")
    } 
  }
  if (PenPan == FALSE & is.null(data$Epan)) { # for using Class-A pan evaporation data
    stop("Required data missing for 'Epan'")
  }
  # check user-input albedo
  if (PenPan == TRUE) {
    if (is.na(as.numeric(alpha))) {
      stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
    }
    if (!is.na(as.numeric(alpha))) {
      if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
        stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
      }
    } 
  }

  # Calculations from data and constants for daily equivalent Penman-Monteith potential evaporation 
  A_p <- 0.17 + 0.011 * abs(constants$lat) # constant (S13.2)
  B_p <- 10 ^ (0.66 - 0.211 * abs(constants$lat)) # constants (S13.3)
  
  if (PenPan == TRUE) {
    # Calculating mean temperature 
    Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
    
    if (!is.null(data$va) & !is.null(data$vs)) {
      vabar <- data$va
      vas <- data$vs
    } else {
      # Saturated vapour pressure
      vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
      vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
      vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
      
      # Vapour pressure
      vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
    }    
    # estimating class-A pan evaporation using PenPan model 
    P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
    delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
    gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
    d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
    delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
    w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
    N <- 24/pi * w_s # calculating daily values
    R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
    R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
    
    if (solar == "data") {
      R_s <- data$Rs
    } else if (solar == "monthly precipitation") {
      # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
      R_s <- (0.85 - 0.047*data$Cd)*R_a 
    } else {
      # calculate R_s from sunshine hours - data or estimation using cloudness
      R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
    }
    
    # Wind speed
    if (is.null(data$u2)) {
      u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
    } else {
      u2 <- data$u2
    }
    
    R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
    # For short grass
    P_rad <- 1.32 + 4 * 10^(-4) * abs(constants$lat) + 8 * 10^(-5) * (constants$lat)^2 # pan radiation factor (S6.6)
    f_dir <- -0.11 + 1.31 * R_s / R_a # fraction of R_S that os direct (S6.5)
    R_span <- (f_dir * P_rad + 1.42 * (1 - f_dir) + 0.42 * alpha) * R_s # total shortwave radiation received (S6.4)
    R_npan <-(1 - constants$alphaA) * R_span - R_nl # net radiation at the pan (S6.3)
    f_pan_u <-1.201 + 1.621 * u2 # (S6.2)
    
    Epan <- delta / (delta + constants$ap * gamma) * R_npan / constants$lambda + constants$ap * gamma / (delta + constants$ap * gamma) * f_pan_u * (vas - vabar) # PenPan estimation of Class-A pan evaporation (S6.1)
    ET_eqPM.Daily <- A_p * Epan + B_p # daily equivalent Penman-Monteith potential evaporation (mm.day^-1)
  } else if (PenPan == FALSE & is.null(data$Epan)) {
    stop("No data available for Class-A pan evaporation ")
  } else {
    ET_eqPM.Daily <- A_p * data$Epan + B_p # daily equivalent Penman-Monteith potential evaporation (mm.day^-1)
  }
 
  ET.Daily <- ET_eqPM.Daily
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
  ET_formulation <- "Chapman"
  ET_type <- "Potential ET"
  if (PenPan == TRUE) {
    Surface <- paste("user-defined, albedo =", alpha)
  } else {
    Surface <- paste("not specified, actual Class-A pan evaporation data is used")
  }
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (PenPan == TRUE) {
    message5 <- "PenPan formulation has been used to estimate Class-A pan evaporation for the calculation of potential evapotranspiration"
  } else {
    message5 <- "Class-A pan evaporation has been used for the calculation of potential evapotranspiration"
  }
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message5=message5)
  if (ts=="daily") {
    res_ts <-  ET.Daily
  } else if (ts=="monthly") {
    res_ts <- ET.Monthly
  } else if (ts=="annual") {
    res_ts <- ET.Annual
  }
  
  if (message == "yes") {
    # Generate summary message for results
    
    message(ET_formulation, " ", ET_type)
    message("Evaporative surface: ", Surface)
    message(message1)
    message(message5)
    
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
      write.table(as.character(namer), file="ET_ChapmanAustralian.csv", dec=".", 
                  quote=FALSE, col.names=FALSE, row.names=F, append=TRUE, sep=",")
      write.table(data.frame(get(namer,results)), file="ET_ChapmanAustralian.csv", col.names=F, append= T, sep=',' )
    }
    invisible(results)
  } else {
    return(results)
  }
 
}
