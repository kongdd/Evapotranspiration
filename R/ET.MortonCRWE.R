ET.MortonCRWE <- function (data, constants, ts = "monthly", est = "potential ET", 
                           solar = "sunshine hours", Tdew = T, alpha = NULL, 
                           message = "yes",save.csv ="yes",...) 
{
  constants$epsilonMo <- 0.97
  constants$fz <- 25
  constants$b0 <- 1.12
  constants$b1 <- 13
  constants$b2 <- 1.12
  
  variables <- Radiation(data, constants, ts, solar, Tdew, 
                         alpha, model = "CRWE")
  R_TC <- as.vector(variables$R_T)
  for (i in 1:length(R_TC)) {
    if (R_TC[i] < 0) {
      R_TC[i] <- 0
    }
    else {
      R_TC[i] <- R_TC[i]
    }
  }
  xiMo <- 1/(0.28 * (1 + variables$vD_Mo/variables$v_Mo) + 
               R_TC * variables$deltaMo/(variables$ptops * constants$gammaps * 
                                           (1/variables$ptops)^0.5 * constants$b0 * constants$fz * 
                                           (variables$v_Mo - variables$vD_Mo)))
  for (i in 1:length(xiMo)) {
    if (xiMo[i] < 1) {
      xiMo[i] <- 1
    }
    else {
      xiMo[i] <- xiMo[i]
    }
  }
  f_T <- (1/variables$ptops)^0.5 * constants$fz/xiMo
  lambdaMo1 <- constants$gammaps * variables$ptops + 4 * constants$epsilonMo * 
    constants$sigmaMo * (variables$T_Mo + 274)^3/f_T
  T_p <- variables$T_Mo
  for (i in 1:99999) {
    v_p <- 6.11 * exp((constants$alphaMo * T_p)/(T_p + constants$betaMo))
    delta_p <- constants$alphaMo * constants$betaMo * v_p/((T_p + 
                                                              constants$betaMo)^2)
    delta_T_p <- (R_TC/f_T + variables$vD_Mo - v_p + lambdaMo1 * 
                    (variables$T_Mo - T_p))/(delta_p + lambdaMo1)
    T_p <- T_p + delta_T_p
    if (abs(max(na.omit(delta_T_p))) < 0.01) 
      break
  }
  v_p <- 6.11 * exp((constants$alphaMo * T_p)/(T_p + constants$betaMo))
  delta_p <- constants$alphaMo * constants$betaMo * v_p/((T_p + 
                                                            constants$betaMo)^2)
  E_P.temp <- R_TC - lambdaMo1 * f_T * (T_p - variables$T_Mo)
  R_P <- E_P.temp + variables$ptops * constants$gammaps * 
    f_T * (T_p - variables$T_Mo)
  E_W.temp <- constants$b1 + constants$b2 * R_P/(1 + variables$ptops * 
                                                   constants$gammaps/delta_p)
  E_T_Mo.temp <- 2 * E_W.temp - E_P.temp
  E_P.temp <- 1/(constants$lambdaMo) * E_P.temp
  E_W.temp <- 1/(constants$lambdaMo) * E_W.temp
  E_T_Mo.temp <- 1/(constants$lambdaMo) * E_T_Mo.temp
  E_P <- E_P.temp * data$Ndays
  E_W <- E_W.temp * data$Ndays
  E_T_Mo <- E_T_Mo.temp * data$Ndays
  if (est == "potential ET") {
    ET_Mo.Monthly <- E_P
    ET_Mo.Average <- E_P.temp
    ET_type <- "Potential ET"
  }
  else if (est == "shallow lake ET") {
    ET_Mo.Monthly <- E_W
    ET_Mo.Average <- E_W.temp
    ET_type <- "Shallow Lake Evaporation"
  }
  ET.Daily <- NULL
  ET.Monthly <- ET_Mo.Monthly
  ET.Annual <- aggregate(ET.Monthly, floor(as.numeric(as.yearmon(data$Date.monthly, 
                                                                 "%m/%y"))), FUN = sum)
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.monthly)$mon):max(as.POSIXlt(data$Date.monthly)$mon)) {
    i = mon - min(as.POSIXlt(data$Date.monthly)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$mon == 
                                             mon])
  }
  for (year in min(as.POSIXlt(data$Date.monthly)$year):max(as.POSIXlt(data$Date.monthly)$year)) {
    i = year - min(as.POSIXlt(data$Date.monthly)$year) + 
      1
    ET.AnnualAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$year == 
                                            year])
  }
  ET_formulation <- "Morton CRWE"
  
  results <- list(ET.Daily = ET.Daily, ET.Monthly = ET.Monthly, 
                  ET.Annual = ET.Annual, ET.MonthlyAve = ET.MonthlyAve, 
                  ET.AnnualAve = ET.AnnualAve, ET_formulation = ET_formulation, 
                  ET_type = ET_type, message1 = variables$message1, message6 = variables$message6)
  if (ts == "monthly") {
    res_ts <- ET.Monthly
  }
  else if (ts == "annual") {
    res_ts <- ET.Annual
  }
  if (message =="yes") {
    message(ET_formulation, " ", ET_type)
    message(variables$message1)
    message(variables$message6)
    
    message("Timestep: ", ts)
    message("Units: mm")
    message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
    if (NA %in% res_ts) {
      message(length(res_ts), " ET estimates obtained; ", 
              length(which(is.na(res_ts))), " NA output entries due to missing data")
      message("Basic stats (NA excluded)")
      message("Mean: ", round(mean(res_ts, na.rm = T), digits = 2))
      message("Max: ", round(max(res_ts, na.rm = T), digits = 2))
      message("Min: ", round(min(res_ts, na.rm = T), digits = 2))
    }
    else {
      message(length(res_ts), " ET estimates obtained")
      message("Basic stats")
      message("Mean: ", round(mean(res_ts), digits = 2))
      message("Max: ", round(max(res_ts), digits = 2))
      message("Min: ", round(min(res_ts), digits = 2))
    }
  }
  if (save.csv == "yes") {
    for (i in 1:length(results)) {
      namer <- names(results[i])
      write.table(as.character(namer), file = "ET_MortonCRWE.csv", 
                  dec = ".", quote = FALSE, col.names = FALSE, row.names = F, 
                  append = TRUE, sep = ",")
      write.table(data.frame(get(namer, results)), file = "ET_MortonCRWE.csv", 
                  col.names = F, append = T, sep = ",")
    }
    invisible(results)
  } else {
    return(results)
  }
  
}