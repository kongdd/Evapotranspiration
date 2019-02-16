ET.MortonCRAE <- function (data, constants, ts = "monthly", est = "potential ET", 
                           solar = "sunshine hours", Tdew = T, alpha = NULL, 
                           message = "yes",save.csv ="yes",...) 
{
  variables <- Radiation(data, constants, ts, solar, Tdew, 
                         alpha)
  R_T <- (1 - variables$alpha_Mo) * variables$G_Mo - variables$B_Mo
  R_TC <- as.vector(R_T)
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
    delta_T_p <- (R_T/f_T + variables$vD_Mo - v_p + lambdaMo1 * 
                    (variables$T_Mo - T_p))/(delta_p + lambdaMo1)
    T_p <- T_p + delta_T_p
    if (abs(max(na.omit(delta_T_p))) < 0.01) 
      break
  }
  v_p <- 6.11 * exp((constants$alphaMo * T_p)/(T_p + constants$betaMo))
  delta_p <- constants$alphaMo * constants$betaMo * v_p/((T_p + 
                                                            constants$betaMo)^2)
  E_TP.temp <- R_T - lambdaMo1 * f_T * (T_p - variables$T_Mo)
  R_TP <- E_TP.temp + variables$ptops * constants$gammaps * 
    f_T * (T_p - variables$T_Mo)
  E_TW.temp <- constants$b1 + constants$b2 * R_TP/(1 + variables$ptops * 
                                                     constants$gammaps/delta_p)
  E_T_Mo.temp <- 2 * E_TW.temp - E_TP.temp
  E_TP.temp <- 1/(constants$lambdaMo) * E_TP.temp
  E_TW.temp <- 1/(constants$lambdaMo) * E_TW.temp
  E_T_Mo.temp <- 1/(constants$lambdaMo) * E_T_Mo.temp
  E_TP <- E_TP.temp * data$Ndays
  E_TW <- E_TW.temp * data$Ndays
  E_T_Mo <- E_T_Mo.temp * data$Ndays
  if (est == "potential ET") {
    ET_Mo.Monthly <- E_TP
    ET_Mo.Average <- E_TP.temp
    ET_type <- "Potential ET"
  }
  else if (est == "wet areal ET") {
    ET_Mo.Monthly <- E_TW
    ET_Mo.Average <- E_TW.temp
    ET_type <- "Wet-environment Areal ET"
  }
  else if (est == "actual areal ET") {
    ET_Mo.Monthly <- E_T_Mo
    ET_Mo.Average <- E_T_Mo.temp
    ET_type <- "Actual Areal ET"
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
  ET_formulation <- "Morton CRAE"
  
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
  
  if (message == "yes") {
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
  if (save.csv== "yes") {
    for (i in 1:length(results)) {
      namer <- names(results[i])
      write.table(as.character(namer), file = "ET_MortonCRAE.csv", 
                  dec = ".", quote = FALSE, col.names = FALSE, row.names = F, 
                  append = TRUE, sep = ",")
      write.table(data.frame(get(namer, results)), file = "ET_MortonCRAE.csv", 
                  col.names = F, append = T, sep = ",")
    }
    invisible(results)
  } else {
    return(results)
  }
}
