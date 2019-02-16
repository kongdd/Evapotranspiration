#'
#' @references Oudin et al., 2005
ET.Romanenko <- function(data, constants = NULL, ts = "daily",
    message = "yes", save.csv = "yes", ...) {
    # Check of specific data requirement
    if (is.null(data$Tmax) | is.null(data$Tmin)) {
        stop("Required data missing for 'Tmax' and 'Tmin', or 'Temp'")
    }
    # Check of specific data requirement
    if (is.null(data$va) | is.null(data$vs)) {
        if (is.null(data$RHmax) | is.null(data$RHmin)) {
            stop("Required data missing: need either 'va' and 'vs', or 'RHmax' and 'RHmin' (or 'RH')")
        }
    }
    Ta <- (data$Tmax + data$Tmin)/2

    if (!is.null(data$va) & !is.null(data$vs)) {
        vabar <- data$va
        vas <- data$vs
    } else {
        # Saturated vapour pressure
        vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax/(data$Tmax +
            237.3))  # Equation S2.5
        vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin/(data$Tmin +
            237.3))  # Equation S2.5
        vas <- (vs_Tmax + vs_Tmin)/2  # Equation S2.6

        # Vapour pressure
        vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2  # Equation S2.7
    }

    ET_Romanenko.Daily <- 4.5 * (1 + Ta/25)^2 * (1 - vabar/vas)
    ET.Daily <- ET_Romanenko.Daily
    ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily,
        "%m/%y"), FUN = sum)
    ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily,
        "%m/%y"))), FUN = sum)

    ET.MonthlyAve <- ET.AnnualAve <- NULL
    for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)) {
        i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
        ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon ==
            mon])
    }
    for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)) {
        i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
        ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year ==
            year])
    }

    # Generate summary message for results
    ET_formulation <- "Romanenko"
    ET_type <- "Actual ET"

    results <- list(ET.Daily = ET.Daily, ET.Monthly = ET.Monthly,
        ET.Annual = ET.Annual, ET.MonthlyAve = ET.MonthlyAve,
        ET.AnnualAve = ET.AnnualAve, ET_formulation = ET_formulation,
        ET_type = ET_type)
    if (ts == "daily") {
        res_ts <- ET.Daily
    } else if (ts == "monthly") {
        res_ts <- ET.Monthly
    } else if (ts == "annual") {
        res_ts <- ET.Annual
    }

    if (message == "yes") {
        message(ET_formulation, " ", ET_type)
        message("Timestep: ", ts)
        message("Units: mm")
        message("Time duration: ", time(res_ts[1]), " to ", time(res_ts[length(res_ts)]))
        if (NA %in% res_ts) {
            message(length(res_ts), " ET estimates obtained; ",
                length(which(is.na(res_ts))), " NA output entries due to missing data")
            message("Basic stats (NA excluded)")
            message("Mean: ", round(mean(res_ts, na.rm = T),
                digits = 2))
            message("Max: ", round(max(res_ts, na.rm = T), digits = 2))
            message("Min: ", round(min(res_ts, na.rm = T), digits = 2))
        } else {
            message(length(res_ts), " ET estimates obtained")
            message("Basic stats")
            message("Mean: ", round(mean(res_ts), digits = 2))
            message("Max: ", round(max(res_ts), digits = 2))
            message("Min: ", round(min(res_ts), digits = 2))
        }
        # class(results) <- funname
    }
    if (save.csv == "yes") {
        # write to csv file
        for (i in 1:length(results)) {
            namer <- names(results[i])
            write.table(as.character(namer), file = "ET_Romanenko.csv",
                dec = ".", quote = FALSE, col.names = FALSE,
                row.names = F, append = TRUE, sep = ",")
            write.table(data.frame(get(namer, results)), file = "ET_Romanenko.csv",
                col.names = F, append = T, sep = ",")
        }
        invisible(results)
    } else {
        return(results)
    }

}