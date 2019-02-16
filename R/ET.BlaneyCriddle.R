ET.BlaneyCriddle <- function(data, constants, ts = "daily", solar = "sunshine hours",
    height = F, message = "yes", save.csv = "yes", ...) 
{
    # class(data) <- funname

    # Check of specific data requirement
    if (is.null(data$Tmax) | is.null(data$Tmin)) {
        stop("Required data missing for 'Tmax' and 'Tmin', or 'Temp'")
    }
    if (solar == "sunshine hours" & is.null(data$n)) {
        # sunshine hour data is required
        stop("Required data missing for 'n'")
    } else if (solar == "cloud") {
        if (is.null(data$n)) {
            # for alternative calculation of sunshine hours using cloud
            # cover
            stop("Required data missing for 'Cd'")
        }
        if (is.null(data$u2) & is.null(data$uz)) {
            stop("Required data missing for 'uz' or 'u2'")
        }
        if (is.null(data$RHmin)) {
            stop("Required data missing for 'RHmin'")
        }
    }
    if (solar == "data" | solar == "monthly precipitation") {
        stop("Only 'sunshine hours' and 'cloud' are accepted because estimations of sunshine hours is required")
    }
    # Calculating mean temperature
    Ta <- (data$Tmax + data$Tmin)/2  # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998.

    # Calculations from data and constants for Blaney and Criddle
    delta2 <- 0.409 * sin(2 * pi/365 * data$J - 1.39)  # solar dedication (S3.7)
    w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
    N <- 24/pi * w_s  # calculating daily values

    # Wind speed
    if (is.null(data$u2)) {
        u2 <- data$uz * 4.87/log(67.8 * constants$z - 5.42)  # Equation S5.20 for PET formulations other than Penman
    } else {
        u2 <- data$u2
    }

    bvar <- constants$e0 + constants$e1 * data$RHmin + constants$e2 *
        data$n/N + constants$e3 * u2 + constants$e4 * data$RHmin *
        data$n/N + constants$e5 * data$RHmin * u2  # undefined working variable (Allena and Pruitt, 1986; Shuttleworth, 1992) (S9.8)
    N.annual <- ave(N, format(time(N), "%y"), FUN = sum)  # Annual sum of maximum sunshine hours
    # chech if data from first/last years is incomplete, and
    # adjust N.annual values for incomplete years first year a
    # normal year
    if (data$J[1] != 1 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[1]/4) ==
        FALSE) {
        N.annual[floor(as.numeric(as.yearmon(data$Date.daily))) ==
            floor(as.numeric(as.yearmon(data$Date.daily)))[1]] <- sum(24/pi *
            acos(-tan(constants$lat_rad) * tan(0.409 * sin(2 *
                pi/365 * c(1:365) - 1.39))))
    }
    if (data$J[1] != 1 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[1]/4) ==
        TRUE) {
        # first year a leap year
        N.annual[floor(as.numeric(as.yearmon(data$Date.daily))) ==
            floor(as.numeric(as.yearmon(data$Date.daily)))[1]] <- sum(24/pi *
            acos(-tan(constants$lat_rad) * tan(0.409 * sin(2 *
                pi/365 * c(1:366) - 1.39))))
    }
    if (data$J[length(data$J)] != 365 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]/4) ==
        FALSE) {
        # last year a normal year
        N.annual[floor(as.numeric(as.yearmon(data$Date.daily))) ==
            floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]] <- sum(24/pi *
            acos(-tan(constants$lat_rad) * tan(0.409 * sin(2 *
                pi/365 * c(1:365) - 1.39))))
    }
    if (data$J[length(data$J)] != 366 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]/4) ==
        TRUE) {
        # first year a leap year
        N.annual[floor(as.numeric(as.yearmon(data$Date.daily))) ==
            floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]] <- sum(24/pi *
            acos(-tan(constants$lat_rad) * tan(0.409 * sin(2 *
                pi/366 * c(1:366) - 1.39))))
    }
    p_y <- 100 * data$n/N.annual  # percentage of actual daytime hours for the day comparing to the annual sum of maximum sunshine hours


    ET_BC.Daily <- (0.0043 * data$RHmin - data$n/N - 1.41) +
        bvar * p_y * (0.46 * Ta + 8.13)  # Blaney-Criddle Reference Crop evapotranspiration (mm.day^-1) (S9.7)

    ET.Daily <- ET_BC.Daily
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

    if (height == T) {
        ET_BC.Daily = ET_BC.Daily * (1 + 0.1 * constants$Elev/1000)  # with adjustment for site elevation by Allen and Pruitt (1986) (S9.9)
    }

    # Generate summary message for results
    ET_formulation <- "Blaney-Criddle"
    ET_type <- "Reference Crop ET"
    if (solar == "sunshine hours") {
        message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
    } else if (solar == "cloud") {
        message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
    } else {
        message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
    }

    if (height == T) {
        message3 <- "Height adjustment has been applied to calculated Blaney-Criddle reference crop evapotranspiration"
    } else {
        message3 <- "No height adjustment has been applied to calculated Blaney-Criddle reference crop evapotranspiration"
    }
    results <- list(ET.Daily = ET.Daily, ET.Monthly = ET.Monthly,
        ET.Annual = ET.Annual, ET.MonthlyAve = ET.MonthlyAve,
        ET.AnnualAve = ET.AnnualAve, ET_formulation = ET_formulation,
        ET_type = ET_type, message1 = message1, message3 = message3)

    if (ts == "daily") {
        res_ts <- ET.Daily
    } else if (ts == "monthly") {
        res_ts <- ET.Monthly
    } else if (ts == "annual") {
        res_ts <- ET.Annual
    }
    if (message == "yes") {
        message(ET_formulation, " ", ET_type)
        message("Evaporative surface: reference crop")
        message(message1)
        message(message3)

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
            write.table(as.character(namer), file = "ET_BlaneyCriddle.csv",
                dec = ".", quote = FALSE, col.names = FALSE,
                row.names = F, append = TRUE, sep = ",")
            write.table(data.frame(get(namer, results)), file = "ET_BlaneyCriddle.csv",
                col.names = F, append = T, sep = ",")
        }
        invisible(results)
    } else {
        return(results)
    }
}
