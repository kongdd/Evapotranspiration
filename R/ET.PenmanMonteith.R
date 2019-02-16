ET.PenmanMonteith <- function(data, constants, ts = "daily",
    solar = "sunshine hours", wind = "yes", crop = "short", message = "yes",
    save.csv = "yes", ...) 
{
    # class(data) <- funname

    # Check of specific data requirement
    if (is.null(data$Tmax) | is.null(data$Tmin)) {
        stop("Required data missing for 'Tmax' and 'Tmin', or 'Temp'")
    }
    if (is.null(data$va) | is.null(data$vs)) {
        if (is.null(data$RHmax) | is.null(data$RHmin)) {
            stop("Required data missing: need either 'va' and 'vs', or 'RHmax' and 'RHmin' (or 'RH')")
        }
    }
    if (wind == "yes") {
        # wind data is required
        if (is.null(data$u2) & is.null(data$uz)) {
            stop("Required data missing for 'uz' or 'u2'")
        }
    }

    if (solar == "data" & is.null(data$Rs)) {
        # solar radiation data is required
        stop("Required data missing for 'Rs'")
    } else if (solar == "sunshine hours" & is.null(data$n)) {
        # for alternative calculation of solar radiation with
        # sunshine hour
        stop("Required data missing for 'n'")
    } else if (solar == "cloud" & is.null(data$Cd)) {
        # for alternative calculation of sunshine hours using cloud
        # cover
        stop("Required data missing for 'Cd'")
    } else if (solar == "monthly precipitation" & is.null(data$Precip)) {
        # for alternative calculation of cloudiness using monthly
        # precipitation
        stop("Required data missing for 'Precip'")
    }

    if (wind != "yes" & wind != "no") {
        stop("Please choose if actual data will be used for wind speed from wind = 'yes' and wind = 'no'")
    }
    # check user-input crop type and specify albedo
    if (wind == "yes") {
        if (crop != "short" & crop != "tall") {
            stop("Please enter 'short' or 'tall' for the desired reference crop type")
        } else {
            alpha <- 0.23  # albedo for both short and tall crop
            if (crop == "short") {
                z0 <- 0.02  # roughness height for short grass
            } else {
                z0 <- 0.1  # roughness height for tall grass
            }
        }
    } else {
        z0 <- 0.02  # roughness height for short grass
        alpha <- 0.25  # semi-desert short grass - will not be used for calculation - just informative
    }

    # Calculating mean temperature
    Ta <- (data$Tmax + data$Tmin)/2  # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998.

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
    # Calculations from data and constants for Penman-Monteith
    # Reference Crop

    P <- 101.3 * ((293 - 0.0065 * constants$Elev)/293)^5.26  # atmospheric pressure (S2.10)
    delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta + 237.3)))/((Ta +
        237.3)^2)  # slope of vapour pressure curve (S2.4)
    gamma <- 0.00163 * P/constants$lambda  # psychrometric constant (S2.9)
    d_r2 <- 1 + 0.033 * cos(2 * pi/365 * data$J)  # dr is the inverse relative distance Earth-Sun (S3.6)
    delta2 <- 0.409 * sin(2 * pi/365 * data$J - 1.39)  # solar dedication (S3.7)
    w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
    N <- 24/pi * w_s  # calculating daily values
    R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) *
        sin(delta2) + cos(constants$lat_rad) * cos(delta2) *
        sin(w_s))  # extraterristrial radiation (S3.5)
    R_so <- (0.75 + (2 * 10^-5) * constants$Elev) * R_a  # clear sky radiation (S3.4)

    if (solar == "data") {
        R_s <- data$Rs
    } else if (solar != "monthly precipitation") {
        # calculate R_s from sunshine hours - data or estimation
        # using cloudness
        R_s <- (constants$as + constants$bs * (data$n/N)) * R_a  # estimated incoming solar radiation (S3.9)
    } else {
        # calculate R_s from cloudness estimated from monthly
        # precipitation (#S3.14)
        R_s <- (0.85 - 0.047 * data$Cd) * R_a
    }

    R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax +
        273.2)^4 + (data$Tmin + 273.2)^4)/2 * (1.35 * R_s/R_so -
        0.35)  # estimated net outgoing longwave radiation (S3.3)
    R_nsg <- (1 - alpha) * R_s  # net incoming shortwave radiation (S3.2)
    R_ng <- R_nsg - R_nl  # net radiation (S3.1)

    if (wind == "yes") {
        # Wind speed
        if (is.null(data$u2)) {
            u2 <- data$uz * 4.87/log(67.8 * constants$z - 5.42)  # Equation S5.20 for PET formulations other than Penman
        } else {
            u2 <- data$u2
        }

        if (crop == "short") {
            r_s <- 70  # will not be used for calculation - just informative
            CH <- 0.12  # will not be used for calculation - just informative
            ET_RC.Daily <- (0.408 * delta * (R_ng - constants$G) +
                gamma * 900 * u2 * (vas - vabar)/(Ta + 273))/(delta +
                gamma * (1 + 0.34 * u2))  # FAO-56 reference crop evapotranspiration from short grass (S5.18)
        } else {
            r_s <- 45  # will not be used for calculation - just informative
            CH <- 0.5  # will not be used for calculation - just informative
            ET_RC.Daily <- (0.408 * delta * (R_ng - constants$G) +
                gamma * 1600 * u2 * (vas - vabar)/(Ta + 273))/(delta +
                gamma * (1 + 0.38 * u2))  # ASCE-EWRI standardised Penman-Monteith for long grass (S5.19)
        }
        ET.Daily <- ET_RC.Daily
        ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily,
            "%m/%y"), FUN = sum)
        ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily,
            "%m/%y"))), FUN = sum)

    } else {
        # mean relative humidity
        RHmean <- (data$RHmax + data$RHmin)/2

        R_s.Monthly <- aggregate(R_s, as.yearmon(data$Date.daily,
            "%m/%y"), mean)
        R_a.Monthly <- aggregate(R_a, as.yearmon(data$Date.daily,
            "%m/%y"), mean)
        Ta.Monthly <- aggregate(Ta, as.yearmon(data$Date.daily,
            "%m/%y"), mean)
        RHmean.Monthly <- aggregate(RHmean, as.yearmon(data$Date.daily,
            "%m/%y"), mean)
        # ET_RC.Daily <- matrix(NA,length(data$date.Daily),1)
        ET_RC.Monthly <- 0.038 * R_s.Monthly * sqrt(Ta.Monthly +
            9.5) - 2.4 * (R_s.Monthly/R_a.Monthly)^2 + 0.075 *
            (Ta.Monthly + 20) * (1 - RHmean.Monthly/100)  # Reference crop evapotranspiration without wind data by Valiantzas (2006) (S5.21)
        ET_RC.Daily <- data$Tmax
        for (cont in 1:length(data$i)) {
            ET_RC.Daily[(((as.numeric(as.yearmon(time(ET_RC.Daily)))) -
                floor(as.numeric(as.yearmon(time(ET_RC.Daily))))) *
                12 + 1) == data$i[cont]] <- ET_RC.Monthly[cont]
        }

        ET.Daily <- ET_RC.Daily
        ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily,
            "%m/%y"), FUN = sum)
        ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.monthly,
            "%m/%y"))), FUN = sum)
    }

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
    if (wind == "no") {
        ET_formulation <- "Penman-Monteith (without wind data)"
        ET_type <- "Reference Crop ET"
        Surface <- paste("short grass, albedo =", alpha, "; roughness height =",
            z0, "m")
    } else {
        if (crop == "short") {
            ET_formulation <- "Penman-Monteith FAO56"
            ET_type <- "Reference Crop ET"
            Surface <- paste("FAO-56 hypothetical short grass, albedo =",
                alpha, "; surface resistance =", r_s, "sm^-1; crop height =",
                CH, " m; roughness height =", z0, "m")
        } else {
            ET_formulation <- "Penman-Monteith ASCE-EWRI Standardised"
            ET_type <- "Reference Crop ET"
            Surface <- paste("ASCE-EWRI hypothetical tall grass, albedo =",
                alpha, "; surface resistance =", r_s, "sm^-1; crop height =",
                CH, " m; roughness height =", z0, "m")
        }
    }
    if (solar == "data") {
        message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
    } else if (solar == "sunshine hours") {
        message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
    } else if (solar == "cloud") {
        message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
    } else {
        message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
    }

    if (wind == "yes") {
        message2 <- "Wind data have been used for calculating the reference crop evapotranspiration"
    } else {
        message2 <- "Alternative calculation for reference crop evapotranspiration without wind data have been performed"
    }

    results <- list(ET.Daily = ET.Daily, ET.Monthly = ET.Monthly,
        ET.Annual = ET.Annual, ET.MonthlyAve = ET.MonthlyAve,
        ET.AnnualAve = ET.AnnualAve, ET_formulation = ET_formulation,
        ET_type = ET_type, message1 = message1, message2 = message2)

    if (ts == "daily") {
        res_ts <- ET.Daily
    } else if (ts == "monthly") {
        res_ts <- ET.Monthly
    } else if (ts == "annual") {
        res_ts <- ET.Annual
    }
    if (message == "yes") {

        message(ET_formulation, " ", ET_type)
        message("Evaporative surface: ", Surface)
        message(message1)
        message(message2)


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
    }

    # class(results) <- funname
    if (save.csv == "yes") {
        # write to csv file
        for (i in 1:length(results)) {
            namer <- names(results[i])
            write.table(as.character(namer), file = "ET_PenmanMonteith.csv",
                dec = ".", quote = FALSE, col.names = FALSE,
                row.names = F, append = TRUE, sep = ",")
            write.table(data.frame(get(namer, results)), file = "ET_PenmanMonteith.csv",
                col.names = F, append = T, sep = ",")
        }
        invisible(results)
    } else {
        return(results)
    }
}
