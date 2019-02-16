ET.PenPan <- function(data, constants, ts = "daily", solar = "sunshine hours",
    alpha = 0.23, est = "potential ET", pan_coeff = 0.71, overest = F,
    message = "yes", save.csv = "yes", ...) {
    # class(data) <- funname
    if (is.null(data$Tmax) | is.null(data$Tmin)) {
        stop("Required data missing for 'Tmax' and 'Tmin', or 'Temp'")
    }
    if (is.null(data$va) | is.null(data$vs)) {
        if (is.null(data$RHmax) | is.null(data$RHmin)) {
            stop("Required data missing: need either 'va' and 'vs', or 'RHmax' and 'RHmin' (or 'RH')")
        }
    }
    if (is.null(data$u2) & is.null(data$uz)) {
        stop("Required data missing for 'uz' or 'u2'")
    }
    if (solar == "data" & is.null(data$Rs)) {
        stop("Required data missing for 'Rs'")
    } else if (solar == "sunshine hours" & is.null(data$n)) {
        stop("Required data missing for 'n'")
    } else if (solar == "cloud" & is.null(data$Cd)) {
        stop("Required data missing for 'Cd'")
    } else if (solar == "monthly precipitation" & is.null(data$Precip)) {
        stop("Required data missing for 'Precip'")
    }
    if (is.na(as.numeric(alpha))) {
        stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
    }
    if (!is.na(as.numeric(alpha))) {
        if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
            stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
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

    P <- 101.3 * ((293 - 0.0065 * constants$Elev)/293)^5.26
    delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta + 237.3)))/((Ta +
        237.3)^2)
    gamma <- 0.00163 * P/constants$lambda
    d_r2 <- 1 + 0.033 * cos(2 * pi/365 * data$J)
    delta2 <- 0.409 * sin(2 * pi/365 * data$J - 1.39)
    w_s <- acos(-tan(constants$lat_rad) * tan(delta2))
    N <- 24/pi * w_s
    R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) *
        sin(delta2) + cos(constants$lat_rad) * cos(delta2) *
        sin(w_s))
    R_so <- (0.75 + (2 * 10^-5) * constants$Elev) * R_a
    if (solar == "data") {
        R_s <- data$Rs
    } else if (solar != "monthly precipitation") {
        R_s <- (constants$as + constants$bs * (data$n/N)) * R_a
    } else {
        R_s <- (0.85 - 0.047 * data$Cd) * R_a
    }
    if (is.null(data$u2)) {
        u2 <- data$uz * 4.87/log(67.8 * constants$z - 5.42)
    } else {
        u2 <- data$u2
    }
    R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax +
        273.2)^4 + (data$Tmin + 273.2)^4)/2 * (1.35 * R_s/R_so -
        0.35)
    P_rad <- 1.32 + 4 * 10^(-4) * abs(constants$lat) + 8 * 10^(-5) *
        (constants$lat)^2
    f_dir <- -0.11 + 1.31 * R_s/R_a
    R_span <- (f_dir * P_rad + 1.42 * (1 - f_dir) + 0.42 * alpha) *
        R_s
    R_npan <- (1 - constants$alphaA) * R_span - R_nl
    f_pan_u <- 1.201 + 1.621 * u2
    Epenpan.Daily <- delta/(delta + constants$ap * gamma) * R_npan/constants$lambda +
        constants$ap * gamma/(delta + constants$ap * gamma) *
            f_pan_u * (vas - vabar)
    if (overest == TRUE) {
        if (est == "potential ET") {
            Epenpan.Daily <- Epenpan.Daily/1.078 * pan_coeff
        } else {
            Epenpan.Daily <- Epenpan.Daily/1.078
        }
    }

    ET.Daily <- Epenpan.Daily
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
    ET_formulation <- "PenPan"
    if (est == "potential ET") {
        ET_type <- "potential ET"
    } else if (est == "pan") {
        ET_type <- "Class-A Pan Evaporation"
    }
    Surface <- paste("user-defined, albedo =", alpha)
    if (solar == "data") {
        message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
    } else if (solar == "sunshine hours") {
        message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
    } else if (solar == "cloud") {
        message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
    } else {
        message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
    }
    results <- list(ET.Daily = ET.Daily, ET.Monthly = ET.Monthly,
        ET.Annual = ET.Annual, ET.MonthlyAve = ET.MonthlyAve,
        ET.AnnualAve = ET.AnnualAve, ET_formulation = ET_formulation,
        ET_type = ET_type, message1 = message1)

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
        if (ET_type == "potential ET") {
            message("Pan coeffcient: ", pan_coeff)
        }
        message(message1)


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

    if (save.csv == "yes") {
        # class(results) <- funname write to csv file
        for (i in 1:length(results)) {
            namer <- names(results[i])
            write.table(as.character(namer), file = "ET_PenPan.csv",
                dec = ".", quote = FALSE, col.names = FALSE,
                row.names = F, append = TRUE, sep = ",")
            write.table(data.frame(get(namer, results)), file = "ET_PenPan.csv",
                col.names = F, append = T, sep = ",")
        }
        invisible(results)
    } else {
        return(results)
    }
}