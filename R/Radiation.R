# Calculate radiation variables
Radiation <- function(data, constants, ts = "monthly", solar = "sunshine hours",
    Tdew = TRUE, alpha = NULL, model = "CRAE") {
    if (model != "CRWE" & model != "CRAE" & model != "CRLE") {
        stop("Error: please choose either 'CRWE' or 'CRAE' for the model to use")
    }
    if (ts == "daily") {
        stop("Error: Morton models are not available for daily time step")
    }
    if (is.null(data$Tmax)) {
        stop("Required data missing for 'Tmax' or 'Temp'")
    }
    if (is.null(data$Tmin)) {
        stop("Required data missing for 'Tmin' or 'Temp'")
    }
    if (Tdew == TRUE & is.null(data$Tdew)) {
        stop("Required data missing for 'Tdew'")
    }
    if (Tdew == FALSE) {
        if (is.null(data$va)) {
            if (is.null(data$RHmax) | is.null(data$RHmin)) {
                stop("Required data missing for 'va' , or 'RHmax' and 'RHmin' (or 'RH')")
            }
        }
    }

    if (solar == "sunshine hours" & is.null(data$n)) {
        stop("Required data missing for 'n'")
    }
    if (solar == "monthly precipitation") {
        stop("Only 'data', 'sunshine hours' and 'cloud' are accepted because estimations of sunshine hours is required")
    }
    if (solar == "sunshine hours" & is.null(data$Precip)) {
        if ("PA" %in% names(constants) == FALSE) {
            stop("Required data missing for 'Precip.daily' or required constant missing for 'PA'")
        }
    }
    T_Mo.temp <- (data$Tmax + data$Tmin)/2
    T_Mo <- aggregate(T_Mo.temp, as.yearmon(data$Date.daily,
        "%m/%y"), FUN = mean)
    # calculate Tdew
    if (Tdew == TRUE) {
        Tdew_Mo <- aggregate(data$Tdew, as.yearmon(data$Date.daily,
            "%d/%m/%y"), FUN = mean)
    } else if (!is.null(data$va)) {
        vabar_Mo <- aggregate(data$va, as.yearmon(data$Date.daily,
            "%m/%y"), FUN = mean)
        Tdew_Mo <- (116.9 + 237.3 * log(vabar_Mo))/(16.78 - log(vabar_Mo))
    } else {
        vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2
        vabar_Mo <- aggregate(vabar, as.yearmon(data$Date.daily,
            "%m/%y"), FUN = mean)
        Tdew_Mo <- (116.9 + 237.3 * log(vabar_Mo))/(16.78 - log(vabar_Mo))
    }

    # calculate vas
    if (is.null(data$vs)) {
        vas <- data$vs
    } else {
        vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax/(data$Tmax +
            237.3))
        vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin/(data$Tmin +
            237.3))
        vas <- (vs_Tmax + vs_Tmin)/2
    }


    delta <- 4098 * (0.6108 * exp(17.27 * T_Mo/(T_Mo + 237.3)))/(T_Mo +
        237.3)^2
    deltas <- 0.409 * sin(2 * pi/365 * data$J - 1.39)
    omegas <- acos(-tan(constants$lat_rad) * tan(deltas))
    if (!is.null(data$Precip)) {
        PA <- mean(aggregate(data$Precip, floor(as.numeric(as.yearmon(data$Date.daily,
            "%m/%y"))), FUN = sum))
    } else {
        PA <- constants$PA
    }
    if (model == "CRLE" | model == "CRWE") {
        constants$epsilonMo <- 0.97
        constants$fz <- 25
        constants$b0 <- 1.12
        constants$b1 <- 13
        constants$b2 <- 1.12
    }
    ptops <- ((288 - 0.0065 * constants$Elev)/288)^5.256
    alpha_zd <- 0.26 - 0.00012 * PA * sqrt(ptops) * (1 + abs(constants$lat/42) +
        (constants$lat/42)^2)
    if (alpha_zd < 0.11) {
        alpha_zd <- 0.11
    } else if (alpha_zd > 0.17) {
        alpha_zd <- 0.17
    }

    if (solar == "sunshine hours") {
        N <- 24/pi * omegas
        S_daily <- data$n/N
        for (i in 1:length(S_daily)) {
            if (S_daily[i] > 1) {
                S_daily[i] <- 1
            }
        }
        S <- mean(S_daily)


        vD_Mo <- 6.11 * exp(constants$alphaMo * Tdew_Mo/(Tdew_Mo +
            constants$betaMo))
        v_Mo <- 6.11 * exp(constants$alphaMo * T_Mo/(T_Mo + constants$betaMo))
        deltaMo <- constants$alphaMo * constants$betaMo * v_Mo/((T_Mo +
            constants$betaMo)^2)
        thetaMo <- (23.2 * sin((29.5 * data$i - 94) * pi/180)) *
            pi/180
        Z_Mo <- acos(cos(constants$lat_rad - thetaMo))
        for (i in 1:length(Z_Mo)) {
            if (cos(Z_Mo[i]) < 0.001) {
                Z_Mo[i] <- acos(0.001)
            }
        }
        omegaMo <- acos(1 - cos(Z_Mo)/(cos(constants$lat_rad) *
            cos(thetaMo)))
        cosz <- cos(Z_Mo) + (sin(omegaMo)/omegaMo - 1) * cos(constants$lat_rad) *
            cos(thetaMo)
        etaMo <- 1 + 1/60 * sin((29.5 * data$i - 106) * pi/180)
        G_E <- 1354/(etaMo^2) * omegaMo/pi * cosz
        alpha_zz <- matrix(NA, length(v_Mo), 1)

        if (model == "CRLE" | model == "CRWE") {
            alpha_zz[1:length(v_Mo)] <- 0.05
        } else {
            alpha_zz[1:length(v_Mo)] <- alpha_zd
            for (i in 1:length(v_Mo)) {
                if (alpha_zz[i] < 0.11) {
                  alpha_zz[i] <- 0.11
                } else {
                  if (alpha_zz[i] > 0.5 * (0.91 - vD_Mo[i]/v_Mo[i])) {
                    alpha_zz[i] <- 0.91 - vD_Mo[i]/v_Mo[i]
                  } else {
                  }
                }
            }
        }
        c_0 <- as.vector(v_Mo - vD_Mo)
        for (i in 1:length(c_0)) {
            if (c_0[i] < 0) {
                c_0[i] <- 0
            } else {
                if (c_0[i] > 1) {
                  c_0[i] <- 1
                } else {
                  c_0[i] <- c_0[i]
                }
            }
        }
        alpha_z <- alpha_zz + (1 - c_0^2) * (0.34 - alpha_zz)
        alpha_0 <- alpha_z * (exp(1.08) - ((2.16 * cos(Z_Mo))/pi +
            sin(Z_Mo)) * exp(0.012 * Z_Mo * 180/pi))/(1.473 *
            (1 - sin(Z_Mo)))
        W_Mo <- vD_Mo/(0.49 + T_Mo/129)
        c_1 <- as.vector(21 - T_Mo)
        for (i in 1:length(c_1)) {
            if (c_1[i] < 0) {
                c_1[i] <- 0
            } else {
                if (c_1[i] > 5) {
                  c_1[i] <- 5
                } else {
                  c_1[i] <- c_1[i]
                }
            }
        }
        j_Mo <- (0.5 + 2.5 * (cosz)^2) * exp(c_1 * (ptops - 1))
        tauMo <- exp(-0.089 * (ptops * 1/cosz)^0.75 - 0.083 *
            (j_Mo/cosz)^0.9 - 0.029 * (W_Mo/cosz)^0.6)
        tauaMo <- as.vector(exp(-0.0415 * (j_Mo/cosz)^0.9 - (0.0029)^0.5 *
            (W_Mo/cosz)^0.3))
        for (i in 1:length(tauaMo)) {
            if (tauaMo[i] < exp(-0.0415 * (as.matrix(j_Mo/cosz)[i])^0.9 -
                0.029 * (as.matrix(W_Mo/cosz)[i])^0.6)) {
                tauaMo[i] <- exp(-0.0415 * (as.matrix(j_Mo/cosz)[i])^0.9 -
                  0.029 * (as.matrix(W_Mo/cosz)[i])^0.6)
            } else {
                tauaMo[i] <- tauaMo[i]
            }
        }
        G_0 <- G_E * tauMo * (1 + (1 - tauMo/tauaMo) * (1 + alpha_0 *
            tauMo))
        G_Mo <- S * G_0 + (0.08 + 0.3 * S) * (1 - S) * G_E
        alpha_Mo <- alpha_0 * (S + (1 - S) * (1 - Z_Mo/330 *
            180/pi))
        c_2 <- as.vector(10 * (vD_Mo/v_Mo - S - 0.42))
        for (i in 1:length(c_2)) {
            if (c_2[i] < 0) {
                c_2[i] <- 0
            } else {
                if (c_2[i] > 1) {
                  c_2[i] <- 1
                } else {
                  c_2[i] <- c_2[i]
                }
            }
        }
        rouMo <- 0.18 * ((1 - c_2) * (1 - S)^2 + c_2 * (1 - S)^0.5) *
            1/ptops
        B_Mo <- as.vector(constants$epsilonMo * constants$sigmaMo *
            (T_Mo + 273)^4 * (1 - (0.71 + 0.007 * vD_Mo * ptops) *
            (1 + rouMo)))
        for (i in 1:length(B_Mo)) {
            if (B_Mo[i] < 0.05 * constants$epsilonMo * constants$sigmaMo *
                (T_Mo[i] + 274)^4) {
                B_Mo[i] <- 0.05 * constants$epsilonMo * constants$sigmaMo *
                  (T_Mo[i] + 274)^4
            } else {
                B_Mo[i] <- B_Mo[i]
            }
        }
        R_T <- (1 - alpha_Mo) * G_Mo - B_Mo
    } else if (solar == "data") {
        vD_Mo <- 6.11 * exp(constants$alphaMo * Tdew_Mo/(Tdew_Mo +
            constants$betaMo))
        v_Mo <- 6.11 * exp(constants$alphaMo * T_Mo/(T_Mo + constants$betaMo))
        ptops <- ((288 - 0.0065 * constants$Elev)/288)^5.256
        deltaMo <- constants$alphaMo * constants$betaMo * v_Mo/((T_Mo +
            constants$betaMo)^2)
        thetaMo <- (23.2 * sin((29.5 * data$i - 94) * pi/180)) *
            pi/180
        Z_Mo <- acos(cos(constants$lat_rad - thetaMo))
        for (i in 1:length(Z_Mo)) {
            if (cos(Z_Mo[i]) < 0.001) {
                Z_Mo[i] <- acos(0.001)
            }
        }
        omegaMo <- acos(1 - cos(Z_Mo)/(cos(constants$lat_rad) *
            cos(thetaMo)))
        cosz <- cos(Z_Mo) + (sin(omegaMo)/omegaMo - 1) * cos(constants$lat_rad) *
            cos(thetaMo)
        etaMo <- 1 + 1/60 * sin((29.5 * data$i - 106) * pi/180)
        G_E <- 1354/(etaMo^2) * omegaMo/pi * cosz
        G_Mo <- aggregate(data$Rs * 10^6/86400, as.yearmon(data$Date.daily,
            "%m/%y"), FUN = mean)
        # S = NULL
        Ta <- (data$Tmax + data$Tmin)/2
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

        alpha_zz <- matrix(NA, length(v_Mo), 1)
        if (model == "CRLE" | model == "CRWE") {
            alpha_zz[1:length(v_Mo)] <- 0.05
        } else {
            alpha_zz[1:length(v_Mo)] <- alpha_zd
            for (i in 1:length(v_Mo)) {
                if (alpha_zz[i] < 0.11) {
                  alpha_zz[i] <- 0.11
                } else {
                  if (alpha_zz[i] > 0.5 * (0.91 - vD_Mo[i]/v_Mo[i])) {
                    alpha_zz[i] <- 0.5 * (0.91 - vD_Mo[i]/v_Mo[i])
                  } else {
                  }
                }
            }
        }
        c_0 <- as.vector(v_Mo - vD_Mo)
        for (i in 1:length(c_0)) {
            if (c_0[i] < 0) {
                c_0[i] <- 0
            } else {
                if (c_0[i] > 1) {
                  c_0[i] <- 1
                } else {
                  c_0[i] <- c_0[i]
                }
            }
        }
        alpha_z <- alpha_zz + (1 - c_0^2) * (0.34 - alpha_zz)
        alpha_0 <- alpha_z * (exp(1.08) - ((2.16 * cos(Z_Mo))/pi +
            sin(Z_Mo)) * exp(0.012 * Z_Mo * 180/pi))/(1.473 *
            (1 - sin(Z_Mo)))
        W_Mo <- vD_Mo/(0.49 + T_Mo/129)
        c_1 <- as.vector(21 - T_Mo)
        for (i in 1:length(c_1)) {
            if (c_1[i] < 0) {
                c_1[i] <- 0
            } else {
                if (c_1[i] > 5) {
                  c_1[i] <- 5
                } else {
                  c_1[i] <- c_1[i]
                }
            }
        }
        j_Mo <- (0.5 + 2.5 * (cosz)^2) * exp(c_1 * (ptops - 1))
        tauMo <- exp(-0.089 * (ptops * 1/cosz)^0.75 - 0.083 *
            (j_Mo/cosz)^0.9 - 0.029 * (W_Mo/cosz)^0.6)
        tauaMo <- as.vector(exp(-0.0415 * (j_Mo/cosz)^0.9 - (0.0029)^0.5 *
            (W_Mo/cosz)^0.3))
        for (i in 1:length(tauaMo)) {
            if (tauaMo[i] < exp(-0.0415 * (as.matrix(j_Mo/cosz)[i])^0.9 -
                0.029 * (as.matrix(W_Mo/cosz)[i])^0.6)) {
                tauaMo[i] <- exp(-0.0415 * (as.matrix(j_Mo/cosz)[i])^0.9 -
                  0.029 * (as.matrix(W_Mo/cosz)[i])^0.6)
            } else {
                tauaMo[i] <- tauaMo[i]
            }
        }
        G_0 <- G_E * tauMo * (1 + (1 - tauMo/tauaMo) * (1 + alpha_0 *
            tauMo))


        S <- 0.53 * G_Mo/(G_0 - 0.47 * G_Mo)
        S[which(S < 0)] <- 0
        S[which(S > 1)] <- 1

        alpha_Mo <- alpha_0 * (S + (1 - S) * (1 - Z_Mo/330 *
            180/pi))
        c_2 <- as.vector(10 * (vD_Mo/v_Mo - S - 0.42))
        for (i in 1:length(c_2)) {
            if (c_2[i] < 0) {
                c_2[i] <- 0
            } else {
                if (c_2[i] > 1) {
                  c_2[i] <- 1
                } else {
                  c_2[i] <- c_2[i]
                }
            }
        }
        rouMo <- 0.18 * ((1 - c_2) * (1 - S)^2 + c_2 * (1 - S)^0.5) *
            1/ptops
        B_Mo <- as.vector(constants$epsilonMo * constants$sigmaMo *
            (T_Mo + 273)^4 * (1 - (0.71 + 0.007 * vD_Mo * ptops) *
            (1 + rouMo)))
        for (i in 1:length(B_Mo)) {
            if (B_Mo[i] < 0.05 * constants$epsilonMo * constants$sigmaMo *
                (T_Mo[i] + 274)^4) {
                B_Mo[i] <- 0.05 * constants$epsilonMo * constants$sigmaMo *
                  (T_Mo[i] + 274)^4
            } else {
                B_Mo[i] <- B_Mo[i]
            }
        }
        R_T <- (1 - alpha_Mo) * G_Mo - B_Mo
    }
    if (solar == "sunshine hours") {
        message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
    } else if (solar == "cloud") {
        message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
    } else {
        message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
    }
    if (Tdew == TRUE) {
        message6 <- "Data of dew point temperature has been used"
    } else {
        message6 <- "Data of average vapour pressure has been used to estimate dew point pressure"
    }
    variables <- list(T_Mo = T_Mo, Tdew_Mo = Tdew_Mo, S = S,
        R_T = R_T, ptops = ptops, vD_Mo = vD_Mo, v_Mo = v_Mo,
        deltaMo = deltaMo, G_E = G_E, G_Mo = G_Mo, alpha_Mo = alpha_Mo,
        B_Mo = B_Mo, message1 = message1, message6 = message6)
    return(variables)
}