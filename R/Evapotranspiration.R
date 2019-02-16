
ET <- function(data, constants, ...) UseMethod("ET")

ET.default <- function(data, constants, crop = NULL, alpha = NULL,
    solar = NULL, wind = NULL, ...) {

    if (is.null(solar)) {
        if (!is.null(data$Rs)) {
            solar = "data"
        } else if (!is.null(data$n)) {
            solar = "sunshine hours"
        } else if (!is.null(data$Cd)) {
            solar = "cloud"
        } else if (!is.null(data$Precip)) {
            solar = "monthly precipitation"
        }
    }

    if (is.null(wind)) {
        if (!is.null(data$u2) | !is.null(data$uz)) {
            wind = "yes"
        } else {
            wind = "no"
        }
    }
    if (all(any(is.null(data$RHmax), is.null(data$RHmin)), is.null(data$Rs),
        is.null(data$n), is.null(data$Cd), is.null(data$Precip),
        all(is.null(data$uz), is.null(data$u2)), any(is.null(data$Tmax),
            is.null(data$Tmin)))) {
        # no data available
        stop("No ET model is suitable according to the data availability")

    } else if (all(any(is.null(data$RHmax), is.null(data$RHmin)),
        is.null(data$Rs), is.null(data$n), is.null(data$Cd),
        is.null(data$Precip)) & all(is.null(data$uz), is.null(data$u2)) &
        all(!is.null(data$Tmax), !is.null(data$Tmin))) {
        # Only Tmax/Tmin available
        message("No ET model specified, choose the Hargreaves-Samani model according to the data availability")
        class(data) = "HargreavesSamani"
        ET(data, constants, ...)
    } else if (all(any(is.null(data$RHmax), is.null(data$RHmin)),
        all(is.null(data$uz), is.null(data$u2))) & any(!is.null(data$Rs),
        !is.null(data$n), !is.null(data$Cd), !is.null(data$Precip)) &
        all(!is.null(data$Tmax), !is.null(data$Tmin))) {
        # Tmax/Tmin & any Rs data available
        message("No ET model specified, choose the Makkink model according to the data availability")
        class(data) = "Makkink"

        ET(data, constants, solar = solar, ...)
    } else if (all(is.null(data$uz), is.null(data$u2)) & any(!is.null(data$Rs),
        !is.null(data$n), !is.null(data$Cd), !is.null(data$Precip)) &
        all(!is.null(data$Tmax), !is.null(data$Tmin), !is.null(data$RHmax),
            !is.null(data$RHmin)) & !is.null(alpha)) {
        # Tmax/Tmin & any Rs & RHmax/RHmin data available
        message("No ET model specified, choose the Priestley-Taylor model according to the data availability")
        class(data) = "PriestleyTaylor"

        ET(data, constants, solar = solar, ...)
    } else if (any(!is.null(data$Rs), !is.null(data$n), !is.null(data$Cd),
        !is.null(data$Precip)) & all(!is.null(data$Tmax), !is.null(data$Tmin),
        !is.null(data$RHmax), !is.null(data$RHmin)) & any(!is.null(data$uz),
        !is.null(data$u2))) {
        # All data available & crop specified
        Flag = 1
    } else {
        "No ET model can be recommended according to the data availability"
    }

    if (exists("Flag")) {
        if (Flag == 1) {
            # Penman-Monteith or Penman
            if (is.null(crop) | all(crop != "short", crop !=
                "tall")) {
                alpha = 0.08
                z0 = 0.001
                message("No ET model specified and no valid evaporative surface specified, choose the Penman model for open-water evaporation according to the data availability")
                class(data) = "Penman"

                ET(data, constants, solar = solar, wind = wind,
                  windfunction_ver = 1948, alpha = alpha, z0 = z0,
                  ...)
            } else {
                message("No ET model specified, choose the Penman-Monteith model according to the data availability")
                class(data) = "PenmanMonteith"
                ET(data, constants, solar = solar, wind = wind,
                  crop = crop, ...)
            }
        }
    }
}
