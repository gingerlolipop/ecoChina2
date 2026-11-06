# Standalone reference implementation and small logic checks.
# This file does not modify existing pipeline scripts or raster outputs.

dual_suitability_preserve_missing <- function(climate, soil = NULL, gate = 0.20) {
  if (length(gate) != 1L || !is.finite(gate) || gate < 0 || gate > 1) {
    stop("gate must be one finite probability in [0,1].")
  }
  raster_input <- inherits(climate, "SpatRaster")
  if (gate == 0) {
    if (raster_input) return(terra::ifel(is.finite(climate), climate, NA_real_))
    answer <- climate
    answer[!is.finite(answer)] <- NA_real_
    return(answer)
  }
  if (raster_input) {
    if (!inherits(soil, "SpatRaster")) stop("soil must also be a SpatRaster.")
    if (!terra::compareGeom(climate, soil, stopOnError = FALSE) ||
        terra::nlyr(climate) != terra::nlyr(soil)) {
      stop("Climate and soil rasters must have identical geometry and layer count.")
    }
    return(terra::ifel(!is.finite(climate) | !is.finite(soil), NA_real_,
                       terra::ifel(soil > gate, climate, 0)))
  }
  if (!is.numeric(climate) || !is.numeric(soil) || length(climate) != length(soil) ||
      !identical(dim(climate), dim(soil))) stop("Provide equally shaped numeric inputs.")
  answer <- climate
  available <- is.finite(climate) & is.finite(soil)
  answer[!available] <- NA_real_
  answer[available & soil <= gate] <- 0
  answer
}

candidate_availability_status <- function(dual, cutoff = 0.40) {
  if (length(cutoff) != 1L || !is.finite(cutoff)) stop("Invalid cutoff.")
  if (is.null(dim(dual))) dual <- matrix(dual, nrow = 1L)
  if (!is.numeric(dual) || !is.matrix(dual) || ncol(dual) == 0L || nrow(dual) == 0L) {
    stop("Use a nonempty numeric cells-by-ecotypes matrix.")
  }
  available <- is.finite(dual)
  n_available <- rowSums(available)
  has_analogue <- rowSums(available & dual >= cutoff) > 0L
  complete <- n_available == ncol(dual)
  # An observed suitable candidate establishes analogue presence. A claim of
  # no suitable candidate requires the complete set of required candidate data.
  status <- ifelse(has_analogue, "analogue_present",
    ifelse(complete, "no_current_analogue", "missing_candidate_data"))
  data.frame(n_available = n_available, n_required = ncol(dual),
             complete_candidate_data = complete, status = status)
}

local({
  c <- c(NA_real_, 0.8, 0.8, 0.8, Inf, 0.6)
  s <- c(0.1, 0.1, 0.3, NA_real_, 0.1, 0.2)
  fixed <- dual_suitability_preserve_missing(c, s)
  stopifnot(is.na(fixed[1]), fixed[2] == 0, fixed[3] == 0.8,
            is.na(fixed[4]), is.na(fixed[5]), fixed[6] == 0)
  climate_only <- dual_suitability_preserve_missing(c, gate = 0)
  stopifnot(climate_only[4] == 0.8, is.na(climate_only[1]), is.na(climate_only[5]))
  status <- candidate_availability_status(rbind(
    c(NA, NA), c(0, 0), c(0.6, NA), c(0.1, NA), c(0.2, 0.3), c(0.4, 0)))
  stopifnot(identical(status$status, c("missing_candidate_data", "no_current_analogue",
    "analogue_present", "missing_candidate_data", "no_current_analogue", "analogue_present")))
})
