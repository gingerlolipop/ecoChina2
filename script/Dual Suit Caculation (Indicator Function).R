# Define the dual suitability function
# Formula: Suitability_dual  =I(Suitability_soil>threshold)* Suitability_climate+ I(Suitability_soil<threshold)  0
# The indicator function I() can be represented using the as.integer() function
# in R which will convert TRUE to 1 and FALSE to 0.

# Define the dual suitability function
calculate_dual_suitability <- function(climate, soil) {
  ifel(
    is.na(soil),
    NA,
    ifel(soil > threshold, climate, 0)
  )
}

# Exaple threshold == 0.2
threshold <- 0.2 

# Apply the function
dual_suitability <- overlay(pclim75, psoil75_masked, fun = calculate_dual_suitability)

# Plot and save the result
plot(dual_suitability)
writeRaster(dual_suitability, "dual_suitability.tif", overwrite = TRUE)


