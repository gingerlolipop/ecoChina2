library(terra)
library(doFuture)
library(foreach)

# ------------------------------------------------------------
# Assign each pixel to a vegetation zone using dual suitability maps.
#
# Main rules:
#   1. The zone with the highest dual suitability is assigned.
#   2. If the maximum suitability is below threshold, assign 99.
#      Here 99 means novel ecosystem / no suitable current vegetation zone.
#   3. If multiple zones have exactly the same highest suitability:
#        - keep the original vegetation zone if it is one of the tied zones;
#        - otherwise randomly choose one of the tied zones.
#
# The original vegetation raster r is also used as the spatial template,
# so output maps can be compared with r pixel by pixel.
# ------------------------------------------------------------

base_dir <- "H:/Jing/ecoChina2"
dual_dir <- file.path(base_dir, "dual suit")
output_dir <- file.path(base_dir, "result maps")

r_file <- file.path(base_dir, "raster/ecosys_ori.tif")
r <- rast(r_file)

zoneID <- c(1:7, 9:30, 31:50, 52:55)

threshold <- 0.1
novel_value <- 99
tie_tol <- 1e-6

set.seed(49)

scenarios <- list.dirs(
  dual_dir,
  recursive = FALSE,
  full.names = FALSE
)

# Keep only scenarios that already have at least one dual suitability tif.
# This is useful because some future scenarios may still be running.
# sapply() checks each scenario folder and returns TRUE/FALSE.
scenarios <- scenarios[sapply(scenarios, function(scenario) {
  any(file.exists(file.path(
    dual_dir,
    scenario,
    paste0("dual_suitability_zone", zoneID, ".tif")
  )))
})]

print(scenarios)

# Parallelization is done by scenario, not by pixel.
# Each worker processes one scenario folder independently.
# Do not set this too high because each scenario loads many raster layers.
n_workers <- min(1, length(scenarios))

registerDoFuture()
plan(multisession, workers = n_workers)

output_files <- foreach(
  scenario = scenarios,
  .options.future = list(
    seed = TRUE,
    packages = "terra"
  )
) %dofuture% {
  
  message("Processing scenario: ", scenario)
  
  input_dir <- file.path(dual_dir, scenario)
  
  files <- file.path(
    input_dir,
    paste0("dual_suitability_zone", zoneID, ".tif")
  )
  
  names(files) <- zoneID
  files <- files[file.exists(files)]
  
  if (length(files) == 0) {
    warning("No dual suitability tif files found for scenario: ", scenario)
    return(NA_character_)
  }
  
  dual_stack <- rast(files)
  
  # zoneID_available is used instead of zoneID because some zone files
  # may be missing. It keeps the layer-to-zone mapping correct.
  zoneID_available <- as.numeric(names(files))
  names(dual_stack) <- zoneID_available
  
  # The assigned map must match the original vegetation map r exactly
  # in extent, resolution, CRS, and number of rows/columns.
  # Otherwise later pixel-by-pixel comparison with r will be invalid.
  if (!compareGeom(dual_stack, r, stopOnError = FALSE)) {
    dual_stack <- resample(
      dual_stack,
      r,
      method = "near"
    )
  }
  
  # Add r as the last layer.
  # In assign_zone():
  #   x[1:length(zoneID_available)] = dual suitability values
  #   x[length(x)]                  = original vegetation zone from r
  stack_with_r <- c(dual_stack, r)
  names(stack_with_r)[nlyr(stack_with_r)] <- "original_zone"
  
  assign_zone <- function(x) {
    
    current_values <- x[seq_along(zoneID_available)]
    original_zone  <- x[length(x)]
    
    # Keep the original raster mask.
    # If the original vegetation map is NA, the assigned map should also be NA.
    if (is.na(original_zone)) {
      return(NA_real_)
    }
    
    if (all(is.na(current_values))) {
      return(NA_real_)
    }
    
    max_value <- max(current_values, na.rm = TRUE)
    
    if (max_value < threshold) {
      return(novel_value)
    }
    
    # Treat nearly equal values as ties.
    # tie_tol = 1e-6 means differences after the 6th decimal place are ignored.
    tied_zones <- zoneID_available[
      !is.na(current_values) & (max_value - current_values <= tie_tol)
    ]
    
    # If the original zone is among the near-best zones, keep it.
    if (original_zone %in% tied_zones) {
      return(original_zone)
    }
    
    sample(tied_zones, 1)
  }  
  output_file <- file.path(
    output_dir,
    paste0("assigned_zone_", scenario, "_threshold", threshold, "_novel99_originalTie.tif")
  )
  
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
  
  # terra::app() applies assign_zone() pixel by pixel across layers.
  # filename writes the result directly to disk, which avoids holding
  # the full output raster in memory.
  app(
    stack_with_r,
    fun = assign_zone,
    filename = output_file,
    overwrite = TRUE
  )
  
  message("Saved: ", output_file)
  
  output_file
}

print(output_files)