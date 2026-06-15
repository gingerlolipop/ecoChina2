library(terra)
library(doFuture)
library(foreach)

# ------------------------------------------------------------
# Assign each pixel to a vegetation zone using dual suitability maps.
#
# Main rules:
#   1. If the original vegetation raster r is NA, return NA.
#      This preserves the original map mask.
#   2. The zone with the highest dual suitability is assigned.
#   3. For future scenarios, if the maximum suitability is below threshold,
#      assign 99. Here 99 means novel ecosystem / no suitable current zone.
#   4. For normal period, do not assign 99.
#      The normal period must be assigned to one of the existing zones.
#   5. If multiple zones are nearly tied:
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

zoneID <- c(1:7, 9:30, 31:50, 52:55)

threshold <- 0.1
novel_value <- 99
tie_tol <- 1e-4

set.seed(49)

scenarios <- list.dirs(
  dual_dir,
  recursive = FALSE,
  full.names = FALSE
)

scenarios <- scenarios[sapply(scenarios, function(scenario) {
  any(file.exists(file.path(
    dual_dir,
    scenario,
    paste0("dual_suitability_zone", zoneID, ".tif")
  )))
})]

print(scenarios)

# Parallelization is done by scenario, not by pixel.
# workers = 1 is more stable for large terra raster operations.
n_workers <- min(1, length(scenarios))

registerDoFuture()
plan(multisession, workers = n_workers)

output_files <- foreach(
  scenario = scenarios,
  .options.future = list(seed = TRUE)
) %dofuture% {
  
  library(terra)
  
  message("Processing scenario: ", scenario)
  
  r <- rast(r_file)
  
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
  
  # Use available zone IDs because some zone files may be missing.
  zoneID_available <- as.numeric(names(files))
  names(dual_stack) <- zoneID_available
  
  if (!compareGeom(dual_stack, r, stopOnError = FALSE)) {
    dual_stack <- resample(
      dual_stack,
      r,
      method = "near"
    )
  }
  
  stack_with_r <- c(dual_stack, r)
  names(stack_with_r)[nlyr(stack_with_r)] <- "original_zone"
  
  # Historical / normal period should not have novel ecosystem.
  is_normal <- scenario == "normal"
  
  assign_zone <- function(x) {
    
    current_values <- x[seq_along(zoneID_available)]
    original_zone  <- x[length(x)]
    
    # Preserve original map mask.
    if (is.na(original_zone)) {
      return(NA_real_)
    }
    
    if (all(is.na(current_values))) {
      return(NA_real_)
    }
    
    max_value <- max(current_values, na.rm = TRUE)
    
    # Novel ecosystem only applies to future scenarios.
    if (!is_normal && max_value < threshold) {
      return(novel_value)
    }
    
    # Near ties are treated as ties.
    # tie_tol = 1e-4 means differences after the 4th decimal place are ignored.
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
    paste0(
      "assigned_zone_",
      scenario,
      "_threshold", threshold,
      "_tol", tie_tol,
      "_novel99_maskNA_noNovelNormal.tif"
    )
  )
  
  app(
    stack_with_r,
    fun = assign_zone,
    filename = output_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "INT2S",
      gdal = c("COMPRESS=LZW")
    )
  )
  
  message("Saved: ", output_file)
  
  output_file
}

print(output_files)