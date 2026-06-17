library(terra)
library(doFuture)
library(foreach)

# ------------------------------------------------------------
# Assign each pixel to a vegetation zone using dual suitability maps.
# This version reads method-specific dual suitability folders:
#   dual suit/optimized_mf/normal
#   dual suit/optimized_rf/normal
#   dual suit/plain_mf/normal
#   dual suit/plain_rf/normal
# and the same structure for future periods.
# ------------------------------------------------------------

base_dir <- "H:/Jing/ecoChina2"
dual_dir <- file.path(base_dir, "dual suit")
output_root <- file.path(base_dir, "result maps")

r_file <- file.path(base_dir, "raster/ecosys_ori.tif")

zoneID <- c(1:7, 9:30, 31:50, 52:55)

threshold <- 0.1
novel_value <- 99
tie_tol <- 1e-4

set.seed(49)

method_order <- c("optimized_mf", "optimized_rf", "plain_mf", "plain_rf")
methods <- method_order[dir.exists(file.path(dual_dir, method_order))]

# Run only the first method if desired.
# methods <- "optimized_mf"

print(methods)

jobs <- do.call(rbind, lapply(methods, function(method) {
  method_dir <- file.path(dual_dir, method)
  scenarios <- list.dirs(method_dir, recursive = FALSE, full.names = FALSE)
  scenarios <- scenarios[sapply(scenarios, function(scenario) {
    any(file.exists(file.path(
      method_dir,
      scenario,
      paste0("dual_suitability_zone", zoneID, ".tif")
    )))
  })]
  data.frame(method = method, scenario = scenarios)
}))

print(jobs)

n_workers <- min(1, nrow(jobs))

registerDoFuture()
plan(multisession, workers = n_workers)

output_files <- foreach(
  j = seq_len(nrow(jobs)),
  .options.future = list(seed = TRUE)
) %dofuture% {
  
  library(terra)
  
  method <- jobs$method[j]
  scenario <- jobs$scenario[j]
  
  message("Processing: ", method, " | ", scenario)
  
  r <- rast(r_file)
  input_dir <- file.path(dual_dir, method, scenario)
  output_dir <- file.path(output_root, method)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  files <- file.path(input_dir, paste0("dual_suitability_zone", zoneID, ".tif"))
  names(files) <- zoneID
  files <- files[file.exists(files)]
  
  if (length(files) == 0) {
    warning("No dual suitability tif files found: ", method, " | ", scenario)
    return(NA_character_)
  }
  
  dual_stack <- rast(files)
  zoneID_available <- as.numeric(names(files))
  names(dual_stack) <- zoneID_available
  
  if (!compareGeom(dual_stack, r, stopOnError = FALSE)) {
    dual_stack <- resample(dual_stack, r, method = "near")
  }
  
  stack_with_r <- c(dual_stack, r)
  names(stack_with_r)[nlyr(stack_with_r)] <- "original_zone"
  
  is_normal <- scenario == "normal"
  
  assign_zone <- function(x) {
    current_values <- x[seq_along(zoneID_available)]
    original_zone  <- x[length(x)]
    
    if (is.na(original_zone)) return(NA_real_)
    if (all(is.na(current_values))) return(NA_real_)
    
    max_value <- max(current_values, na.rm = TRUE)
    
    if (!is_normal && max_value < threshold) {
      return(novel_value)
    }
    
    tied_zones <- zoneID_available[
      !is.na(current_values) & (max_value - current_values <= tie_tol)
    ]
    
    if (original_zone %in% tied_zones) {
      return(original_zone)
    }
    
    sample(tied_zones, 1)
  }
  
  output_file <- file.path(
    output_dir,
    paste0(
      "assigned_zone_", scenario,
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
    wopt = list(datatype = "INT2S", gdal = c("COMPRESS=LZW"))
  )
  
  message("Saved: ", output_file)
  output_file
}

print(output_files)
