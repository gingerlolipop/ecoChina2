library(terra)

# Assign vegetation zones from dual suitability maps
# A map is generated only when all zone-level rasters exist.

base_dir <- "H:/Jing/ecoChina2"
dual_dir <- file.path(base_dir, "dual suit")
output_root <- file.path(base_dir, "result maps")
r_file <- file.path(base_dir, "raster/ecosys_ori.tif")

zoneID <- c(1:7, 9:30, 31:50, 52:55)

threshold <- 0.2
novel_value <- 99
tie_tol <- 1e-4

set.seed(49)

method_order <- c(
  "optimized_mf",
  "optimized_rf",
  "plain_mf",
  "plain_rf"
)

future_order <- c(
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

scenario_order <- c("normal", future_order)

# Build jobs only from complete scenario folders.
jobs <- list()

for (method in method_order) {
  for (scenario in scenario_order) {
    input_dir <- file.path(dual_dir, method, scenario)
    if (!dir.exists(input_dir)) next
    
    files <- file.path(
      input_dir,
      paste0("dual_suitability_zone", zoneID, ".tif")
    )
    
    missing_zone <- zoneID[!file.exists(files)]
    
    if (length(missing_zone) > 0) {
      cat(
        "[SKIP INCOMPLETE]", method, "|", scenario,
        "| missing zones:", paste(missing_zone, collapse = ", "), "\n"
      )
      next
    }
    
    jobs[[length(jobs) + 1L]] <- data.frame(
      method = method,
      scenario = scenario
    )
  }
}

if (length(jobs) == 0) {
  stop("No complete method-scenario folders were found.")
}

jobs <- do.call(rbind, jobs)
print(jobs)

r <- rast(r_file)

for (j in seq_len(nrow(jobs))) {
  method <- jobs$method[j]
  scenario <- jobs$scenario[j]
  
  cat("\n==============================\n")
  cat("ASSIGN:", method, "|", scenario, "\n")
  cat("==============================\n")
  
  input_dir <- file.path(dual_dir, method, scenario)
  output_dir <- file.path(output_root, method)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  files <- file.path(
    input_dir,
    paste0("dual_suitability_zone", zoneID, ".tif")
  )
  names(files) <- zoneID
  
  dual_stack <- rast(files)
  names(dual_stack) <- zoneID
  
  if (!compareGeom(dual_stack, r, stopOnError = FALSE)) {
    cat("[RESAMPLE] dual suitability rasters to original map geometry\n")
    dual_stack <- resample(dual_stack, r, method = "bilinear")
    names(dual_stack) <- zoneID
  }
  
  stack_with_r <- c(dual_stack, r)
  names(stack_with_r)[nlyr(stack_with_r)] <- "original_zone"
  
  is_normal <- scenario == "normal"
  
  assign_zone <- function(x) {
    suit <- x[seq_along(zoneID)]
    original_zone <- x[length(x)]
    
    if (is.na(original_zone) || original_zone == 8) {
      return(NA_real_)
    }
    
    if (all(is.na(suit))) {
      return(NA_real_)
    }
    
    max_suit <- max(suit, na.rm = TRUE)
    
    if (!is_normal && max_suit < threshold) {
      return(novel_value)
    }
    
    tied_zone <- zoneID[
      !is.na(suit) & (max_suit - suit <= tie_tol)
    ]
    
    if (length(tied_zone) == 0) {
      return(NA_real_)
    }
    
    if (original_zone %in% tied_zone) {
      return(original_zone)
    }
    
    tied_zone[sample.int(length(tied_zone), 1L)]
  }
  
  output_file <- file.path(
    output_dir,
    paste0(
      "assigned_zone_", scenario,
      "_threshold", threshold,
      "_tol", tie_tol,
      "_novel99_maskNA8_noNovelNormal.tif"
    )
  )
  
  assigned_raw <- app(
    stack_with_r,
    fun = assign_zone
  )
  
  assigned <- ifel(
    is.na(r) | r == 8,
    NA,
    assigned_raw
  )
  
  allowed_value <- if (is_normal) {
    zoneID
  } else {
    c(zoneID, novel_value)
  }
  
  assigned_freq <- freq(assigned)
  bad_value <- setdiff(assigned_freq$value, allowed_value)
  
  if (length(bad_value) > 0) {
    stop(
      "Assigned raster contains unexpected values: ",
      paste(bad_value, collapse = ", ")
    )
  }
  
  writeRaster(
    assigned,
    output_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "INT2S",
      gdal = "COMPRESS=LZW"
    )
  )
  
  assigned_check <- rast(output_file)
  saved_freq <- freq(assigned_check)
  bad_value <- setdiff(saved_freq$value, allowed_value)
  
  if (length(bad_value) > 0) {
    stop(
      "Saved raster contains unexpected values: ",
      paste(bad_value, collapse = ", ")
    )
  }
  
  cat("[SAVED]", output_file, "\n")
  
  rm(
    dual_stack,
    stack_with_r,
    assigned_raw,
    assigned,
    assigned_check
  )
  gc()
}