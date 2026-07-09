# 8. Rank dual suitability by pixel
# ============================================================
# For each method x scenario:
#   1. Read all zone-level dual suitability rasters.
#   2. For each pixel, rank zones by dual suitability.
#   3. Save top-K ranked zones and their dual suitability values.
#
# Main outputs:
#   dual suit ranking/{method}/{scenario}/ranked_zone.tif
#   dual suit ranking/{method}/{scenario}/ranked_suitability.tif
#   dual suit ranking/{method}/{scenario}/ranked_summary.tif
#
# Resume behavior:
#   - If all three output rasters already exist, have correct geometry,
#     and have the expected number of layers, the job is reused.
#   - Existing valid outputs are not recalculated when the script is rerun.
#   - Missing/incomplete jobs are skipped and recorded, not allowed to stop
#     the whole script.
#
# Notes:
#   - Zone 8 and zone 51 were not modeled.
#   - NA cells and original zone 8 cells are masked out.
#   - Only positive dual suitability values are ranked by default.
#   - Temporary raster filenames are unique to avoid overwrite errors.

library(terra)
library(data.table)

rm(list = ls())
gc()

# 0. Paths and settings =========================================================

base_dir <- "H:/Jing/ecoChina2"

dual_root <- file.path(base_dir, "dual suit")
output_root <- file.path(base_dir, "dual suit ranking")
table_dir <- file.path(output_root, "tables")
tmp_dir <- file.path(base_dir, "tmp_dual_suit_ranking")

reference_file <- file.path(base_dir, "raster/ecosys_ori.tif")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

terraOptions(tempdir = tmp_dir, memfrac = 0.15)

# Zones that were modeled. Zone 8 and zone 51 are excluded.
zoneID <- c(1:7, 9:50, 52:55)

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

# Run all by default. Restrict these only for testing.
methods_to_run <- method_order
scenarios_to_run <- scenario_order

# Threshold used for summary only.
dual_threshold <- 0.2

# Only rank values above this.
# 0 means dual suitability = 0 is ignored.
rank_min_suitability <- 0

# Keep only top K zones per pixel.
top_k_rank <- 10

# terra::app parallelization can be unstable on Windows.
# If worker errors occur, set this to 1L.
rank_cores <- min(4L, max(1L, parallel::detectCores() - 1L))

# Resume control.
# TRUE = do not rerun jobs whose output rasters are already valid.
reuse_existing_outputs <- TRUE

# If a job is invalid or incomplete, old outputs are overwritten before recomputing.
overwrite_outputs <- TRUE

# 1. Helper functions ===========================================================

cat0 <- function(...) {
  cat(..., "\n", sep = "")
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "unnamed", x)
}

regex_escape <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x, perl = TRUE)
}

remove_raster_files <- function(filepath, must_remove = FALSE) {
  base <- tools::file_path_sans_ext(filepath)
  folder <- dirname(filepath)
  base_nm <- basename(base)
  
  side_files <- character(0)
  
  if (dir.exists(folder)) {
    side_files <- list.files(
      folder,
      pattern = paste0("^", regex_escape(base_nm), "\\."),
      full.names = TRUE
    )
  }
  
  files <- unique(c(
    filepath,
    paste0(filepath, ".aux.xml"),
    paste0(filepath, ".ovr"),
    paste0(filepath, ".msk"),
    paste0(base, ".aux.xml"),
    paste0(base, ".ovr"),
    paste0(base, ".tfw"),
    side_files
  ))
  
  files <- files[file.exists(files)]
  
  if (length(files) > 0) {
    gc()
    suppressWarnings(unlink(files, force = TRUE, recursive = TRUE))
    Sys.sleep(0.2)
    gc()
  }
  
  still_exists <- files[file.exists(files)]
  
  if (must_remove && length(still_exists) > 0) {
    stop(
      "Cannot remove existing raster file(s):\n",
      paste(still_exists, collapse = "\n"),
      "\n\nClose these files in QGIS, ArcGIS, RStudio Viewer, or another R session, then rerun."
    )
  }
  
  invisible(TRUE)
}

dual_file <- function(method, scenario, zone) {
  file.path(
    dual_root,
    method,
    scenario,
    paste0("dual_suitability_zone", zone, ".tif")
  )
}

rank_layer_names <- function(k, suffix) {
  paste0("rank", seq_len(k), "_", suffix)
}

rank_output_files <- function(method, scenario) {
  out_dir <- file.path(output_root, method, scenario)
  list(
    out_dir = out_dir,
    ranked_zone_file = file.path(out_dir, "ranked_zone.tif"),
    ranked_suit_file = file.path(out_dir, "ranked_suitability.tif"),
    ranked_summary_file = file.path(out_dir, "ranked_summary.tif")
  )
}

valid_output <- function(file, expected_nlyr, template) {
  if (!file.exists(file)) return(FALSE)
  
  x <- tryCatch(
    rast(file),
    error = function(e) NULL
  )
  
  if (is.null(x)) return(FALSE)
  if (nlyr(x) != expected_nlyr) return(FALSE)
  
  isTRUE(compareGeom(x, template, stopOnError = FALSE))
}

valid_rank_outputs <- function(method, scenario, template, k_rank, n_summary_layers) {
  out <- rank_output_files(method, scenario)
  
  valid_zone <- valid_output(out$ranked_zone_file, k_rank, template)
  valid_suit <- valid_output(out$ranked_suit_file, k_rank, template)
  valid_summary <- valid_output(out$ranked_summary_file, n_summary_layers, template)
  
  list(
    all_valid = valid_zone && valid_suit && valid_summary,
    valid_zone = valid_zone,
    valid_suit = valid_suit,
    valid_summary = valid_summary,
    ranked_zone_file = out$ranked_zone_file,
    ranked_suitability_file = out$ranked_suit_file,
    ranked_summary_file = out$ranked_summary_file
  )
}

make_layer_index <- function(method, scenario, k_rank, summary_layer_names) {
  rbindlist(list(
    data.table(
      method = method,
      scenario = scenario,
      raster_type = "ranked_zone",
      layer = seq_len(k_rank),
      layer_name = rank_layer_names(k_rank, "zone"),
      rank = seq_len(k_rank),
      meaning = "zone ID at this rank"
    ),
    data.table(
      method = method,
      scenario = scenario,
      raster_type = "ranked_suitability",
      layer = seq_len(k_rank),
      layer_name = rank_layer_names(k_rank, "suit"),
      rank = seq_len(k_rank),
      meaning = "dual suitability at this rank"
    ),
    data.table(
      method = method,
      scenario = scenario,
      raster_type = "ranked_summary",
      layer = seq_along(summary_layer_names),
      layer_name = summary_layer_names,
      rank = NA_integer_,
      meaning = c(
        "number of zones with dual suitability above rank_min_suitability",
        "number of zones with dual suitability above dual_threshold",
        "rank1_suit minus rank2_suit",
        "highest dual suitability",
        "second-highest dual suitability",
        "1 if no zone is above dual_threshold, otherwise 0"
      )
    )
  ))
}

read_dual_stack <- function(files, template, valid_mask) {
  s <- rast(files)
  names(s) <- as.character(zoneID)
  
  if (!compareGeom(s, template, stopOnError = FALSE)) {
    cat0("[RESAMPLE] dual suitability rasters to reference geometry")
    s <- resample(s, template, method = "bilinear")
    names(s) <- as.character(zoneID)
  }
  
  # Keep only valid original-map cells.
  s <- mask(s, valid_mask)
  names(s) <- as.character(zoneID)
  
  s
}

# Function factory.
# The returned function is self-contained, so terra workers do not need to find
# external helper functions.
make_rank_fun <- function(zoneID, rank_min, threshold, top_k) {
  force(zoneID)
  force(rank_min)
  force(threshold)
  force(top_k)
  
  function(x) {
    k <- min(top_k, length(zoneID))
    
    out_zone <- rep(NA_real_, k)
    out_suit <- rep(NA_real_, k)
    
    valid_rank <- which(
      !is.na(x) &
        is.finite(x) &
        x > rank_min
    )
    
    valid_suit <- which(
      !is.na(x) &
        is.finite(x) &
        x > threshold
    )
    
    if (length(valid_rank) > 0) {
      # Higher dual suitability first.
      # Exact ties are ordered by smaller zone ID for reproducibility.
      ord <- valid_rank[order(-x[valid_rank], zoneID[valid_rank])]
      ord <- ord[seq_len(min(length(ord), k))]
      
      out_zone[seq_along(ord)] <- zoneID[ord]
      out_suit[seq_along(ord)] <- x[ord]
    }
    
    top1_suit <- if (!is.na(out_suit[1])) out_suit[1] else NA_real_
    
    top2_suit <- if (
      k >= 2 &&
      !is.na(out_suit[2])
    ) {
      out_suit[2]
    } else {
      NA_real_
    }
    
    top1_minus_top2 <- if (
      is.finite(top1_suit) &&
      is.finite(top2_suit)
    ) {
      top1_suit - top2_suit
    } else {
      NA_real_
    }
    
    novel_by_threshold <- as.numeric(length(valid_suit) == 0)
    
    c(
      out_zone,
      out_suit,
      n_zone_ranked = length(valid_rank),
      n_zone_above_threshold = length(valid_suit),
      top1_minus_top2 = top1_minus_top2,
      top1_suit = top1_suit,
      top2_suit = top2_suit,
      novel_by_threshold = novel_by_threshold
    )
  }
}

# 2. Reference map ==============================================================

if (!file.exists(reference_file)) {
  stop("Missing reference raster: ", reference_file)
}

r <- rast(reference_file)
names(r) <- "original_zone"

# Mask out NA and unmodelled original zone 8.
valid_mask <- ifel(
  !is.na(r) & r != 8,
  1,
  NA
)

k_rank <- min(top_k_rank, length(zoneID))
n_summary_layers <- 6L

expected_zone_layers <- k_rank
expected_suit_layers <- k_rank
expected_summary_layers <- n_summary_layers

summary_layer_names <- c(
  "n_zone_ranked",
  "n_zone_above_threshold",
  "top1_minus_top2",
  "top1_suit",
  "top2_suit",
  "novel_by_threshold"
)

cat(
  "\n[RANKING SETTINGS]\n",
  "Zones: ", length(zoneID), "\n",
  "Top K: ", k_rank, "\n",
  "Rank minimum suitability: ", rank_min_suitability, "\n",
  "Dual threshold: ", dual_threshold, "\n",
  "Cores: ", rank_cores, "\n",
  "Reuse existing valid outputs: ", reuse_existing_outputs, "\n",
  sep = ""
)

# 3. Build complete job table ===================================================

job_list <- list()

for (method in methods_to_run) {
  for (scenario in scenarios_to_run) {
    files <- dual_file(method, scenario, zoneID)
    missing_zone <- zoneID[!file.exists(files)]
    out_valid <- valid_rank_outputs(
      method = method,
      scenario = scenario,
      template = r,
      k_rank = k_rank,
      n_summary_layers = n_summary_layers
    )
    
    job_list[[length(job_list) + 1L]] <- data.table(
      method = method,
      scenario = scenario,
      n_required_zones = length(zoneID),
      n_missing_zones = length(missing_zone),
      missing_zones = paste(missing_zone, collapse = ","),
      input_complete = length(missing_zone) == 0,
      output_valid = out_valid$all_valid,
      valid_ranked_zone = out_valid$valid_zone,
      valid_ranked_suitability = out_valid$valid_suit,
      valid_ranked_summary = out_valid$valid_summary,
      ranked_zone_file = out_valid$ranked_zone_file,
      ranked_suitability_file = out_valid$ranked_suitability_file,
      ranked_summary_file = out_valid$ranked_summary_file
    )
  }
}

jobs <- rbindlist(job_list)
setorder(jobs, method, scenario)

fwrite(
  jobs,
  file.path(table_dir, "ranking_jobs_all.csv")
)

fwrite(
  jobs[input_complete == TRUE],
  file.path(table_dir, "ranking_jobs_input_complete.csv")
)

fwrite(
  jobs[input_complete == FALSE & output_valid == FALSE],
  file.path(table_dir, "ranking_jobs_missing_input_and_output.csv")
)

cat(
  "\n[JOB STATUS]\n",
  "Total method-scenario jobs: ", nrow(jobs), "\n",
  "Already valid outputs: ", jobs[, sum(output_valid)], "\n",
  "Input-complete jobs: ", jobs[, sum(input_complete)], "\n",
  "Missing-input and no valid output jobs: ", jobs[, sum(!input_complete & !output_valid)], "\n",
  sep = ""
)

# 4. Rank dual suitability ======================================================

output_index <- list()
layer_index <- list()

for (j in seq_len(nrow(jobs))) {
  method <- jobs$method[j]
  scenario <- jobs$scenario[j]
  
  cat(
    "\n==============================\n",
    "RANK DUAL SUITABILITY: ",
    method,
    " | ",
    scenario,
    "\n",
    "==============================\n",
    sep = ""
  )
  
  out <- rank_output_files(method, scenario)
  dir.create(out$out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Reuse completed outputs first. This allows a full rerun to resume safely,
  # even if the original dual-suitability inputs have since been moved.
  if (reuse_existing_outputs && isTRUE(jobs$output_valid[j])) {
    cat0("[REUSE] valid outputs already exist: ", method, " | ", scenario)
    
    output_index[[length(output_index) + 1L]] <- data.table(
      method = method,
      scenario = scenario,
      status = "reused",
      message = NA_character_,
      n_input_zones = length(zoneID),
      n_missing_zones = jobs$n_missing_zones[j],
      missing_zones = jobs$missing_zones[j],
      top_k_rank = k_rank,
      rank_min_suitability = rank_min_suitability,
      dual_threshold = dual_threshold,
      ranked_zone_file = out$ranked_zone_file,
      ranked_suitability_file = out$ranked_suit_file,
      ranked_summary_file = out$ranked_summary_file
    )
    
    layer_index[[length(layer_index) + 1L]] <- make_layer_index(
      method,
      scenario,
      k_rank,
      summary_layer_names
    )
    
    next
  }
  
  if (!isTRUE(jobs$input_complete[j])) {
    msg <- paste0(
      "Missing input dual suitability rasters for zones: ",
      jobs$missing_zones[j]
    )
    
    cat0("[SKIP MISSING INPUT] ", method, " | ", scenario)
    cat0("  ", msg)
    
    output_index[[length(output_index) + 1L]] <- data.table(
      method = method,
      scenario = scenario,
      status = "skipped_missing_input",
      message = msg,
      n_input_zones = length(zoneID) - jobs$n_missing_zones[j],
      n_missing_zones = jobs$n_missing_zones[j],
      missing_zones = jobs$missing_zones[j],
      top_k_rank = k_rank,
      rank_min_suitability = rank_min_suitability,
      dual_threshold = dual_threshold,
      ranked_zone_file = out$ranked_zone_file,
      ranked_suitability_file = out$ranked_suit_file,
      ranked_summary_file = out$ranked_summary_file
    )
    
    next
  }
  
  tryCatch({
    files <- dual_file(method, scenario, zoneID)
    names(files) <- zoneID
    
    # Use a unique temporary filename each time.
    tmp_rank_all <- tempfile(
      pattern = paste0(
        "tmp_rank_all_",
        safe_name(method),
        "_",
        safe_name(scenario),
        "_top",
        k_rank,
        "_"
      ),
      tmpdir = tmp_dir,
      fileext = ".tif"
    )
    
    dual_stack <- read_dual_stack(
      files = files,
      template = r,
      valid_mask = valid_mask
    )
    
    if (overwrite_outputs) {
      remove_raster_files(out$ranked_zone_file, must_remove = TRUE)
      remove_raster_files(out$ranked_suit_file, must_remove = TRUE)
      remove_raster_files(out$ranked_summary_file, must_remove = TRUE)
    }
    
    remove_raster_files(tmp_rank_all, must_remove = TRUE)
    
    rank_fun <- make_rank_fun(
      zoneID = zoneID,
      rank_min = rank_min_suitability,
      threshold = dual_threshold,
      top_k = k_rank
    )
    
    # Main computation.
    # Temporary raster is uncompressed for speed.
    rank_all <- app(
      dual_stack,
      fun = rank_fun,
      filename = tmp_rank_all,
      overwrite = TRUE,
      cores = rank_cores,
      wopt = list(
        datatype = "FLT4S",
        gdal = "COMPRESS=NONE"
      )
    )
    
    zone_layers <- seq_len(k_rank)
    suit_layers <- k_rank + seq_len(k_rank)
    summary_layers <- 2 * k_rank + seq_len(n_summary_layers)
    
    ranked_zone <- rank_all[[zone_layers]]
    names(ranked_zone) <- rank_layer_names(k_rank, "zone")
    
    ranked_suit <- rank_all[[suit_layers]]
    names(ranked_suit) <- rank_layer_names(k_rank, "suit")
    
    ranked_summary <- rank_all[[summary_layers]]
    names(ranked_summary) <- summary_layer_names
    
    remove_raster_files(out$ranked_zone_file, must_remove = TRUE)
    writeRaster(
      ranked_zone,
      out$ranked_zone_file,
      overwrite = TRUE,
      wopt = list(
        datatype = "INT2S",
        gdal = "COMPRESS=LZW"
      )
    )
    
    remove_raster_files(out$ranked_suit_file, must_remove = TRUE)
    writeRaster(
      ranked_suit,
      out$ranked_suit_file,
      overwrite = TRUE,
      wopt = list(
        datatype = "FLT4S",
        gdal = "COMPRESS=LZW"
      )
    )
    
    remove_raster_files(out$ranked_summary_file, must_remove = TRUE)
    writeRaster(
      ranked_summary,
      out$ranked_summary_file,
      overwrite = TRUE,
      wopt = list(
        datatype = "FLT4S",
        gdal = "COMPRESS=LZW"
      )
    )
    
    output_index[[length(output_index) + 1L]] <- data.table(
      method = method,
      scenario = scenario,
      status = "created",
      message = NA_character_,
      n_input_zones = length(zoneID),
      n_missing_zones = 0L,
      missing_zones = "",
      top_k_rank = k_rank,
      rank_min_suitability = rank_min_suitability,
      dual_threshold = dual_threshold,
      ranked_zone_file = out$ranked_zone_file,
      ranked_suitability_file = out$ranked_suit_file,
      ranked_summary_file = out$ranked_summary_file
    )
    
    layer_index[[length(layer_index) + 1L]] <- make_layer_index(
      method,
      scenario,
      k_rank,
      summary_layer_names
    )
    
    cat(
      "[SAVED]\n",
      "  ", out$ranked_zone_file, "\n",
      "  ", out$ranked_suit_file, "\n",
      "  ", out$ranked_summary_file, "\n",
      sep = ""
    )
    
    rm(
      dual_stack,
      rank_fun,
      rank_all,
      ranked_zone,
      ranked_suit,
      ranked_summary
    )
    
    gc()
    remove_raster_files(tmp_rank_all)
  }, error = function(e) {
    msg <- conditionMessage(e)
    
    cat0("[ERROR] ", method, " | ", scenario)
    cat0("  ", msg)
    
    output_index[[length(output_index) + 1L]] <<- data.table(
      method = method,
      scenario = scenario,
      status = "error",
      message = msg,
      n_input_zones = length(zoneID),
      n_missing_zones = 0L,
      missing_zones = "",
      top_k_rank = k_rank,
      rank_min_suitability = rank_min_suitability,
      dual_threshold = dual_threshold,
      ranked_zone_file = out$ranked_zone_file,
      ranked_suitability_file = out$ranked_suit_file,
      ranked_summary_file = out$ranked_summary_file
    )
    
    gc()
  })
}

# 5. Save index tables ==========================================================

output_index <- rbindlist(output_index, fill = TRUE)

if (length(layer_index) > 0) {
  layer_index <- rbindlist(layer_index, fill = TRUE)
} else {
  layer_index <- data.table(
    method = character(),
    scenario = character(),
    raster_type = character(),
    layer = integer(),
    layer_name = character(),
    rank = integer(),
    meaning = character()
  )
}

setorder(output_index, method, scenario)
setorder(layer_index, method, scenario, raster_type, layer)

fwrite(
  output_index,
  file.path(table_dir, "ranking_output_index.csv")
)

fwrite(
  layer_index,
  file.path(table_dir, "ranking_layer_index.csv")
)

cat(
  "\nCOMPLETE\n",
  "Output root: ", output_root, "\n",
  "Jobs total: ", nrow(output_index), "\n",
  "Created: ", output_index[, sum(status == "created")], "\n",
  "Reused: ", output_index[, sum(status == "reused")], "\n",
  "Skipped missing input: ", output_index[, sum(status == "skipped_missing_input")], "\n",
  "Errors: ", output_index[, sum(status == "error")], "\n",
  "Top K ranks saved: ", k_rank, "\n",
  "Rank minimum suitability: ", rank_min_suitability, "\n",
  "Dual threshold: ", dual_threshold, "\n",
  "Summary tables saved to: ", table_dir, "\n",
  sep = ""
)
