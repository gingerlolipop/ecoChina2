# 5.2 Quick diagnostic: raster alignment + normal-map agreement
# ==============================================================================
# This script does NOT refit or reproject models.
#
# It checks:
#   1. Raster geometry/alignment for vegetation, raw climate, raw soil, and
#      normal-period suitability/dual/assigned outputs.
#   2. Whether climate training rows still match the raw climate rasters at the
#      stored x/y coordinates.
#   3. Whether soil training rows match BOTH the vegetation zone and raw soil
#      rasters at their stored x/y coordinates.
#   4. On a large random map sample, whether low exact agreement is explained by
#      missing original-zone dual suitability, incomplete 53-zone competition,
#      or genuine rank competition among ecotypes.
#
# Important:
#   Different raster origins/resolutions are NOT automatically an error when
#   values are extracted by geographic coordinates. The training-value checks
#   below test the important question directly: are the values attached to each
#   training row actually the values at that row's geographic location?
# ==============================================================================

library(terra)
library(data.table)

rm(list = ls())
gc()

set.seed(49)


# 0. Paths and settings =========================================================

base_dir <- "H:/Jing/ecoChina2"

reference_file <- file.path(
  base_dir,
  "raster/ecosys_ori.tif"
)

climate_raw_dir <- paste0(
  "H:/Jing/ecoChina/play/China/ClimateData/CN/800m/",
  "Normal_1961_1990"
)

soil_raw_dir <- "H:/Jing/soil rasters/tif2"

climate_train_file <- file.path(
  base_dir,
  "results/train_data.csv"
)

soil_train_file <- file.path(
  base_dir,
  "results/soil_train_data.csv"
)

soil_coord_file <- file.path(
  base_dir,
  "results/soil_train_coords.csv"
)

out_dir <- file.path(
  base_dir,
  "assessment_var",
  "normal_map_diagnostic"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

zoneID <- c(
  1:7,
  9:50,
  52:55
)

method_order <- c(
  "rf_var",
  "mf_var"
)

method_labels <- c(
  rf_var = "Plain RF",
  mf_var = "Plain MF RF"
)

soil_source <- c(
  rf_var = "plain_rf",
  mf_var = "plain_mf"
)

tie_tol <- 1e-4

# Training-coordinate checks are sample based and should run quickly.
TRAIN_SAMPLE_N <- 10000L

# Random map sample for the 53-zone ranking diagnostic.
MAP_SAMPLE_N <- 100000L


# 1. Helpers ====================================================================

safe_div <- function(a, b) {
  ifelse(
    is.finite(b) & b > 0,
    a / b,
    NA_real_
  )
}


extract_lonlat <- function(x, xy) {
  
  xy <- as.data.frame(xy)
  
  if (!all(c("x", "y") %in% names(xy))) {
    stop("Coordinate table must contain x and y.")
  }
  
  pts <- vect(
    xy,
    geom = c("x", "y"),
    crs = "EPSG:4326"
  )
  
  x_crs <- crs(x)
  
  if (!is.na(x_crs) && nzchar(x_crs)) {
    pts <- project(
      pts,
      x_crs
    )
  }
  
  as.data.table(
    extract(
      x,
      pts,
      ID = FALSE
    )
  )
}


geometry_row <- function(
    file,
    group,
    reference) {
  
  if (!file.exists(file)) {
    return(
      data.table(
        group = group,
        file = file,
        exists = FALSE
      )
    )
  }
  
  x <- rast(file)
  
  same_geometry <- tryCatch(
    compareGeom(
      reference,
      x,
      stopOnError = FALSE
    ),
    error = function(e) FALSE
  )
  
  ref_res <- res(reference)
  x_res <- res(x)
  
  ref_origin <- origin(reference)
  x_origin <- origin(x)
  
  data.table(
    group = group,
    file = file,
    exists = TRUE,
    same_geometry_as_vegetation = same_geometry,
    
    nrow = nrow(x),
    ncol = ncol(x),
    
    res_x = x_res[1],
    res_y = x_res[2],
    
    ref_res_x = ref_res[1],
    ref_res_y = ref_res[2],
    
    origin_x = x_origin[1],
    origin_y = x_origin[2],
    
    ref_origin_x = ref_origin[1],
    ref_origin_y = ref_origin[2],
    
    origin_offset_ref_cells_x =
      (x_origin[1] - ref_origin[1]) /
      ref_res[1],
    
    origin_offset_ref_cells_y =
      (x_origin[2] - ref_origin[2]) /
      ref_res[2],
    
    xmin = xmin(x),
    xmax = xmax(x),
    ymin = ymin(x),
    ymax = ymax(x),
    
    ref_xmin = xmin(reference),
    ref_xmax = xmax(reference),
    ref_ymin = ymin(reference),
    ref_ymax = ymax(reference)
  )
}


compare_training_predictors <- function(
    dat,
    xy,
    raster_dir,
    dataset_name,
    sample_n) {
  
  if (!dir.exists(raster_dir)) {
    return(
      data.table(
        dataset = dataset_name,
        issue = paste(
          "Raster directory not found:",
          raster_dir
        )
      )
    )
  }
  
  raster_files <- list.files(
    raster_dir,
    pattern = "\\.tif$",
    full.names = TRUE
  )
  
  raster_names <- tools::file_path_sans_ext(
    basename(raster_files)
  )
  
  common_vars <- intersect(
    names(dat),
    raster_names
  )
  
  if (!length(common_vars)) {
    return(
      data.table(
        dataset = dataset_name,
        issue = "No training columns matched raster filenames."
      )
    )
  }
  
  n <- nrow(dat)
  
  sample_index <- sample.int(
    n,
    min(sample_n, n)
  )
  
  xy_sample <- xy[
    sample_index,
    c("x", "y"),
    drop = FALSE
  ]
  
  matched_files <- raster_files[
    match(
      common_vars,
      raster_names
    )
  ]
  
  r <- rast(matched_files)
  names(r) <- common_vars
  
  raster_values <- extract_lonlat(
    r,
    xy_sample
  )
  
  results <- lapply(
    common_vars,
    function(v) {
      
      a <- as.numeric(
        dat[[v]][sample_index]
      )
      
      b <- as.numeric(
        raster_values[[v]]
      )
      
      finite <- (
        is.finite(a) &
          is.finite(b)
      )
      
      n_finite <- sum(finite)
      
      if (n_finite > 0) {
        
        difference <- abs(
          a[finite] -
            b[finite]
        )
        
        tolerance <- 1e-6 *
          pmax(
            1,
            abs(a[finite]),
            abs(b[finite])
          )
        
        match_tolerance <- mean(
          difference <= tolerance
        )
        
        correlation <- if (
          n_finite > 2 &&
          sd(a[finite]) > 0 &&
          sd(b[finite]) > 0
        ) {
          cor(
            a[finite],
            b[finite]
          )
        } else {
          NA_real_
        }
        
        mean_abs_difference <- mean(
          difference
        )
        
        max_abs_difference <- max(
          difference
        )
        
      } else {
        
        match_tolerance <- NA_real_
        correlation <- NA_real_
        mean_abs_difference <- NA_real_
        max_abs_difference <- NA_real_
      }
      
      data.table(
        dataset = dataset_name,
        variable = v,
        n_sample = length(a),
        n_training_finite = sum(
          is.finite(a)
        ),
        n_raster_finite = sum(
          is.finite(b)
        ),
        n_both_finite = n_finite,
        finite_mask_match = mean(
          is.finite(a) ==
            is.finite(b)
        ),
        value_match_with_tolerance =
          match_tolerance,
        correlation = correlation,
        mean_absolute_difference =
          mean_abs_difference,
        max_absolute_difference =
          max_abs_difference
      )
    }
  )
  
  rbindlist(
    results,
    fill = TRUE
  )
}


zone_coordinate_check <- function(
    zone,
    xy,
    reference,
    dataset_name,
    sample_n) {
  
  n <- length(zone)
  
  sample_index <- sample.int(
    n,
    min(sample_n, n)
  )
  
  ref_value <- extract_lonlat(
    reference,
    xy[
      sample_index,
      c("x", "y"),
      drop = FALSE
    ]
  )[[1]]
  
  observed_zone <- as.integer(
    as.character(
      zone[sample_index]
    )
  )
  
  valid <- (
    is.finite(observed_zone) &
      is.finite(ref_value)
  )
  
  data.table(
    dataset = dataset_name,
    n_sample = length(sample_index),
    n_both_valid = sum(valid),
    zone_match_share = if (
      any(valid)
    ) {
      mean(
        observed_zone[valid] ==
          as.integer(ref_value[valid])
      )
    } else {
      NA_real_
    },
    n_zone_mismatch = if (
      any(valid)
    ) {
      sum(
        observed_zone[valid] !=
          as.integer(ref_value[valid])
      )
    } else {
      NA_integer_
    }
  )
}


assigned_map_file <- function(method) {
  
  file.path(
    base_dir,
    "result maps",
    method,
    paste0(
      "assigned_zone_normal_threshold0.4_tol",
      tie_tol,
      "_novel99_maskNA8_noNovelNormal.tif"
    )
  )
}


dual_files <- function(method) {
  
  file.path(
    base_dir,
    "dual suit",
    method,
    "normal",
    paste0(
      "dual_suitability_zone",
      zoneID,
      ".tif"
    )
  )
}


# 2. Geometry and raster alignment ==============================================

reference_map <- rast(
  reference_file
)

names(reference_map) <- "original_zone"

geometry_results <- list()

# Raw climate rasters
climate_files <- list.files(
  climate_raw_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

for (f in climate_files) {
  geometry_results[[
    length(geometry_results) + 1L
  ]] <- geometry_row(
    f,
    "raw_climate",
    reference_map
  )
}

# Raw soil rasters
soil_files <- list.files(
  soil_raw_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

for (f in soil_files) {
  geometry_results[[
    length(geometry_results) + 1L
  ]] <- geometry_row(
    f,
    "raw_soil",
    reference_map
  )
}

# Final normal-period products
for (method in method_order) {
  
  geometry_results[[
    length(geometry_results) + 1L
  ]] <- geometry_row(
    assigned_map_file(method),
    paste0(
      method,
      "_assigned_normal"
    ),
    reference_map
  )
  
  for (z in zoneID) {
    
    geometry_results[[
      length(geometry_results) + 1L
    ]] <- geometry_row(
      file.path(
        base_dir,
        "clim suitability",
        method,
        "normal",
        paste0(
          "clim_suit_zone",
          z,
          ".tif"
        )
      ),
      paste0(
        method,
        "_climate_suitability"
      ),
      reference_map
    )
    
    geometry_results[[
      length(geometry_results) + 1L
    ]] <- geometry_row(
      file.path(
        base_dir,
        "soil suitability",
        soil_source[[method]],
        "normal",
        paste0(
          "soil_suit_zone",
          z,
          ".tif"
        )
      ),
      paste0(
        method,
        "_soil_suitability"
      ),
      reference_map
    )
    
    geometry_results[[
      length(geometry_results) + 1L
    ]] <- geometry_row(
      file.path(
        base_dir,
        "dual suit",
        method,
        "normal",
        paste0(
          "dual_suitability_zone",
          z,
          ".tif"
        )
      ),
      paste0(
        method,
        "_dual_suitability"
      ),
      reference_map
    )
  }
}

geometry_table <- rbindlist(
  geometry_results,
  fill = TRUE
)

fwrite(
  geometry_table,
  file.path(
    out_dir,
    "diagnostic_raster_geometry.csv"
  )
)

geometry_summary <- geometry_table[
  ,
  .(
    n_files = .N,
    n_missing = sum(
      !exists
    ),
    n_same_geometry = sum(
      exists &
        same_geometry_as_vegetation,
      na.rm = TRUE
    ),
    n_different_geometry = sum(
      exists &
        !same_geometry_as_vegetation,
      na.rm = TRUE
    )
  ),
  by = group
]

fwrite(
  geometry_summary,
  file.path(
    out_dir,
    "diagnostic_raster_geometry_summary.csv"
  )
)

cat(
  "\n[RASTER GEOMETRY SUMMARY]\n"
)

print(
  geometry_summary
)


# 3. Training-data coordinate checks ============================================

zone_checks <- list()
predictor_checks <- list()


# 3.1 Climate training -----------------------------------------------------------

if (
  file.exists(climate_train_file)
) {
  
  climate_train <- fread(
    climate_train_file
  )
  
  if (
    all(
      c(
        "zoneID",
        "x",
        "y"
      ) %in%
      names(climate_train)
    )
  ) {
    
    climate_xy <- as.data.frame(
      climate_train[
        ,
        .(
          x,
          y
        )
      ]
    )
    
    zone_checks[[
      length(zone_checks) + 1L
    ]] <- zone_coordinate_check(
      zone = climate_train$zoneID,
      xy = climate_xy,
      reference = reference_map,
      dataset_name = "climate_training",
      sample_n = TRAIN_SAMPLE_N
    )
    
    predictor_checks[[
      length(predictor_checks) + 1L
    ]] <- compare_training_predictors(
      dat = climate_train,
      xy = climate_xy,
      raster_dir = climate_raw_dir,
      dataset_name = "climate_training",
      sample_n = TRAIN_SAMPLE_N
    )
  }
}


# 3.2 Soil training --------------------------------------------------------------

if (
  file.exists(soil_train_file) &&
  file.exists(soil_coord_file)
) {
  
  soil_train <- fread(
    soil_train_file
  )
  
  soil_xy <- fread(
    soil_coord_file
  )
  
  if (
    nrow(soil_train) !=
    nrow(soil_xy)
  ) {
    stop(
      "soil_train_data.csv and soil_train_coords.csv ",
      "have different row counts."
    )
  }
  
  if (
    !all(
      c(
        "x",
        "y"
      ) %in%
      names(soil_xy)
    )
  ) {
    stop(
      "soil_train_coords.csv must contain x and y."
    )
  }
  
  soil_xy_df <- as.data.frame(
    soil_xy[
      ,
      .(
        x,
        y
      )
    ]
  )
  
  if (
    "zoneID" %in%
    names(soil_train)
  ) {
    
    zone_checks[[
      length(zone_checks) + 1L
    ]] <- zone_coordinate_check(
      zone = soil_train$zoneID,
      xy = soil_xy_df,
      reference = reference_map,
      dataset_name = "soil_training",
      sample_n = TRAIN_SAMPLE_N
    )
  }
  
  predictor_checks[[
    length(predictor_checks) + 1L
  ]] <- compare_training_predictors(
    dat = soil_train,
    xy = soil_xy_df,
    raster_dir = soil_raw_dir,
    dataset_name = "soil_training",
    sample_n = TRAIN_SAMPLE_N
  )
}


zone_check_table <- rbindlist(
  zone_checks,
  fill = TRUE
)

predictor_check_table <- rbindlist(
  predictor_checks,
  fill = TRUE
)

fwrite(
  zone_check_table,
  file.path(
    out_dir,
    "diagnostic_training_zone_at_coordinates.csv"
  )
)

fwrite(
  predictor_check_table,
  file.path(
    out_dir,
    "diagnostic_training_predictors_at_coordinates.csv"
  )
)

cat(
  "\n[TRAINING ZONE / COORDINATE CHECK]\n"
)

print(
  zone_check_table
)

cat(
  "\n[TRAINING PREDICTOR CHECK - WORST MATCHES]\n"
)

if (
  nrow(predictor_check_table) > 0 &&
  "value_match_with_tolerance" %in%
  names(predictor_check_table)
) {
  
  print(
    predictor_check_table[
      order(
        value_match_with_tolerance,
        correlation
      )
    ][
      1:min(
        .N,
        20L
      )
    ]
  )
}


# 4. Fast random-sample map diagnostic ==========================================

# Oversample, then retain only the 53 modeled zones.
sample_reference <- spatSample(
  reference_map,
  size = as.integer(
    MAP_SAMPLE_N * 1.25
  ),
  method = "random",
  na.rm = TRUE,
  xy = TRUE,
  cells = TRUE,
  as.df = TRUE
)

sample_reference <- as.data.table(
  sample_reference
)

value_column <- setdiff(
  names(sample_reference),
  c(
    "cell",
    "x",
    "y"
  )
)[1]

setnames(
  sample_reference,
  value_column,
  "original_zone"
)

sample_reference[
  ,
  original_zone :=
    as.integer(
      original_zone
    )
]

sample_reference <- sample_reference[
  original_zone %in%
    zoneID
]

if (
  nrow(sample_reference) >
  MAP_SAMPLE_N
) {
  sample_reference <- sample_reference[
    sample.int(
      .N,
      MAP_SAMPLE_N
    )
  ]
}

sample_xy <- as.data.frame(
  sample_reference[
    ,
    .(
      x,
      y
    )
  ]
)

sample_overall_results <- list()
sample_zone_results <- list()


for (method in method_order) {
  
  cat(
    "\n========================================\n",
    "MAP SAMPLE: ",
    method_labels[[method]],
    "\n========================================\n",
    sep = ""
  )
  
  prediction_file <- assigned_map_file(
    method
  )
  
  files_dual <- dual_files(
    method
  )
  
  if (
    !file.exists(prediction_file) ||
    !all(
      file.exists(
        files_dual
      )
    )
  ) {
    stop(
      "Missing normal-period output for ",
      method
    )
  }
  
  predicted_map <- rast(
    prediction_file
  )
  
  dual_stack <- rast(
    files_dual
  )
  
  names(
    dual_stack
  ) <- as.character(
    zoneID
  )
  
  predicted <- extract_lonlat(
    predicted_map,
    sample_xy
  )[[1]]
  
  suitability <- as.matrix(
    extract_lonlat(
      dual_stack,
      sample_xy
    )
  )
  
  original <- sample_reference$original_zone
  
  original_column <- match(
    original,
    zoneID
  )
  
  n_sample <- length(
    original
  )
  
  predicted_valid <- (
    !is.na(predicted) &
      predicted %in%
      zoneID
  )
  
  finite_suitability <- is.finite(
    suitability
  )
  
  n_valid_dual <- rowSums(
    finite_suitability
  )
  
  reference_dual <- suitability[
    cbind(
      seq_len(
        n_sample
      ),
      original_column
    )
  ]
  
  reference_dual_valid <- is.finite(
    reference_dual
  )
  
  exact <- (
    predicted_valid &
      predicted ==
      original
  )
  
  # Rank of the original zone. Rank 1 includes ties within tie_tol.
  suitability_rank <- suitability
  
  suitability_rank[
    !is.finite(
      suitability_rank
    )
  ] <- -Inf
  
  original_rank <- rep(
    NA_integer_,
    n_sample
  )
  
  rank_rows <- which(
    reference_dual_valid
  )
  
  if (
    length(
      rank_rows
    ) > 0
  ) {
    
    higher <- sweep(
      suitability_rank[
        rank_rows,
        ,
        drop = FALSE
      ],
      1L,
      reference_dual[
        rank_rows
      ] + tie_tol,
      FUN = ">"
    )
    
    original_rank[
      rank_rows
    ] <- 1L +
      rowSums(
        higher
      )
  }
  
  current_mask <- predicted_valid
  
  reference_dual_mask <- (
    predicted_valid &
      reference_dual_valid
  )
  
  all_53_mask <- (
    predicted_valid &
      n_valid_dual ==
      length(
        zoneID
      )
  )
  
  # Under script 4.1, rank 1 means the original zone should be retained.
  consistency_violation <- (
    reference_dual_mask &
      (
        (
          original_rank == 1L
        ) !=
          exact
      )
  )
  
  sample_overall <- data.table(
    method = method,
    method_label =
      method_labels[[method]],
    
    sample_pixels = n_sample,
    
    current_prediction_valid_share =
      mean(
        predicted_valid
      ),
    
    current_exact_share =
      safe_div(
        sum(
          current_mask &
            exact
        ),
        sum(
          current_mask
        )
      ),
    
    reference_zone_dual_valid_share =
      mean(
        reference_dual_valid
      ),
    
    predicted_valid_but_reference_dual_NA_share =
      mean(
        predicted_valid &
          !reference_dual_valid
      ),
    
    reference_zone_dual_exact_share =
      safe_div(
        sum(
          reference_dual_mask &
            exact
        ),
        sum(
          reference_dual_mask
        )
      ),
    
    all_53_dual_valid_share =
      mean(
        n_valid_dual ==
          length(
            zoneID
          )
      ),
    
    all_53_exact_share =
      safe_div(
        sum(
          all_53_mask &
            exact
        ),
        sum(
          all_53_mask
        )
      ),
    
    original_zone_top1_share =
      safe_div(
        sum(
          reference_dual_mask &
            original_rank <= 1L,
          na.rm = TRUE
        ),
        sum(
          reference_dual_mask
        )
      ),
    
    original_zone_top2_share =
      safe_div(
        sum(
          reference_dual_mask &
            original_rank <= 2L,
          na.rm = TRUE
        ),
        sum(
          reference_dual_mask
        )
      ),
    
    original_zone_top3_share =
      safe_div(
        sum(
          reference_dual_mask &
            original_rank <= 3L,
          na.rm = TRUE
        ),
        sum(
          reference_dual_mask
        )
      ),
    
    original_zone_top5_share =
      safe_div(
        sum(
          reference_dual_mask &
            original_rank <= 5L,
          na.rm = TRUE
        ),
        sum(
          reference_dual_mask
        )
      ),
    
    assignment_consistency_violations =
      sum(
        consistency_violation,
        na.rm = TRUE
      )
  )
  
  sample_overall_results[[
    method
  ]] <- sample_overall
  
  
  sample_zone <- data.table(
    original_zone = original,
    predicted_valid = predicted_valid,
    exact = exact,
    reference_dual_valid =
      reference_dual_valid,
    n_valid_dual =
      n_valid_dual,
    original_rank =
      original_rank
  )[
    ,
    .(
      n_sample = .N,
      
      prediction_valid_share =
        mean(
          predicted_valid
        ),
      
      exact_share =
        safe_div(
          sum(
            predicted_valid &
              exact
          ),
          sum(
            predicted_valid
          )
        ),
      
      reference_dual_valid_share =
        mean(
          reference_dual_valid
        ),
      
      predicted_valid_but_reference_dual_NA_share =
        mean(
          predicted_valid &
            !reference_dual_valid
        ),
      
      exact_when_reference_dual_valid =
        safe_div(
          sum(
            predicted_valid &
              reference_dual_valid &
              exact
          ),
          sum(
            predicted_valid &
              reference_dual_valid
          )
        ),
      
      original_zone_top3_share =
        safe_div(
          sum(
            predicted_valid &
              reference_dual_valid &
              original_rank <= 3L,
            na.rm = TRUE
          ),
          sum(
            predicted_valid &
              reference_dual_valid
          )
        ),
      
      original_zone_top5_share =
        safe_div(
          sum(
            predicted_valid &
              reference_dual_valid &
              original_rank <= 5L,
            na.rm = TRUE
          ),
          sum(
            predicted_valid &
              reference_dual_valid
          )
        )
    ),
    by = original_zone
  ]
  
  sample_zone[
    ,
    `:=`(
      method = method,
      method_label =
        method_labels[[method]]
    )
  ]
  
  setcolorder(
    sample_zone,
    c(
      "method",
      "method_label",
      "original_zone",
      setdiff(
        names(sample_zone),
        c(
          "method",
          "method_label",
          "original_zone"
        )
      )
    )
  )
  
  sample_zone_results[[
    method
  ]] <- sample_zone
  
  print(
    sample_overall
  )
  
  rm(
    predicted_map,
    dual_stack,
    suitability,
    suitability_rank
  )
  
  gc()
}


sample_overall_table <- rbindlist(
  sample_overall_results,
  fill = TRUE
)

sample_zone_table <- rbindlist(
  sample_zone_results,
  fill = TRUE
)

setorder(
  sample_zone_table,
  method,
  original_zone
)

fwrite(
  sample_overall_table,
  file.path(
    out_dir,
    "diagnostic_sample_map_overall.csv"
  )
)

fwrite(
  sample_zone_table,
  file.path(
    out_dir,
    "diagnostic_sample_map_by_original_zone.csv"
  )
)


# 5. Console interpretation guide ===============================================

cat(
  "\n============================================================\n",
  "WHAT TO CHECK FIRST\n",
  "============================================================\n",
  "1. diagnostic_training_zone_at_coordinates.csv\n",
  "   zone_match_share should be ~1 for both climate and soil training.\n\n",
  "2. diagnostic_training_predictors_at_coordinates.csv\n",
  "   value_match_with_tolerance and correlation should be ~1.\n",
  "   If soil is poor here, the soil training table is spatially mismatched\n",
  "   to the raw soil rasters even if RF test metrics are high.\n\n",
  "3. diagnostic_raster_geometry_summary.csv\n",
  "   Final climate/soil/dual/assigned products should match vegetation.\n",
  "   Raw climate/soil rasters may differ in origin/resolution; that alone is\n",
  "   not an error because training values can be extracted by coordinates.\n\n",
  "4. diagnostic_sample_map_overall.csv\n",
  "   assignment_consistency_violations should be 0.\n",
  "   Compare current_exact_share with reference_zone_dual_exact_share and\n",
  "   all_53_exact_share, then inspect top-3/top-5 shares.\n",
  "============================================================\n",
  sep = ""
)

cat(
  "\nCOMPLETE\n",
  "Output folder: ",
  out_dir,
  "\n",
  sep = ""
)
