# 5.3 Feature-importance analysis for the final ecotype RF models
# Compatibility revision v4: fixes terra::freq() compatibility, data.table SE
# assignment, and grouped result types when category medoid values are absent.
# ==============================================================================
# Purpose
# -------
# Extract and summarize the stored feature importance from the final climate
# and soil models without refitting any model or recalculating projections.
#
# Models included
# ---------------
#   climate: Plain RF (rf_var) and Plain MF RF (mf_var)
#   soil:    Plain RF (rf_var) and Plain MF RF (mf_var)
#
# Primary importance metric
# -------------------------
# The primary metric is the unscaled global MeanDecreaseAccuracy (MDA) stored
# by randomForest. Positive MDA values are normalized to sum to one within each
# zone x niche x method model. Raw negative values are retained in the long
# table, but contribute zero to normalized shares. MeanDecreaseGini is exported
# only as a secondary diagnostic because it is more sensitive to predictor
# distributions and available split points.
#
# Interpretation
# --------------
# Importance identifies variables that contribute to prediction. It does not
# provide response direction, response shape, mechanism, or causality. Climate
# predictors are correlated, so category-level inference should emphasize
# recurring variables, ecological factor groups, and agreement across RF and
# MF RF rather than a single rank in a single model.
#
# Main outputs
# ------------
#   - complete model and importance audit
#   - zone-level importance and Top-10 factors
#   - RF versus MF RF rank agreement
#   - common variables and factor groups within broad vegetation categories
#   - contrasts among categories
#   - category-profile medoids and observed-area extremes
#   - category-specific figures for representative, largest, and smallest zones
#
# Run after scripts 2.16 and 2.4. The script reads saved .Rdata models only.
# ==============================================================================

library(data.table)
library(terra)
library(randomForest)
library(ggplot2)

rm(list = ls())
gc()


# 0. Paths and settings =========================================================

base_dir <- "H:/Jing/ecoChina2"

climate_model_dir <- file.path(
  base_dir,
  "rf"
)

soil_model_dir <- file.path(
  base_dir,
  "rf_soil"
)

palette_file <- file.path(
  base_dir,
  "color_palette_China.csv"
)

reference_raster_file <- file.path(
  base_dir,
  "raster",
  "ecosys_ori.tif"
)

output_root <- file.path(
  base_dir,
  "assessment_var",
  "feature_importance_analysis"
)

table_dir <- file.path(
  output_root,
  "tables"
)

figure_dir <- file.path(
  output_root,
  "figures"
)

representative_figure_dir <- file.path(
  figure_dir,
  "category representative zones"
)

for (directory in c(
  output_root,
  table_dir,
  figure_dir,
  representative_figure_dir
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

model_zoneID <- c(
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

niche_order <- c(
  "climate",
  "soil"
)

top_n_zone_factors <- 10L
top_n_category_factors <- 6L
top_n_theory_groups <- 3L
n_global_area_extremes <- 5L

# A category-level variable is flagged as common only when it is selected in at
# least half of the category's zones and enters the consensus Top-5 in at least
# one quarter. The full continuous frequencies are always exported so these
# descriptive cutoffs can be changed without rerunning the models.
common_selection_frequency <- 0.50
common_top5_frequency <- 0.25

# Categories represented by fewer zones remain in all outputs, but their
# summaries are marked descriptive_only rather than within-category evidence.
minimum_zones_for_category_inference <- 3L

# Stop before aggregation if any of the 53 x 2 x 2 expected model files fails.
strict_complete_models <- TRUE


# 1. Helpers ====================================================================

require_file <- function(file) {
  if (!file.exists(file)) {
    stop(
      "Missing required file: ",
      file
    )
  }
  
  file
}


safe_name <- function(x) {
  x <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    x
  )
  
  x <- gsub(
    "^_+|_+$",
    "",
    x
  )
  
  ifelse(
    nzchar(x),
    x,
    "unnamed"
  )
}


wrap_text <- function(
    x,
    width = 34L) {
  vapply(
    x,
    function(value) {
      paste(
        strwrap(
          value,
          width = width
        ),
        collapse = "\n"
      )
    },
    character(1)
  )
}


percent_label <- function(x) {
  paste0(
    round(
      100 * x,
      1
    ),
    "%"
  )
}


save_plot <- function(
    plot_object,
    filename,
    width,
    height) {
  ggsave(
    filename = filename,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}


safe_cor <- function(
    x,
    y,
    method = "spearman") {
  keep <- is.finite(x) &
    is.finite(y)
  
  if (sum(keep) < 3L ||
      length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) {
    return(NA_real_)
  }
  
  suppressWarnings(
    cor(
      x[keep],
      y[keep],
      method = method
    )
  )
}


cosine_similarity <- function(
    x,
    y) {
  denominator <- sqrt(
    sum(x ^ 2)
  ) * sqrt(
    sum(y ^ 2)
  )
  
  if (!is.finite(denominator) ||
      denominator <= 0) {
    return(NA_real_)
  }
  
  sum(x * y) /
    denominator
}


get_importance_column <- function(
    importance_matrix,
    column_name) {
  if (column_name %in%
      colnames(importance_matrix)) {
    return(
      as.numeric(
        importance_matrix[
          ,
          column_name
        ]
      )
    )
  }
  
  rep(
    NA_real_,
    nrow(importance_matrix)
  )
}


annotate_features <- function(
    niche,
    feature) {
  niche <- as.character(niche)
  feature <- as.character(feature)
  feature_standard <- toupper(
    gsub(
      "-",
      "_",
      feature,
      fixed = TRUE
    )
  )
  
  feature_group <- rep(
    "Other",
    length(feature)
  )
  
  climate <- niche == "climate"
  soil <- niche == "soil"
  
  feature_group[
    climate &
      grepl(
        "^(MAT|MWMT|MCMT|TD$|EMT|EXT|TMAX|TMIN|TAVE)",
        feature_standard
      )
  ] <- "Temperature"
  
  feature_group[
    climate &
      grepl(
        "^(MAP|MSP|PPT)",
        feature_standard
      )
  ] <- "Precipitation"
  
  feature_group[
    climate &
      grepl(
        "^(AHM|SHM|CMD|CMI)",
        feature_standard
      )
  ] <- "Climatic moisture balance"
  
  feature_group[
    climate &
      grepl(
        "^(BFFP|EFFP|FFP|NFFD|DD)",
        feature_standard
      )
  ] <- "Growing season and degree-days"
  
  feature_group[
    climate &
      grepl(
        "^(EREF|RSDS)",
        feature_standard
      )
  ] <- "Radiation and evaporative demand"
  
  feature_group[
    climate &
      grepl(
        "^RH",
        feature_standard
      )
  ] <- "Humidity"
  
  feature_group[
    climate &
      grepl(
        "^PAS",
        feature_standard
      )
  ] <- "Snowfall"
  
  feature_group[
    soil &
      grepl(
        "CEC|TEB|BASE|(^|_)BS($|_)",
        feature_standard
      )
  ] <- "Nutrient retention and base status"
  
  feature_group[
    soil &
      grepl(
        "BULK|DENS",
        feature_standard
      )
  ] <- "Bulk density"
  
  feature_group[
    soil &
      grepl(
        "GRAVEL|COARSE",
        feature_standard
      )
  ] <- "Coarse fragments"
  
  feature_group[
    soil &
      feature_group == "Other" &
      grepl(
        "SAND|SILT|CLAY|TEXT",
        feature_standard
      )
  ] <- "Texture"
  
  feature_group[
    soil &
      grepl(
        "(^|_)OC($|_)|ORGANIC",
        feature_standard
      )
  ] <- "Organic carbon"
  
  feature_group[
    soil &
      grepl(
        "PH",
        feature_standard
      )
  ] <- "Soil reaction"
  
  feature_group[
    soil &
      grepl(
        "CACO3|CARBONATE",
        feature_standard
      )
  ] <- "Carbonates"
  
  feature_group[
    soil &
      grepl(
        "CASO4|GYPS",
        feature_standard
      )
  ] <- "Gypsum"
  
  feature_group[
    soil &
      grepl(
        "ESP|SOD",
        feature_standard
      )
  ] <- "Sodicity"
  
  feature_group[
    soil &
      grepl(
        "(^|_)ECE($|_)|SALIN|ELECTRICAL",
        feature_standard
      )
  ] <- "Salinity"
  
  seasonal_suffix <- tolower(
    sub(
      "^.*_(wt|sp|sm|at)$",
      "\\1",
      feature
    )
  )
  
  has_seasonal_suffix <- grepl(
    "_(wt|sp|sm|at)$",
    tolower(feature)
  )
  
  climate_period <- rep(
    NA_character_,
    length(feature)
  )
  
  climate_period[climate] <- "Annual"
  
  climate_period[
    climate &
      has_seasonal_suffix &
      seasonal_suffix == "wt"
  ] <- "Winter"
  
  climate_period[
    climate &
      has_seasonal_suffix &
      seasonal_suffix == "sp"
  ] <- "Spring"
  
  climate_period[
    climate &
      has_seasonal_suffix &
      seasonal_suffix == "sm"
  ] <- "Summer"
  
  climate_period[
    climate &
      has_seasonal_suffix &
      seasonal_suffix == "at"
  ] <- "Autumn"
  
  data.table(
    niche = niche,
    feature = feature,
    feature_group = feature_group,
    climate_period = climate_period
  )
}


extract_model_importance <- function(
    model_file,
    expected_object,
    zoneID,
    niche,
    method) {
  model_environment <- new.env(
    parent = baseenv()
  )
  
  loaded_objects <- load(
    model_file,
    envir = model_environment
  )
  
  object_name <- expected_object
  
  if (!(object_name %in%
        loaded_objects)) {
    random_forest_objects <- loaded_objects[
      vapply(
        loaded_objects,
        function(name) {
          inherits(
            get(
              name,
              envir = model_environment
            ),
            "randomForest"
          )
        },
        logical(1)
      )
    ]
    
    if (length(random_forest_objects) != 1L) {
      stop(
        "Could not identify one randomForest object in ",
        model_file
      )
    }
    
    object_name <- random_forest_objects[[1]]
  }
  
  fitted_model <- get(
    object_name,
    envir = model_environment
  )
  
  if (!inherits(
    fitted_model,
    "randomForest"
  )) {
    stop(
      "Object is not a randomForest model: ",
      model_file,
      " | ",
      object_name
    )
  }
  
  importance_unscaled <- randomForest::importance(
    fitted_model,
    scale = FALSE
  )
  
  importance_scaled <- randomForest::importance(
    fitted_model,
    scale = TRUE
  )
  
  if (is.null(dim(importance_unscaled))) {
    importance_unscaled <- matrix(
      importance_unscaled,
      ncol = 1L,
      dimnames = list(
        names(importance_unscaled),
        "MeanDecreaseAccuracy"
      )
    )
  }
  
  if (is.null(dim(importance_scaled))) {
    importance_scaled <- matrix(
      importance_scaled,
      ncol = 1L,
      dimnames = list(
        names(importance_scaled),
        "MeanDecreaseAccuracy"
      )
    )
  }
  
  feature <- rownames(
    importance_unscaled
  )
  
  if (is.null(feature) ||
      !length(feature) ||
      anyNA(feature) ||
      any(!nzchar(feature))) {
    stop(
      "Importance matrix has invalid feature names: ",
      model_file
    )
  }
  
  if (!("MeanDecreaseAccuracy" %in%
        colnames(importance_unscaled))) {
    stop(
      "MeanDecreaseAccuracy is absent from: ",
      model_file
    )
  }
  
  scaled_match <- match(
    feature,
    rownames(importance_scaled)
  )
  
  if (anyNA(scaled_match)) {
    stop(
      "Scaled and unscaled importance rows differ: ",
      model_file
    )
  }
  
  importance_scaled <- importance_scaled[
    scaled_match,
    ,
    drop = FALSE
  ]
  
  result <- data.table(
    zoneID = as.integer(zoneID),
    niche = niche,
    method = method,
    method_label = method_labels[[method]],
    feature = feature,
    mda_unscaled = get_importance_column(
      importance_unscaled,
      "MeanDecreaseAccuracy"
    ),
    mda_scaled = get_importance_column(
      importance_scaled,
      "MeanDecreaseAccuracy"
    ),
    presence_mda_unscaled = get_importance_column(
      importance_unscaled,
      "1"
    ),
    presence_mda_scaled = get_importance_column(
      importance_scaled,
      "1"
    ),
    mean_decrease_gini = get_importance_column(
      importance_unscaled,
      "MeanDecreaseGini"
    ),
    model_file = model_file,
    model_object = object_name
  )
  
  if (any(!is.finite(
    result$mda_unscaled
  ))) {
    stop(
      "Non-finite MeanDecreaseAccuracy in: ",
      model_file
    )
  }
  
  result[
    ,
    `:=`(
      positive_mda = pmax(
        mda_unscaled,
        0
      ),
      importance_rank = frank(
        -mda_unscaled,
        ties.method = "min"
      )
    )
  ]
  
  positive_total <- sum(
    result$positive_mda
  )
  
  if (!is.finite(positive_total) ||
      positive_total <= 0) {
    stop(
      "No positive MDA values in: ",
      model_file
    )
  }
  
  result[
    ,
    `:=`(
      importance_share =
        positive_mda /
        positive_total,
      selected = TRUE,
      top5 = positive_mda > 0 &
        importance_rank <= 5L,
      top10 = positive_mda > 0 &
        importance_rank <= 10L
    )
  ]
  
  model_varlist <- fitted_model$varlist
  
  if (is.null(model_varlist)) {
    model_varlist <- feature
  }
  
  model_varlist <- as.character(
    model_varlist
  )
  
  importance_sd_available <-
    !is.null(fitted_model$importanceSD) &&
    any(
      is.finite(
        fitted_model$importanceSD
      )
    )
  
  audit <- data.table(
    zoneID = as.integer(zoneID),
    niche = niche,
    method = method,
    method_label = method_labels[[method]],
    model_file = model_file,
    expected_object = expected_object,
    loaded_object = object_name,
    file_exists = TRUE,
    load_ok = TRUE,
    has_importance = TRUE,
    has_importance_sd = importance_sd_available,
    has_presence_class_importance =
      "1" %in%
      colnames(importance_unscaled),
    n_predictors = nrow(result),
    n_trees = as.integer(
      fitted_model$ntree
    ),
    n_negative_mda = sum(
      result$mda_unscaled < 0
    ),
    positive_mda_total = positive_total,
    normalized_share_sum = sum(
      result$importance_share
    ),
    varlist_matches_importance = setequal(
      model_varlist,
      feature
    ),
    error_message = NA_character_
  )
  
  list(
    importance = result,
    audit = audit
  )
}


# 2. Expected model inventory ===================================================

model_specs <- rbindlist(
  list(
    data.table(
      zoneID = model_zoneID,
      niche = "climate",
      method = "rf_var",
      method_label = method_labels[["rf_var"]],
      model_file = file.path(
        climate_model_dir,
        paste0(
          "clm_rfVar_zone",
          model_zoneID,
          ".Rdata"
        )
      ),
      expected_object = "clm_rfVar"
    ),
    data.table(
      zoneID = model_zoneID,
      niche = "climate",
      method = "mf_var",
      method_label = method_labels[["mf_var"]],
      model_file = file.path(
        climate_model_dir,
        paste0(
          "clm_mfVar_zone",
          model_zoneID,
          ".Rdata"
        )
      ),
      expected_object = "clm_mfVar"
    ),
    data.table(
      zoneID = model_zoneID,
      niche = "soil",
      method = "rf_var",
      method_label = method_labels[["rf_var"]],
      model_file = file.path(
        soil_model_dir,
        paste0(
          "soil_plain_zone",
          model_zoneID,
          ".Rdata"
        )
      ),
      expected_object = "soil_plain"
    ),
    data.table(
      zoneID = model_zoneID,
      niche = "soil",
      method = "mf_var",
      method_label = method_labels[["mf_var"]],
      model_file = file.path(
        soil_model_dir,
        paste0(
          "soil_mf_zone",
          model_zoneID,
          ".Rdata"
        )
      ),
      expected_object = "soil_mf"
    )
  ),
  fill = TRUE
)

setorder(
  model_specs,
  niche,
  method,
  zoneID
)

stopifnot(
  nrow(model_specs) ==
    length(model_zoneID) *
    length(niche_order) *
    length(method_order),
  !anyDuplicated(
    model_specs[
      ,
      .(
        zoneID,
        niche,
        method
      )
    ]
  )
)


# 3. Extract stored importance ==================================================

importance_results <- list()
model_audit_results <- list()

for (model_index in seq_len(
  nrow(model_specs)
)) {
  specification <- model_specs[
    model_index
  ]
  
  cat(
    "[IMPORTANCE] ",
    specification$niche,
    " | ",
    specification$method,
    " | zone ",
    specification$zoneID,
    "\n",
    sep = ""
  )
  
  if (!file.exists(
    specification$model_file
  )) {
    model_audit_results[[
      length(model_audit_results) + 1L
    ]] <- data.table(
      zoneID = specification$zoneID,
      niche = specification$niche,
      method = specification$method,
      method_label = specification$method_label,
      model_file = specification$model_file,
      expected_object = specification$expected_object,
      loaded_object = NA_character_,
      file_exists = FALSE,
      load_ok = FALSE,
      has_importance = FALSE,
      has_importance_sd = FALSE,
      has_presence_class_importance = FALSE,
      n_predictors = NA_integer_,
      n_trees = NA_integer_,
      n_negative_mda = NA_integer_,
      positive_mda_total = NA_real_,
      normalized_share_sum = NA_real_,
      varlist_matches_importance = FALSE,
      error_message = "model file not found"
    )
    
    next
  }
  
  extracted <- tryCatch(
    extract_model_importance(
      model_file = specification$model_file,
      expected_object = specification$expected_object,
      zoneID = specification$zoneID,
      niche = specification$niche,
      method = specification$method
    ),
    error = function(error_condition) {
      list(
        importance = NULL,
        audit = data.table(
          zoneID = specification$zoneID,
          niche = specification$niche,
          method = specification$method,
          method_label = specification$method_label,
          model_file = specification$model_file,
          expected_object = specification$expected_object,
          loaded_object = NA_character_,
          file_exists = TRUE,
          load_ok = FALSE,
          has_importance = FALSE,
          has_importance_sd = FALSE,
          has_presence_class_importance = FALSE,
          n_predictors = NA_integer_,
          n_trees = NA_integer_,
          n_negative_mda = NA_integer_,
          positive_mda_total = NA_real_,
          normalized_share_sum = NA_real_,
          varlist_matches_importance = FALSE,
          error_message = conditionMessage(
            error_condition
          )
        )
      )
    }
  )
  
  model_audit_results[[
    length(model_audit_results) + 1L
  ]] <- extracted$audit
  
  if (!is.null(extracted$importance)) {
    importance_results[[
      length(importance_results) + 1L
    ]] <- extracted$importance
  }
}

model_audit <- rbindlist(
  model_audit_results,
  fill = TRUE
)

setorder(
  model_audit,
  niche,
  method,
  zoneID
)

fwrite(
  model_audit,
  file.path(
    table_dir,
    "FI_01_model_importance_audit.csv"
  )
)

failed_models <- model_audit[
  load_ok != TRUE |
    has_importance != TRUE |
    varlist_matches_importance != TRUE |
    !is.finite(normalized_share_sum) |
    abs(normalized_share_sum - 1) > 1e-10
]

if (nrow(failed_models) > 0L) {
  print(
    failed_models
  )
  
  if (strict_complete_models) {
    stop(
      "Feature-importance extraction failed for ",
      nrow(failed_models),
      " expected model(s). See FI_01_model_importance_audit.csv."
    )
  }
}

if (!length(importance_results)) {
  stop(
    "No model importance was extracted."
  )
}

importance_long <- rbindlist(
  importance_results,
  fill = TRUE
)

feature_annotation <- unique(
  annotate_features(
    importance_long$niche,
    importance_long$feature
  )
)

importance_long <- merge(
  importance_long,
  feature_annotation,
  by = c(
    "niche",
    "feature"
  ),
  all.x = TRUE,
  sort = FALSE
)

setorder(
  importance_long,
  niche,
  method,
  zoneID,
  importance_rank,
  feature
)

fwrite(
  importance_long,
  file.path(
    table_dir,
    "FI_02_zone_model_feature_importance_long.csv"
  )
)


# 4. Zone names, categories, and observed reference area ========================

palette <- fread(
  require_file(
    palette_file
  )
)

required_palette_columns <- c(
  "zoneID",
  "zone",
  "category",
  "category2"
)

missing_palette_columns <- setdiff(
  required_palette_columns,
  names(palette)
)

if (length(missing_palette_columns) > 0L) {
  stop(
    "Palette is missing columns: ",
    paste(
      missing_palette_columns,
      collapse = ", "
    )
  )
}

palette[
  ,
  zoneID := as.integer(zoneID)
]

palette <- palette[
  zoneID %in% model_zoneID,
  .(
    zoneID,
    zone_name = as.character(zone),
    detailed_category = as.character(category),
    category2 = as.character(category2),
    zone_color = if (
      "COLOR" %in% names(palette)
    ) {
      as.character(COLOR)
    } else {
      NA_character_
    }
  )
]

if (nrow(palette) !=
    length(model_zoneID) ||
    anyDuplicated(palette$zoneID) ||
    length(setdiff(
      model_zoneID,
      palette$zoneID
    )) > 0L) {
  stop(
    "Palette does not contain exactly one row for every modeled zone."
  )
}

reference_zone <- rast(
  require_file(
    reference_raster_file
  )
)[[1]]

names(reference_zone) <- "zoneID"

reference_cell_area <- cellSize(
  reference_zone,
  unit = "km"
)

names(reference_cell_area) <-
  "reference_area_km2"

reference_area <- as.data.table(
  zonal(
    reference_cell_area,
    reference_zone,
    fun = "sum",
    na.rm = TRUE
  )
)

if (ncol(reference_area) < 2L) {
  stop(
    "Could not calculate reference area by zone."
  )
}

setnames(
  reference_area,
  names(reference_area)[1:2],
  c(
    "zoneID",
    "reference_area_km2"
  )
)

reference_area <- reference_area[
  ,
  .(
    zoneID = as.integer(zoneID),
    reference_area_km2 = as.numeric(
      reference_area_km2
    )
  )
]

reference_frequency <- as.data.table(
  freq(
    reference_zone,
    bylayer = FALSE
  )
)

frequency_value_column <- intersect(
  c(
    "value",
    "zoneID"
  ),
  names(reference_frequency)
)

frequency_count_column <- intersect(
  c(
    "count",
    "frequency"
  ),
  names(reference_frequency)
)

if (length(frequency_value_column) != 1L ||
    length(frequency_count_column) != 1L) {
  stop(
    "Could not identify value/count columns returned by terra::freq."
  )
}

reference_frequency <- reference_frequency[
  !is.na(
    get(
      frequency_value_column
    )
  )
]

reference_frequency <- reference_frequency[
  ,
  .(
    zoneID = as.integer(
      get(
        frequency_value_column
      )
    ),
    reference_pixels = as.numeric(
      get(
        frequency_count_column
      )
    )
  )
]

zone_metadata <- merge(
  palette,
  reference_area,
  by = "zoneID",
  all.x = TRUE,
  sort = FALSE
)

zone_metadata <- merge(
  zone_metadata,
  reference_frequency,
  by = "zoneID",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(
  zone_metadata[
    ,
    .(
      reference_area_km2,
      reference_pixels
    )
  ]
)) {
  stop(
    "Observed reference area or pixel count is missing for a modeled zone."
  )
}

zone_metadata[
  ,
  `:=`(
    category_label = gsub(
      "_",
      " ",
      category2,
      fixed = TRUE
    ),
    global_area_rank_largest = frank(
      -reference_area_km2,
      ties.method = "min"
    ),
    global_area_rank_smallest = frank(
      reference_area_km2,
      ties.method = "min"
    )
  )
]

zone_metadata[
  ,
  `:=`(
    category_area_rank_largest = frank(
      -reference_area_km2,
      ties.method = "min"
    ),
    category_area_rank_smallest = frank(
      reference_area_km2,
      ties.method = "min"
    )
  ),
  by = category2
]

setorder(
  zone_metadata,
  zoneID
)

fwrite(
  zone_metadata,
  file.path(
    table_dir,
    "FI_03_zone_metadata_and_reference_area.csv"
  )
)


# 5. RF versus MF RF agreement and consensus profiles ===========================

agreement_grid_list <- list()
agreement_grid_index <- 0L

for (zone_value in model_zoneID) {
  for (niche_value in niche_order) {
    subset_importance <- importance_long[
      zoneID == zone_value &
        niche == niche_value
    ]
    
    feature_universe <- sort(
      unique(
        subset_importance$feature
      )
    )
    
    agreement_grid_index <-
      agreement_grid_index + 1L
    
    agreement_grid <- CJ(
      method = method_order,
      feature = feature_universe,
      unique = TRUE
    )
    
    agreement_grid[
      ,
      `:=`(
        zoneID = zone_value,
        niche = niche_value
      )
    ]
    
    agreement_grid <- merge(
      agreement_grid,
      subset_importance[
        ,
        .(
          zoneID,
          niche,
          method,
          feature,
          importance_share,
          mda_unscaled,
          importance_rank,
          selected,
          top5
        )
      ],
      by = c(
        "zoneID",
        "niche",
        "method",
        "feature"
      ),
      all.x = TRUE,
      sort = FALSE
    )
    
    agreement_grid[
      is.na(importance_share),
      `:=`(
        importance_share = 0,
        mda_unscaled = 0,
        selected = FALSE,
        top5 = FALSE
      )
    ]
    
    agreement_grid_list[[
      agreement_grid_index
    ]] <- agreement_grid
  }
}

agreement_grid <- rbindlist(
  agreement_grid_list,
  fill = TRUE
)

importance_share_wide <- dcast(
  agreement_grid,
  zoneID + niche + feature ~ method,
  value.var = "importance_share",
  fill = 0
)

mda_wide <- dcast(
  agreement_grid,
  zoneID + niche + feature ~ method,
  value.var = "mda_unscaled",
  fill = 0
)

top5_wide <- dcast(
  agreement_grid,
  zoneID + niche + feature ~ method,
  value.var = "top5",
  fill = FALSE
)

selected_wide <- dcast(
  agreement_grid,
  zoneID + niche + feature ~ method,
  value.var = "selected",
  fill = FALSE
)

method_agreement <- importance_share_wide[
  ,
  .(
    importance_spearman = safe_cor(
      rf_var,
      mf_var,
      method = "spearman"
    ),
    importance_pearson = safe_cor(
      rf_var,
      mf_var,
      method = "pearson"
    )
  ),
  by = .(
    zoneID,
    niche
  )
]

rank_agreement <- top5_wide[
  ,
  {
    rf_top5 <- rf_var
    mf_top5 <- mf_var
    
    top5_union <- sum(
      rf_top5 |
        mf_top5
    )
    
    .(
      top5_overlap_n = sum(
        rf_top5 &
          mf_top5
      ),
      top5_union_n = top5_union,
      top5_jaccard = if (
        top5_union > 0
      ) {
        sum(
          rf_top5 &
            mf_top5
        ) /
          top5_union
      } else {
        NA_real_
      }
    )
  },
  by = .(
    zoneID,
    niche
  )
]

selected_agreement <- selected_wide[
  ,
  {
    selected_union <- sum(
      rf_var |
        mf_var
    )
    
    .(
      selected_overlap_n = sum(
        rf_var &
          mf_var
      ),
      selected_union_n = selected_union,
      selected_jaccard = if (
        selected_union > 0
      ) {
        sum(
          rf_var &
            mf_var
        ) /
          selected_union
      } else {
        NA_real_
      }
    )
  },
  by = .(
    zoneID,
    niche
  )
]

method_agreement <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = c(
        "zoneID",
        "niche"
      ),
      all = TRUE,
      sort = FALSE
    )
  },
  list(
    method_agreement,
    rank_agreement,
    selected_agreement
  )
)

method_agreement <- merge(
  method_agreement,
  zone_metadata,
  by = "zoneID",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  method_agreement,
  niche,
  category2,
  zoneID
)

fwrite(
  method_agreement,
  file.path(
    table_dir,
    "FI_04_rf_vs_mf_importance_agreement.csv"
  )
)

consensus_profile <- importance_share_wide[
  ,
  .(
    zoneID,
    niche,
    feature,
    rf_importance_share = rf_var,
    mf_importance_share = mf_var,
    consensus_importance_share =
      (
        rf_var +
          mf_var
      ) /
      2
  )
]

consensus_mda <- mda_wide[
  ,
  .(
    zoneID,
    niche,
    feature,
    rf_mda_unscaled = rf_var,
    mf_mda_unscaled = mf_var,
    mean_mda_unscaled =
      (
        rf_var +
          mf_var
      ) /
      2
  )
]

consensus_selected <- selected_wide[
  ,
  .(
    zoneID,
    niche,
    feature,
    selected_in_rf = rf_var,
    selected_in_mf = mf_var,
    selected_method_count =
      as.integer(rf_var) +
      as.integer(mf_var)
  )
]

consensus_profile <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = c(
        "zoneID",
        "niche",
        "feature"
      ),
      all = TRUE,
      sort = FALSE
    )
  },
  list(
    consensus_profile,
    consensus_mda,
    consensus_selected
  )
)

consensus_profile <- merge(
  consensus_profile,
  feature_annotation,
  by = c(
    "niche",
    "feature"
  ),
  all.x = TRUE,
  sort = FALSE
)

consensus_profile[
  ,
  `:=`(
    selected = selected_method_count > 0L,
    consensus_rank = frank(
      -consensus_importance_share,
      ties.method = "min"
    )
  ),
  by = .(
    zoneID,
    niche
  )
]

consensus_profile[
  ,
  `:=`(
    top5 = selected &
      consensus_importance_share > 0 &
      consensus_rank <= 5L,
    top10 = selected &
      consensus_importance_share > 0 &
      consensus_rank <= 10L
  )
]

consensus_share_check <- consensus_profile[
  ,
  .(
    share_sum = sum(
      consensus_importance_share
    )
  ),
  by = .(
    zoneID,
    niche
  )
]

if (any(
  abs(
    consensus_share_check$share_sum -
    1
  ) > 1e-10
)) {
  print(
    consensus_share_check[
      abs(share_sum - 1) > 1e-10
    ]
  )
  
  stop(
    "Consensus importance shares do not sum to one."
  )
}


# 6. Complete zone profiles and category-level common factors ===================

profile_grid_list <- lapply(
  niche_order,
  function(niche_value) {
    feature_universe <- sort(
      unique(
        consensus_profile[
          niche == niche_value,
          feature
        ]
      )
    )
    
    grid <- CJ(
      zoneID = model_zoneID,
      feature = feature_universe,
      unique = TRUE
    )
    
    grid[
      ,
      niche := niche_value
    ]
    
    grid
  }
)

zone_profile <- rbindlist(
  profile_grid_list,
  fill = TRUE
)

zone_profile <- merge(
  zone_profile,
  consensus_profile,
  by = c(
    "zoneID",
    "niche",
    "feature"
  ),
  all.x = TRUE,
  sort = FALSE
)

zone_profile[
  is.na(consensus_importance_share),
  `:=`(
    rf_importance_share = 0,
    mf_importance_share = 0,
    consensus_importance_share = 0,
    rf_mda_unscaled = 0,
    mf_mda_unscaled = 0,
    mean_mda_unscaled = 0,
    selected_in_rf = FALSE,
    selected_in_mf = FALSE,
    selected_method_count = 0L,
    selected = FALSE,
    consensus_rank = NA_real_,
    top5 = FALSE,
    top10 = FALSE
  )
]

missing_annotation <- is.na(
  zone_profile$feature_group
)

if (any(missing_annotation)) {
  new_annotation <- annotate_features(
    zone_profile$niche[missing_annotation],
    zone_profile$feature[missing_annotation]
  )
  
  zone_profile$feature_group[missing_annotation] <-
    new_annotation$feature_group
  
  zone_profile$climate_period[missing_annotation] <-
    new_annotation$climate_period
}

zone_profile <- merge(
  zone_profile,
  zone_metadata,
  by = "zoneID",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  zone_profile,
  niche,
  category2,
  zoneID,
  -consensus_importance_share,
  feature
)

fwrite(
  zone_profile,
  file.path(
    table_dir,
    "FI_05_consensus_zone_feature_profiles_complete.csv"
  )
)

zone_top_factors <- zone_profile[
  selected == TRUE &
    consensus_importance_share > 0 &
    consensus_rank <= top_n_zone_factors
]

setorder(
  zone_top_factors,
  niche,
  category2,
  zoneID,
  consensus_rank,
  feature
)

fwrite(
  zone_top_factors,
  file.path(
    table_dir,
    "FI_06_zone_top10_consensus_factors.csv"
  )
)

category_variable_summary <- zone_profile[
  ,
  .(
    n_zones = uniqueN(zoneID),
    mean_importance_share = mean(
      consensus_importance_share
    ),
    sd_importance_share = sd(
      consensus_importance_share
    ),
    median_importance_share = median(
      consensus_importance_share
    ),
    selected_zones = sum(
      selected
    ),
    top5_zones = sum(
      top5
    ),
    top10_zones = sum(
      top10
    )
  ),
  by = .(
    category2,
    category_label,
    niche,
    feature,
    feature_group,
    climate_period
  )
]

category_variable_summary[
  ,
  `:=`(
    se_importance_share =
      sd_importance_share /
      sqrt(n_zones),
    selection_frequency =
      selected_zones /
      n_zones,
    top5_frequency =
      top5_zones /
      n_zones,
    top10_frequency =
      top10_zones /
      n_zones,
    category_variable_rank = frank(
      -mean_importance_share,
      ties.method = "min"
    ),
    inference_scope = ifelse(
      n_zones >=
        minimum_zones_for_category_inference,
      "within-category pattern",
      "descriptive only"
    )
  ),
  by = .(
    category2,
    niche
  )
]

category_variable_summary[
  ,
  common_factor :=
    n_zones >=
    minimum_zones_for_category_inference &
    selection_frequency >=
    common_selection_frequency &
    top5_frequency >=
    common_top5_frequency
]

setorder(
  category_variable_summary,
  niche,
  category2,
  category_variable_rank,
  feature
)

fwrite(
  category_variable_summary,
  file.path(
    table_dir,
    "FI_07_category_variable_summary.csv"
  )
)

category_common_factors <- category_variable_summary[
  category_variable_rank <=
    top_n_category_factors |
    common_factor == TRUE
]

fwrite(
  category_common_factors,
  file.path(
    table_dir,
    "FI_08_category_common_and_top_factors.csv"
  )
)

overall_variable_summary <- zone_profile[
  ,
  .(
    overall_mean_importance_share = mean(
      consensus_importance_share
    ),
    overall_selection_frequency = mean(
      selected
    ),
    overall_top5_frequency = mean(
      top5
    )
  ),
  by = .(
    niche,
    feature
  )
]

category_contrast <- merge(
  category_variable_summary,
  overall_variable_summary,
  by = c(
    "niche",
    "feature"
  ),
  all.x = TRUE,
  sort = FALSE
)

category_contrast[
  ,
  `:=`(
    importance_difference_from_all_zones =
      mean_importance_share -
      overall_mean_importance_share,
    importance_ratio_to_all_zones = fifelse(
      overall_mean_importance_share > 0,
      mean_importance_share /
        overall_mean_importance_share,
      NA_real_
    ),
    selection_frequency_difference =
      selection_frequency -
      overall_selection_frequency,
    top5_frequency_difference =
      top5_frequency -
      overall_top5_frequency
  )
]

setorder(
  category_contrast,
  niche,
  category2,
  -importance_difference_from_all_zones,
  feature
)

fwrite(
  category_contrast,
  file.path(
    table_dir,
    "FI_09_category_contrasts_vs_all_zones.csv"
  )
)


# 7. Ecological factor-group and seasonal summaries =============================

zone_factor_group <- zone_profile[
  ,
  .(
    group_importance_share = sum(
      consensus_importance_share
    )
  ),
  by = .(
    zoneID,
    zone_name,
    category2,
    category_label,
    niche,
    feature_group,
    reference_area_km2
  )
]

category_factor_group <- zone_factor_group[
  ,
  .(
    n_zones = uniqueN(zoneID),
    mean_group_importance_share = mean(
      group_importance_share
    ),
    sd_group_importance_share = sd(
      group_importance_share
    ),
    median_group_importance_share = median(
      group_importance_share
    )
  ),
  by = .(
    category2,
    category_label,
    niche,
    feature_group
  )
]

category_factor_group[
  ,
  `:=`(
    se_group_importance_share =
      sd_group_importance_share /
      sqrt(n_zones),
    group_rank = frank(
      -mean_group_importance_share,
      ties.method = "min"
    ),
    inference_scope = ifelse(
      n_zones >=
        minimum_zones_for_category_inference,
      "within-category pattern",
      "descriptive only"
    )
  ),
  by = .(
    category2,
    niche
  )
]

setorder(
  category_factor_group,
  niche,
  category2,
  group_rank,
  feature_group
)

fwrite(
  category_factor_group,
  file.path(
    table_dir,
    "FI_10_category_ecological_factor_groups.csv"
  )
)

zone_climate_period <- zone_profile[
  niche == "climate" &
    !is.na(climate_period),
  .(
    period_importance_share = sum(
      consensus_importance_share
    )
  ),
  by = .(
    zoneID,
    zone_name,
    category2,
    category_label,
    climate_period,
    reference_area_km2
  )
]

category_climate_period <- zone_climate_period[
  ,
  .(
    n_zones = uniqueN(zoneID),
    mean_period_importance_share = mean(
      period_importance_share
    ),
    sd_period_importance_share = sd(
      period_importance_share
    ),
    median_period_importance_share = median(
      period_importance_share
    )
  ),
  by = .(
    category2,
    category_label,
    climate_period
  )
]

category_climate_period[
  ,
  se_period_importance_share :=
    sd_period_importance_share /
    sqrt(n_zones)
]

fwrite(
  category_climate_period,
  file.path(
    table_dir,
    "FI_11_category_climate_period_summary.csv"
  )
)

theory_review_shortlist <- category_factor_group[
  group_rank <= top_n_theory_groups,
  .(
    category2,
    category_label,
    niche,
    n_zones,
    inference_scope,
    group_rank,
    feature_group,
    mean_group_importance_share,
    se_group_importance_share
  )
]

top_variable_strings <- category_variable_summary[
  category_variable_rank <=
    top_n_category_factors,
  .(
    top_variables = paste0(
      feature,
      " [",
      round(
        100 *
          mean_importance_share,
        1
      ),
      "%; selected ",
      round(
        100 *
          selection_frequency
      ),
      "%]",
      collapse = "; "
    )
  ),
  by = .(
    category2,
    niche
  )
]

theory_review_shortlist <- merge(
  theory_review_shortlist,
  top_variable_strings,
  by = c(
    "category2",
    "niche"
  ),
  all.x = TRUE,
  sort = FALSE
)

setorder(
  theory_review_shortlist,
  niche,
  category2,
  group_rank
)

fwrite(
  theory_review_shortlist,
  file.path(
    table_dir,
    "FI_12_theory_review_shortlist.csv"
  )
)


# 8. Similarity among broad categories ==========================================

category_profile <- category_variable_summary[
  ,
  .(
    category2,
    category_label,
    niche,
    feature,
    mean_importance_share
  )
]

category_similarity_list <- list()

for (niche_value in niche_order) {
  profile_wide <- dcast(
    category_profile[
      niche == niche_value
    ],
    category2 + category_label ~ feature,
    value.var = "mean_importance_share",
    fill = 0
  )
  
  category_columns <- c(
    "category2",
    "category_label"
  )
  
  profile_matrix <- as.matrix(
    profile_wide[
      ,
      setdiff(
        names(profile_wide),
        category_columns
      ),
      with = FALSE
    ]
  )
  
  categories <- profile_wide$category2
  category_labels <- profile_wide$category_label
  
  for (first_index in seq_along(categories)) {
    for (second_index in seq_along(categories)) {
      category_similarity_list[[
        length(category_similarity_list) + 1L
      ]] <- data.table(
        niche = niche_value,
        category1 = categories[first_index],
        category1_label = category_labels[first_index],
        category2 = categories[second_index],
        category2_label = category_labels[second_index],
        cosine_similarity = cosine_similarity(
          profile_matrix[first_index, ],
          profile_matrix[second_index, ]
        ),
        spearman_similarity = safe_cor(
          profile_matrix[first_index, ],
          profile_matrix[second_index, ],
          method = "spearman"
        )
      )
    }
  }
}

category_similarity <- rbindlist(
  category_similarity_list,
  fill = TRUE
)

fwrite(
  category_similarity,
  file.path(
    table_dir,
    "FI_13_category_importance_profile_similarity.csv"
  )
)


# 9. Category medoids and observed-area extremes ================================

joint_profile <- zone_profile[
  ,
  .(
    category2,
    category_label,
    joint_feature = paste(
      niche,
      feature,
      sep = "::"
    ),
    joint_importance_share =
      consensus_importance_share /
      uniqueN(niche)
  ),
  by = zoneID
]

category_medoid_list <- list()

for (category_value in sort(
  unique(zone_metadata$category2)
)) {
  category_joint <- joint_profile[
    category2 == category_value
  ]
  
  category_joint_wide <- dcast(
    category_joint,
    zoneID ~ joint_feature,
    value.var = "joint_importance_share",
    fill = 0
  )
  
  category_zoneID <- category_joint_wide$zoneID
  
  category_matrix <- as.matrix(
    category_joint_wide[
      ,
      setdiff(
        names(category_joint_wide),
        "zoneID"
      ),
      with = FALSE
    ]
  )
  
  profile_norm <- sqrt(
    rowSums(
      category_matrix ^ 2
    )
  )
  
  if (any(
    !is.finite(profile_norm) |
    profile_norm <= 0
  )) {
    stop(
      "A category contains an invalid joint importance profile: ",
      category_value
    )
  }
  
  pairwise_similarity <-
    tcrossprod(category_matrix) /
    outer(
      profile_norm,
      profile_norm
    )
  
  mean_similarity <- if (
    nrow(pairwise_similarity) > 1L
  ) {
    (
      rowSums(pairwise_similarity) -
        diag(pairwise_similarity)
    ) /
      (
        nrow(pairwise_similarity) - 1L
      )
  } else {
    1
  }
  
  medoid_index <- which.max(
    mean_similarity
  )
  
  category_medoid_list[[
    length(category_medoid_list) + 1L
  ]] <- data.table(
    category2 = category_value,
    zoneID = category_zoneID[medoid_index],
    category_profile_similarity = mean_similarity[medoid_index],
    category_n_zones = length(category_zoneID)
  )
}

category_medoids <- rbindlist(
  category_medoid_list,
  fill = TRUE
)

category_medoids[
  ,
  selection_role := "category profile medoid"
]

category_largest <- zone_metadata[
  ,
  .SD[
    which.max(
      reference_area_km2
    )
  ],
  by = category2
][
  ,
  selection_role := "category largest observed area"
]

category_smallest <- zone_metadata[
  ,
  .SD[
    which.min(
      reference_area_km2
    )
  ],
  by = category2
][
  ,
  selection_role := "category smallest observed area"
]

global_largest <- head(
  zone_metadata[
    order(
      -reference_area_km2
    )
  ],
  n_global_area_extremes
)[
  ,
  selection_role := "global largest observed area"
]

global_smallest <- head(
  zone_metadata[
    order(
      reference_area_km2
    )
  ],
  n_global_area_extremes
)[
  ,
  selection_role := "global smallest observed area"
]

zone_selection_long <- rbindlist(
  list(
    category_medoids[
      ,
      .(
        category2,
        zoneID,
        selection_role,
        category_profile_similarity,
        category_n_zones
      )
    ],
    category_largest[
      ,
      .(
        category2,
        zoneID,
        selection_role,
        category_profile_similarity = NA_real_,
        category_n_zones = NA_integer_
      )
    ],
    category_smallest[
      ,
      .(
        category2,
        zoneID,
        selection_role,
        category_profile_similarity = NA_real_,
        category_n_zones = NA_integer_
      )
    ],
    global_largest[
      ,
      .(
        category2,
        zoneID,
        selection_role,
        category_profile_similarity = NA_real_,
        category_n_zones = NA_integer_
      )
    ],
    global_smallest[
      ,
      .(
        category2,
        zoneID,
        selection_role,
        category_profile_similarity = NA_real_,
        category_n_zones = NA_integer_
      )
    ]
  ),
  fill = TRUE
)

zone_selection <- zone_selection_long[
  ,
  .(
    selection_roles = paste(
      sort(
        unique(selection_role)
      ),
      collapse = "; "
    ),
    category_profile_similarity = {
      values <- category_profile_similarity[
        is.finite(category_profile_similarity)
      ]
      
      if (length(values) > 0L) {
        as.numeric(max(values))
      } else {
        NA_real_
      }
    },
    category_n_zones = {
      values <- category_n_zones[
        !is.na(category_n_zones)
      ]
      
      if (length(values) > 0L) {
        as.integer(max(values))
      } else {
        NA_integer_
      }
    }
  ),
  by = .(
    zoneID,
    category2
  )
]

zone_selection <- merge(
  zone_selection,
  zone_metadata,
  by = c(
    "zoneID",
    "category2"
  ),
  all.x = TRUE,
  sort = FALSE
)

setorder(
  zone_selection,
  category2,
  -reference_area_km2,
  zoneID
)

fwrite(
  zone_selection,
  file.path(
    table_dir,
    "FI_14_representative_and_area_extreme_zones.csv"
  )
)

selected_zone_top_factors <- merge(
  zone_top_factors,
  zone_selection[
    ,
    .(
      zoneID,
      selection_roles,
      category_profile_similarity
    )
  ],
  by = "zoneID",
  all = FALSE,
  sort = FALSE
)

setorder(
  selected_zone_top_factors,
  category2,
  zoneID,
  niche,
  consensus_rank,
  feature
)

fwrite(
  selected_zone_top_factors,
  file.path(
    table_dir,
    "FI_15_selected_zone_top_factors.csv"
  )
)


# 10. Zone area and concentration of importance =================================

zone_importance_concentration <- zone_profile[
  ,
  {
    positive_share <- consensus_importance_share[
      consensus_importance_share > 0
    ]
    
    shannon <- -sum(
      positive_share *
        log(positive_share)
    )
    
    .(
      n_selected_features = sum(
        selected
      ),
      top1_importance_share = max(
        consensus_importance_share
      ),
      top5_importance_share = sum(
        sort(
          consensus_importance_share,
          decreasing = TRUE
        )[
          seq_len(
            min(
              5L,
              length(
                consensus_importance_share
              )
            )
          )
        ]
      ),
      importance_shannon = shannon,
      effective_number_of_features = exp(
        shannon
      )
    )
  },
  by = .(
    zoneID,
    zone_name,
    category2,
    category_label,
    niche,
    reference_area_km2,
    reference_pixels
  )
]

zone_importance_concentration <- merge(
  zone_importance_concentration,
  zone_selection[
    ,
    .(
      zoneID,
      selection_roles
    )
  ],
  by = "zoneID",
  all.x = TRUE,
  sort = FALSE
)

zone_importance_concentration[
  ,
  selected_for_detail :=
    !is.na(selection_roles)
]

fwrite(
  zone_importance_concentration,
  file.path(
    table_dir,
    "FI_16_zone_area_and_importance_concentration.csv"
  )
)

area_concentration_correlations <- rbindlist(
  list(
    zone_importance_concentration[
      ,
      .(
        scope = "all modeled zones",
        n_zones = uniqueN(zoneID),
        area_vs_effective_n_spearman = safe_cor(
          reference_area_km2,
          effective_number_of_features
        ),
        area_vs_top1_share_spearman = safe_cor(
          reference_area_km2,
          top1_importance_share
        ),
        area_vs_top5_share_spearman = safe_cor(
          reference_area_km2,
          top5_importance_share
        )
      ),
      by = niche
    ],
    zone_importance_concentration[
      ,
      .(
        scope = "within broad category",
        n_zones = uniqueN(zoneID),
        area_vs_effective_n_spearman = safe_cor(
          reference_area_km2,
          effective_number_of_features
        ),
        area_vs_top1_share_spearman = safe_cor(
          reference_area_km2,
          top1_importance_share
        ),
        area_vs_top5_share_spearman = safe_cor(
          reference_area_km2,
          top5_importance_share
        )
      ),
      by = .(
        niche,
        category2
      )
    ]
  ),
  fill = TRUE
)

fwrite(
  area_concentration_correlations,
  file.path(
    table_dir,
    "FI_17_area_importance_concentration_correlations.csv"
  )
)


# 11. Figures: common factors within categories =================================

category_top_for_plot <- category_variable_summary[
  category_variable_rank <=
    top_n_category_factors
]

category_top_for_plot[
  ,
  `:=`(
    lower_share = pmax(
      0,
      mean_importance_share -
        se_importance_share
    ),
    upper_share =
      mean_importance_share +
      se_importance_share
  )
]

for (niche_value in niche_order) {
  figure_data <- copy(
    category_top_for_plot[
      niche == niche_value
    ]
  )
  
  setorder(
    figure_data,
    category_label,
    mean_importance_share,
    feature
  )
  
  figure_data[
    ,
    feature_panel := factor(
      paste(
        category2,
        feature,
        sep = "|||"
      ),
      levels = unique(
        paste(
          category2,
          feature,
          sep = "|||"
        )
      )
    )
  ]
  
  category_factor_figure <- ggplot(
    figure_data,
    aes(
      x = mean_importance_share,
      y = feature_panel
    )
  ) +
    geom_errorbarh(
      aes(
        xmin = lower_share,
        xmax = upper_share
      ),
      height = 0.16,
      linewidth = 0.45,
      color = "grey45"
    ) +
    geom_point(
      aes(
        size = selection_frequency
      ),
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.5
    ) +
    facet_wrap(
      ~ category_label,
      scales = "free_y",
      ncol = 2
    ) +
    scale_y_discrete(
      labels = function(x) {
        sub(
          "^.*\\|\\|\\|",
          "",
          x
        )
      }
    ) +
    scale_x_continuous(
      labels = percent_label,
      expand = expansion(
        mult = c(
          0,
          0.08
        )
      )
    ) +
    scale_size_continuous(
      range = c(
        1.8,
        4.2
      ),
      limits = c(
        0,
        1
      ),
      labels = percent_label
    ) +
    labs(
      x = "Mean normalized permutation importance",
      y = NULL,
      size = "Zones selecting\nthe variable",
      title = paste0(
        "Common ",
        niche_value,
        " factors within broad vegetation categories"
      ),
      subtitle = paste(
        "Points are category means across zones; bars are +/-1 SE.",
        "Unselected variables contribute zero to the category mean."
      )
    ) +
    theme_bw(
      base_size = 9.5
    ) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text = element_text(
        face = "bold"
      ),
      panel.spacing = grid::unit(
        0.75,
        "lines"
      )
    )
  
  save_plot(
    category_factor_figure,
    file.path(
      figure_dir,
      paste0(
        "Figure_FI_1_category_common_factors_",
        niche_value,
        ".png"
      )
    ),
    10.5,
    10.5
  )
}


# 12. Figures: ecological factor-group composition ==============================

for (niche_value in niche_order) {
  group_figure_data <- copy(
    category_factor_group[
      niche == niche_value
    ]
  )
  
  group_levels <- sort(
    unique(
      group_figure_data$feature_group
    )
  )
  
  group_colors <- setNames(
    hcl.colors(
      length(group_levels),
      palette = "Dark 3"
    ),
    group_levels
  )
  
  group_figure_data[
    ,
    feature_group := factor(
      feature_group,
      levels = group_levels
    )
  ]
  
  group_composition_figure <- ggplot(
    group_figure_data,
    aes(
      x = category_label,
      y = mean_group_importance_share,
      fill = feature_group
    )
  ) +
    geom_col(
      width = 0.68,
      color = "white",
      linewidth = 0.15
    ) +
    coord_flip() +
    scale_fill_manual(
      values = group_colors
    ) +
    scale_y_continuous(
      limits = c(
        0,
        1
      ),
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      labels = percent_label,
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    labs(
      x = NULL,
      y = "Mean share of normalized importance",
      fill = NULL,
      title = paste0(
        "Category-level ",
        niche_value,
        " importance by ecological factor group"
      ),
      subtitle =
        "Each bar sums to 100%; RF and MF RF importance shares are averaged within zone first."
    ) +
    theme_bw(
      base_size = 10
    ) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )
  
  save_plot(
    group_composition_figure,
    file.path(
      figure_dir,
      paste0(
        "Figure_FI_2_category_factor_groups_",
        niche_value,
        ".png"
      )
    ),
    10.5,
    6.5
  )
}


# 13. Figure: annual versus seasonal climate contribution =======================

climate_period_levels <- c(
  "Annual",
  "Winter",
  "Spring",
  "Summer",
  "Autumn"
)

category_climate_period[
  ,
  climate_period := factor(
    climate_period,
    levels = climate_period_levels
  )
]

climate_period_colors <- c(
  Annual = "#4D4D4D",
  Winter = "#5B8DB8",
  Spring = "#78A96B",
  Summer = "#D58A45",
  Autumn = "#9A6FB0"
)

climate_period_figure <- ggplot(
  category_climate_period,
  aes(
    x = category_label,
    y = mean_period_importance_share,
    fill = climate_period
  )
) +
  geom_col(
    width = 0.68,
    color = "white",
    linewidth = 0.15
  ) +
  coord_flip() +
  scale_fill_manual(
    values = climate_period_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    labels = percent_label,
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  labs(
    x = NULL,
    y = "Mean share of normalized climate importance",
    fill = NULL,
    title = "Annual and seasonal contributions to climate-model importance",
    subtitle =
      "Season denotes the ClimateNA variable suffix, not the season in which vegetation was observed."
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

save_plot(
  climate_period_figure,
  file.path(
    figure_dir,
    "Figure_FI_3_category_climate_period_composition.png"
  ),
  9.2,
  6.2
)


# 14. Figure: similarity among categories =======================================

# Plot each unordered category pair once. A common 0--1 axis makes the climate
# and soil comparisons directly comparable without encoding exact values only
# through a heat-map colour scale.
category_similarity_plot <- category_similarity[
  category1 < category2
]

category_similarity_plot[
  ,
  pair_label := paste(
    category1_label,
    category2_label,
    sep = " -- "
  )
]

setorder(
  category_similarity_plot,
  niche,
  cosine_similarity,
  pair_label
)

category_similarity_plot[
  ,
  `:=`(
    pair_panel = factor(
      paste(
        niche,
        pair_label,
        sep = "|||"
      ),
      levels = unique(
        paste(
          niche,
          pair_label,
          sep = "|||"
        )
      )
    ),
    niche = factor(
      niche,
      levels = niche_order
    )
  )
]

category_similarity_figure <- ggplot(
  category_similarity_plot,
  aes(
    x = cosine_similarity,
    y = pair_panel
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = cosine_similarity,
      yend = pair_panel
    ),
    linewidth = 0.45,
    color = "grey72"
  ) +
  geom_point(
    shape = 21,
    size = 2.2,
    stroke = 0.45,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(
    ~ niche,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub(
        "^.*\\|\\|\\|",
        "",
        x
      )
    }
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  labs(
    x = "Cosine similarity",
    y = NULL,
    title = "Pairwise similarity of category-level feature-importance profiles",
    subtitle =
      "Complete consensus profiles include zeros for variables not selected in a zone."
  ) +
  theme_bw(
    base_size = 8.8
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    panel.spacing = grid::unit(
      1,
      "lines"
    )
  )

save_plot(
  category_similarity_figure,
  file.path(
    figure_dir,
    "Figure_FI_4_category_pairwise_importance_similarity.png"
  ),
  13.5,
  10.0
)


# 15. Figure: RF versus MF RF robustness ========================================

method_agreement_long <- rbindlist(
  list(
    method_agreement[
      ,
      .(
        zoneID,
        zone_name,
        category2,
        category_label,
        niche,
        agreement_metric = "Spearman correlation",
        agreement_value = importance_spearman
      )
    ],
    method_agreement[
      ,
      .(
        zoneID,
        zone_name,
        category2,
        category_label,
        niche,
        agreement_metric = "Top-5 Jaccard overlap",
        agreement_value = top5_jaccard
      )
    ]
  ),
  fill = TRUE
)

method_agreement_figure <- ggplot(
  method_agreement_long,
  aes(
    x = category_label,
    y = agreement_value
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(
    position = position_jitter(
      width = 0.12,
      height = 0,
      seed = 49
    ),
    shape = 21,
    fill = "grey75",
    color = "black",
    size = 1.6,
    stroke = 0.35,
    alpha = 0.75
  ) +
  facet_grid(
    agreement_metric ~ niche
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +
  labs(
    x = NULL,
    y = "Agreement",
    title = "Agreement of feature importance between Plain RF and Plain MF RF",
    subtitle = "Each point is one ecotype model."
  ) +
  theme_bw(
    base_size = 9
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 40,
      hjust = 1
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

save_plot(
  method_agreement_figure,
  file.path(
    figure_dir,
    "Figure_FI_5_rf_vs_mf_importance_agreement.png"
  ),
  12.0,
  7.2
)


# 16. Figure: observed zone area and importance concentration ===================

category_values <- sort(
  unique(
    zone_metadata$category2
  )
)

category_colors <- setNames(
  hcl.colors(
    length(category_values),
    palette = "Dark 3"
  ),
  category_values
)

area_concentration_figure <- ggplot(
  zone_importance_concentration,
  aes(
    x = reference_area_km2,
    y = effective_number_of_features,
    color = category2,
    shape = selected_for_detail
  )
) +
  geom_point(
    size = 2.2,
    alpha = 0.78,
    stroke = 0.45
  ) +
  facet_wrap(
    ~ niche,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_x_log10() +
  scale_color_manual(
    values = category_colors
  ) +
  scale_shape_manual(
    values = c(
      `FALSE` = 16,
      `TRUE` = 17
    )
  ) +
  labs(
    x = expression("Observed reference area (km"^2*", log scale)"),
    y = "Effective number of important variables",
    color = "Broad category",
    shape = "Selected for\ndetailed review",
    title = "Observed ecotype area and concentration of feature importance",
    subtitle =
      "The effective number is exp(Shannon entropy) of the normalized positive MDA profile."
  ) +
  theme_bw(
    base_size = 9.5
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      face = "bold"
    )
  )

save_plot(
  area_concentration_figure,
  file.path(
    figure_dir,
    "Figure_FI_6_zone_area_vs_importance_concentration.png"
  ),
  11.2,
  6.0
)


# 17. Category-specific representative and area-extreme figures =================

representative_figure_index <- list()

for (category_value in sort(
  unique(zone_metadata$category2)
)) {
  category_roles <- zone_selection_long[
    category2 == category_value &
      selection_role %in% c(
        "category profile medoid",
        "category largest observed area",
        "category smallest observed area"
      )
  ]
  
  category_zone_selection <- category_roles[
    ,
    .(
      category_role = paste(
        sort(
          unique(selection_role)
        ),
        collapse = "; "
      )
    ),
    by = zoneID
  ]
  
  category_figure_data <- merge(
    zone_top_factors[
      zoneID %in%
        category_zone_selection$zoneID
    ],
    category_zone_selection,
    by = "zoneID",
    all = FALSE,
    sort = FALSE
  )
  
  category_figure_data[
    ,
    role_short := gsub(
      "category ",
      "",
      category_role,
      fixed = TRUE
    )
  ]
  
  category_figure_data[
    ,
    panel_label := paste0(
      "Zone ",
      zoneID,
      " | ",
      role_short,
      " | ",
      niche,
      "\n",
      wrap_text(
        zone_name,
        width = 42L
      )
    )
  ]
  
  setorder(
    category_figure_data,
    panel_label,
    consensus_importance_share,
    feature
  )
  
  category_figure_data[
    ,
    feature_panel := factor(
      paste(
        panel_label,
        feature,
        sep = "|||"
      ),
      levels = unique(
        paste(
          panel_label,
          feature,
          sep = "|||"
        )
      )
    )
  ]
  
  representative_figure <- ggplot(
    category_figure_data,
    aes(
      x = consensus_importance_share,
      y = feature_panel
    )
  ) +
    geom_col(
      width = 0.68,
      fill = "grey35"
    ) +
    facet_wrap(
      ~ panel_label,
      scales = "free_y",
      ncol = 2
    ) +
    scale_y_discrete(
      labels = function(x) {
        sub(
          "^.*\\|\\|\\|",
          "",
          x
        )
      }
    ) +
    scale_x_continuous(
      labels = percent_label,
      expand = expansion(
        mult = c(
          0,
          0.06
        )
      )
    ) +
    labs(
      x = "Consensus normalized permutation importance",
      y = NULL,
      title = paste0(
        "Representative and area-extreme ecotypes: ",
        gsub(
          "_",
          " ",
          category_value,
          fixed = TRUE
        )
      ),
      subtitle =
        "The profile medoid has the highest mean cosine similarity to other zones in its category; area extremes use the observed reference raster."
    ) +
    theme_bw(
      base_size = 8.8
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 7.8
      ),
      panel.spacing = grid::unit(
        0.8,
        "lines"
      )
    )
  
  representative_figure_file <- file.path(
    representative_figure_dir,
    paste0(
      "Figure_FI_category_zones_",
      safe_name(category_value),
      ".png"
    )
  )
  
  save_plot(
    representative_figure,
    representative_figure_file,
    12.5,
    10.0
  )
  
  representative_figure_index[[
    length(representative_figure_index) + 1L
  ]] <- data.table(
    category2 = category_value,
    n_unique_zones = uniqueN(
      category_figure_data$zoneID
    ),
    figure_file = representative_figure_file
  )
}

representative_figure_index <- rbindlist(
  representative_figure_index,
  fill = TRUE
)

fwrite(
  representative_figure_index,
  file.path(
    table_dir,
    "FI_18_category_representative_figure_index.csv"
  )
)


# 18. Analysis settings and completion report ===================================

analysis_settings <- data.table(
  setting = c(
    "primary_importance",
    "normalization",
    "negative_MDA_handling",
    "category_mean_denominator",
    "consensus_definition",
    "representative_profile_niche_weighting",
    "representative_zone_definition",
    "area_definition",
    "common_selection_frequency",
    "common_top5_frequency",
    "minimum_zones_for_category_inference"
  ),
  value = c(
    "unscaled global MeanDecreaseAccuracy stored by randomForest",
    "positive MDA divided by the within-model sum of positive MDA",
    "retained raw; set to zero only when calculating normalized shares",
    "all modeled zones in the category; unselected variables enter as zero",
    "arithmetic mean of Plain RF and Plain MF RF normalized shares within zone",
    "climate and soil each contribute 0.5 to the joint profile",
    "zone with the highest mean pairwise cosine similarity within the category's joint climate-soil profiles",
    "cell-specific km2 summed from terra::cellSize on ecosys_ori.tif",
    as.character(common_selection_frequency),
    as.character(common_top5_frequency),
    as.character(minimum_zones_for_category_inference)
  )
)

fwrite(
  analysis_settings,
  file.path(
    table_dir,
    "FI_00_analysis_settings.csv"
  )
)

cat(
  "\nCOMPLETE\n",
  "Expected models: ",
  nrow(model_specs),
  "\n",
  "Successfully extracted models: ",
  model_audit[load_ok == TRUE, .N],
  "\n",
  "Modeled zones: ",
  length(model_zoneID),
  "\n",
  "Broad categories: ",
  uniqueN(zone_metadata$category2),
  "\n",
  "Primary metric: unscaled global MeanDecreaseAccuracy\n",
  "Tables: ",
  table_dir,
  "\n",
  "Figures: ",
  figure_dir,
  "\n",
  sep = ""
)
