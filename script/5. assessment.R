# Assessment of the selected-variable Multi-Forest workflow
# ==============================================================================
# Run after script 4. This script writes analysis
# tables only; all manuscript figures are created by script 11.

library(data.table)
library(terra)
library(randomForest)

rm(list = ls())
gc()


# 0. Paths and parameters =====================================================

find_project_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "script", "1. data.R"))) return(path)
    parent <- dirname(path)
    if (parent == path) stop("Run inside the repository or set ECOCHINA2_DIR.")
    path <- parent
  }
}

env_root <- Sys.getenv("ECOCHINA2_DIR", unset = "")
project_dir <- if (nzchar(env_root)) {
  normalizePath(env_root, winslash = "/", mustWork = TRUE)
} else {
  find_project_root()
}

climate_test_file <- file.path(project_dir, "results", "test_data.csv")
soil_test_file <- file.path(project_dir, "results", "soil_test_data.csv")
climate_model_dir <- file.path(project_dir, "rf")
soil_model_dir <- file.path(project_dir, "rf_soil")

reference_candidates <- c(
  file.path(project_dir, "raster", "ecosys_ori.tif"),
  file.path(project_dir, "data", "ecosys_ori.tif")
)
reference_file <- reference_candidates[file.exists(reference_candidates)][1]

normal_map_file <- file.path(
  project_dir,
  "result maps",
  "mf_var",
  "assigned_zone_normal_threshold0.4_tol1e-04_novel99_maskNA8_noNovelNormal.tif"
)
assessment_dir <- file.path(project_dir, "assessment")
legacy_assessment_dir <- file.path(project_dir, "assessment_var")
dir.create(assessment_dir, recursive = TRUE, showWarnings = FALSE)

palette_candidates <- c(
  file.path(project_dir, "color_palette_China.csv"),
  file.path(project_dir, "data", "zone_palette.csv")
)
palette_file <- palette_candidates[file.exists(palette_candidates)][1]

zoneID <- c(1:7, 9:50, 52:55)
prob_threshold <- 0.50
base_seed <- 49L
method_name <- "mf_var"
force_assessment <- tolower(trimws(
  Sys.getenv("ECOCHINA2_FORCE_ASSESSMENT", "false")
)) %in% c("1", "true", "yes", "y")

canonical_names <- c(
  "rf_test_zone_metrics.csv",
  "rf_test_summary.csv",
  "climate_test_zone_metrics.csv",
  "soil_test_zone_metrics.csv",
  "climate_test_model_summary.csv",
  "soil_test_model_summary.csv",
  "normal_map_confusion_long.csv",
  "normal_map_confusion_matrix.csv",
  "normal_map_zone_metrics.csv",
  "normal_map_category_confusion_long.csv",
  "normal_map_category_confusion_matrix.csv",
  "normal_map_category_metrics.csv",
  "normal_map_overall_metrics.csv",
  "normal_map_errors_from_original_zone.csv",
  "normal_map_errors_into_assigned_zone.csv"
)
canonical_files <- file.path(assessment_dir, canonical_names)
assessment_sources <- c(
  climate_test_file,
  soil_test_file,
  reference_file,
  normal_map_file,
  palette_file,
  file.path(climate_model_dir, paste0("clm_mfVar_zone", zoneID, ".Rdata")),
  file.path(soil_model_dir, paste0("soil_mf_zone", zoneID, ".Rdata"))
)
assessment_sources <- assessment_sources[
  !is.na(assessment_sources) & file.exists(assessment_sources)
]

files_are_current <- function(outputs, inputs) {
  if (!all(file.exists(outputs))) return(FALSE)
  if (!length(inputs)) return(TRUE)
  min(file.info(outputs)$mtime) >= max(file.info(inputs)$mtime)
}

map_sources <- c(reference_file, normal_map_file)
map_sources <- map_sources[!is.na(map_sources) & file.exists(map_sources)]
map_cache_pairs <- list(
  c(
    file.path(assessment_dir, "normal_map_confusion_long.csv"),
    file.path(assessment_dir, "normal_map_overall_metrics.csv")
  ),
  c(
    file.path(legacy_assessment_dir, "normal_map_confusion_long_var.csv"),
    file.path(legacy_assessment_dir, "normal_map_overall_metrics_var.csv")
  )
)

valid_map_count_cache <- function(files) {
  if (force_assessment || !files_are_current(files, map_sources)) return(FALSE)
  tryCatch({
    confusion_names <- names(fread(files[1], nrows = 0))
    overall_names <- names(fread(files[2], nrows = 0))
    all(c("original_zone") %in% confusion_names) &&
      any(c("assigned_zone", "predicted_zone") %in% confusion_names) &&
      any(c("pixels", "n") %in% confusion_names) &&
      all(c("valid_original_pixels", "compared_pixels") %in% overall_names)
  }, error = function(e) FALSE)
}

map_cache_pair <- NULL
for (candidate in map_cache_pairs) {
  if (valid_map_count_cache(candidate)) {
    map_cache_pair <- candidate
    break
  }
}

# Assessment tables are small but their map confusion table requires a full
# raster pass. Reuse current short outputs before loading models or rasters.
assessment_ready <- !force_assessment &&
  files_are_current(canonical_files, assessment_sources)

if (assessment_ready) {
  cat(
    "[USE EXISTING] Assessment tables are complete.\n",
    "Tables: ", assessment_dir, "\n",
    sep = ""
  )
} else {

required_files <- c(
  climate_test_file,
  soil_test_file,
  reference_file,
  normal_map_file
)
required_files <- required_files[!is.na(required_files)]
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files)) {
  stop("Missing required file(s):\n", paste(missing_files, collapse = "\n"))
}

if (is.na(reference_file)) {
  stop("Missing raster/ecosys_ori.tif.")
}

if (is.na(palette_file)) {
  stop("Missing zone palette. Run script/color_palette.R first.")
}


# 1. Metric helpers ===========================================================

safe_divide <- function(numerator, denominator) {
  ifelse(is.finite(denominator) & denominator > 0,
         numerator / denominator, NA_real_)
}

mean_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

sd_na <- function(x) {
  if (sum(!is.na(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE)
}

auc_rank <- function(observed, probability) {
  n_presence <- sum(observed == 1)
  n_absence <- sum(observed == 0)

  if (!n_presence || !n_absence) return(NA_real_)

  rank_sum <- sum(
    rank(probability, ties.method = "average")[observed == 1]
  )

  (rank_sum - n_presence * (n_presence + 1) / 2) /
    (n_presence * n_absence)
}

metrics_from_counts <- function(TP, TN, FP, FN, auc = NA_real_) {
  sensitivity <- safe_divide(TP, TP + FN)
  specificity <- safe_divide(TN, TN + FP)
  precision <- safe_divide(TP, TP + FP)

  data.table(
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    accuracy = safe_divide(TP + TN, TP + TN + FP + FN),
    balanced_accuracy = (sensitivity + specificity) / 2,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1 = safe_divide(2 * TP, 2 * TP + FP + FN),
    tss = sensitivity + specificity - 1,
    auc = auc
  )
}

binary_metrics <- function(observed, predicted, probability = NULL) {
  metrics_from_counts(
    TP = sum(observed == 1 & predicted == 1),
    TN = sum(observed == 0 & predicted == 0),
    FP = sum(observed == 0 & predicted == 1),
    FN = sum(observed == 1 & predicted == 0),
    auc = if (is.null(probability)) NA_real_ else
      auc_rank(observed, probability)
  )
}

get_varlist <- function(model, model_file) {
  vars <- model$varlist

  if (is.null(vars) || !length(vars)) {
    vars <- rownames(model$importance)
  }

  if (is.null(vars) || !length(vars)) {
    stop("No predictor names found in: ", model_file)
  }

  if (length(vars) != 9L) {
    stop("Expected 9 selected predictors in ", model_file,
         "; found ", length(vars), ".")
  }

  as.character(vars)
}

load_model <- function(model_file, object_name) {
  environment <- new.env(parent = emptyenv())
  load(model_file, envir = environment)

  if (!exists(object_name, envir = environment, inherits = FALSE)) {
    stop("Object '", object_name, "' not found in: ", model_file)
  }

  get(object_name, envir = environment, inherits = FALSE)
}

balanced_rows <- function(data, zone, seed) {
  presence <- which(data$zoneID == zone)
  absence <- which(data$zoneID %in% zoneID & data$zoneID != zone)

  if (!length(presence) || !length(absence)) return(integer())

  set.seed(seed)
  n_each <- min(length(presence), length(absence))

  if (length(presence) > n_each) {
    presence <- sample(presence, n_each)
  }

  if (length(absence) > n_each) {
    absence <- sample(absence, n_each)
  }

  c(presence, absence)
}


# 2. Balanced held-out binary assessment =====================================

climate_test <- as.data.frame(fread(climate_test_file))
soil_test <- as.data.frame(fread(soil_test_file))
climate_test$zoneID <- as.integer(as.character(climate_test$zoneID))
soil_test$zoneID <- as.integer(as.character(soil_test$zoneID))

model_table <- data.table(
  niche = c("climate", "soil"),
  model_dir = c(climate_model_dir, soil_model_dir),
  model_prefix = c("clm_mfVar_zone", "soil_mf_zone"),
  object_name = c("clm_mfVar", "soil_mf")
)

zone_results <- list()

for (model_row in seq_len(nrow(model_table))) {
  niche <- model_table$niche[model_row]
  test_data <- if (niche == "climate") climate_test else soil_test

  for (zone in zoneID) {
    model_file <- file.path(
      model_table$model_dir[model_row],
      paste0(model_table$model_prefix[model_row], zone, ".Rdata")
    )

    if (!file.exists(model_file)) {
      stop("Missing model: ", model_file)
    }

    model <- load_model(model_file, model_table$object_name[model_row])
    varlist <- get_varlist(model, model_file)
    missing_vars <- setdiff(varlist, names(test_data))

    if (length(missing_vars)) {
      stop(
        "Held-out ", niche, " data lack predictor(s) for Zone ", zone, ": ",
        paste(missing_vars, collapse = ", ")
      )
    }

    complete <- complete.cases(test_data[, varlist, drop = FALSE]) &
      test_data$zoneID %in% zoneID
    test_complete <- test_data[complete, , drop = FALSE]
    seed_offset <- if (niche == "climate") zone else 1000L + zone
    rows <- balanced_rows(
      test_complete,
      zone,
      base_seed + seed_offset
    )

    if (!length(rows)) {
      stop("No balanced held-out sample for ", niche, " Zone ", zone, ".")
    }

    x <- test_complete[rows, varlist, drop = FALSE]
    observed <- as.integer(test_complete$zoneID[rows] == zone)
    probability_matrix <- predict(model, x, type = "prob")

    if (!("1" %in% colnames(probability_matrix))) {
      stop("Model lacks class '1': ", model_file)
    }

    probability <- as.numeric(probability_matrix[, "1"])
    finite_probability <- is.finite(probability)
    probability <- probability[finite_probability]
    observed <- observed[finite_probability]

    if (length(unique(observed)) != 2L) {
      stop("Held-out probabilities do not retain both classes for ",
           niche, " Zone ", zone, ".")
    }

    predicted <- as.integer(probability >= prob_threshold)
    metrics <- binary_metrics(observed, predicted, probability)

    zone_results[[length(zone_results) + 1L]] <- cbind(
      data.table(
        method = method_name,
        niche = niche,
        zone = zone,
        threshold = prob_threshold,
        sampling = "equal presence and absence",
        n_test = length(observed),
        presence = sum(observed == 1),
        absence = sum(observed == 0),
        n_predictors = length(varlist)
      ),
      metrics
    )

    rm(model, x, probability_matrix)
  }
}

rf_test <- rbindlist(zone_results, use.names = TRUE)
setorder(rf_test, niche, zone)

rf_summary <- rf_test[, .(
  zones_assessed = .N,
  mean_accuracy = mean_na(accuracy),
  sd_accuracy = sd_na(accuracy),
  mean_balanced_accuracy = mean_na(balanced_accuracy),
  sd_balanced_accuracy = sd_na(balanced_accuracy),
  mean_sensitivity = mean_na(sensitivity),
  sd_sensitivity = sd_na(sensitivity),
  mean_specificity = mean_na(specificity),
  sd_specificity = sd_na(specificity),
  mean_precision = mean_na(precision),
  sd_precision = sd_na(precision),
  mean_f1 = mean_na(f1),
  sd_f1 = sd_na(f1),
  mean_tss = mean_na(tss),
  sd_tss = sd_na(tss),
  mean_auc = mean_na(auc),
  sd_auc = sd_na(auc)
), by = .(method, niche)]

fwrite(rf_test, file.path(assessment_dir, "rf_test_zone_metrics.csv"))
fwrite(rf_summary, file.path(assessment_dir, "rf_test_summary.csv"))
fwrite(
  rf_test[niche == "climate"],
  file.path(assessment_dir, "climate_test_zone_metrics.csv")
)
fwrite(
  rf_test[niche == "soil"],
  file.path(assessment_dir, "soil_test_zone_metrics.csv")
)
fwrite(
  rf_summary[niche == "climate"],
  file.path(assessment_dir, "climate_test_model_summary.csv")
)
fwrite(
  rf_summary[niche == "soil"],
  file.path(assessment_dir, "soil_test_model_summary.csv")
)


# 3. Normal assigned-map assessment ==========================================

if (!is.null(map_cache_pair)) {
  confusion <- fread(map_cache_pair[1])
  if ("method" %in% names(confusion)) confusion <- confusion[method == method_name]
  if ("predicted_zone" %in% names(confusion) &&
      !"assigned_zone" %in% names(confusion)) {
    setnames(confusion, "predicted_zone", "assigned_zone")
  }
  if ("n" %in% names(confusion) && !"pixels" %in% names(confusion)) {
    setnames(confusion, "n", "pixels")
  }
  confusion <- confusion[, .(
    pixels = sum(as.numeric(pixels))
  ), by = .(
    original_zone = as.integer(original_zone),
    assigned_zone = as.integer(assigned_zone)
  )]
  confusion[, method := method_name]
  setcolorder(confusion, c("method", "original_zone", "assigned_zone", "pixels"))

  cached_overall <- fread(map_cache_pair[2])
  if ("method" %in% names(cached_overall)) {
    cached_overall <- cached_overall[method == method_name]
  }
  if (!nrow(cached_overall)) stop("Map cache lacks mf_var overall metrics.")
  valid_original_pixels <- as.numeric(cached_overall$valid_original_pixels[[1]])
  compared_pixels <- as.numeric(cached_overall$compared_pixels[[1]])
  cat("[REUSE MAP COUNTS]", map_cache_pair[1], "\n")
} else {
  reference_map <- rast(reference_file)
  assigned_map <- rast(normal_map_file)

  if (!compareGeom(reference_map, assigned_map, stopOnError = FALSE)) {
    stop("The normal assigned map does not match the reference-map geometry.")
  }

  # Only modeled original zones and valid modeled predictions enter agreement.
  # Thus missing climate/soil cells and all unmodeled values (8, 51 and 56) do
  # not inflate accuracy.
  original <- subst(reference_map, from = zoneID, to = zoneID, others = NA)
  assigned <- subst(assigned_map, from = zoneID, to = zoneID, others = NA)
  names(original) <- "original_zone"
  names(assigned) <- "assigned_zone"

  valid_original_pixels <- as.numeric(
    global(!is.na(original), "sum", na.rm = TRUE)[1, 1]
  )
  compared_pixels <- as.numeric(
    global(!is.na(original) & !is.na(assigned), "sum", na.rm = TRUE)[1, 1]
  )

  confusion <- as.data.table(
    crosstab(c(original, assigned), long = TRUE, useNA = FALSE)
  )
  if (ncol(confusion) != 3L) {
    stop("Unexpected terra::crosstab output for the normal map.")
  }
  setnames(confusion, names(confusion), c("original_zone", "assigned_zone", "pixels"))
  confusion[, `:=`(
    original_zone = as.integer(original_zone),
    assigned_zone = as.integer(assigned_zone),
    pixels = as.numeric(pixels),
    method = method_name
  )]
  setcolorder(confusion, c("method", "original_zone", "assigned_zone", "pixels"))
}

if (sum(confusion$pixels) != compared_pixels) {
  stop("Crosstab total does not equal the common-mask pixel count.")
}

total <- sum(confusion$pixels)
exact_pixels <- confusion[
  original_zone == assigned_zone,
  sum(pixels)
]

map_zone_metrics <- rbindlist(lapply(zoneID, function(zone) {
  TP <- confusion[
    original_zone == zone & assigned_zone == zone,
    sum(pixels)
  ]
  FN <- confusion[
    original_zone == zone & assigned_zone != zone,
    sum(pixels)
  ]
  FP <- confusion[
    original_zone != zone & assigned_zone == zone,
    sum(pixels)
  ]
  TN <- total - TP - FN - FP

  cbind(
    data.table(
      method = method_name,
      zone = zone,
      original_pixels = TP + FN,
      assigned_pixels = TP + FP
    ),
    metrics_from_counts(TP, TN, FP, FN)
  )
}))
map_zone_metrics[, auc := NULL]

# Category agreement uses the same common pixel mask as exact-zone agreement.
palette <- fread(palette_file)

if (!all(c("zoneID", "category2") %in% names(palette))) {
  stop("The zone palette must contain zoneID and category2 columns.")
}

zone_category <- setNames(
  as.character(palette$category2),
  as.character(as.integer(palette$zoneID))
)

confusion[, `:=`(
  original_category = unname(zone_category[as.character(original_zone)]),
  assigned_category = unname(zone_category[as.character(assigned_zone)])
)]

if (confusion[, anyNA(original_category) || anyNA(assigned_category)]) {
  stop("At least one modeled zone lacks a category2 palette entry.")
}

category_confusion <- confusion[, .(
  pixels = sum(pixels)
), by = .(method, original_category, assigned_category)]
setorder(category_confusion, original_category, assigned_category)

category_pixels <- category_confusion[
  original_category == assigned_category,
  sum(pixels)
]

categories <- sort(unique(c(
  category_confusion$original_category,
  category_confusion$assigned_category
)))

map_category_metrics <- rbindlist(lapply(categories, function(category) {
  TP <- category_confusion[
    original_category == category & assigned_category == category,
    sum(pixels)
  ]
  FN <- category_confusion[
    original_category == category & assigned_category != category,
    sum(pixels)
  ]
  FP <- category_confusion[
    original_category != category & assigned_category == category,
    sum(pixels)
  ]
  TN <- total - TP - FN - FP

  cbind(
    data.table(
      method = method_name,
      category = category,
      original_pixels = TP + FN,
      assigned_pixels = TP + FP
    ),
    metrics_from_counts(TP, TN, FP, FN)
  )
}))
map_category_metrics[, auc := NULL]

map_overall <- data.table(
  method = method_name,
  valid_original_pixels = valid_original_pixels,
  compared_pixels = compared_pixels,
  missing_predictions = valid_original_pixels - compared_pixels,
  coverage = safe_divide(compared_pixels, valid_original_pixels),
  accuracy = safe_divide(exact_pixels, compared_pixels),
  exact_zone_agreement = safe_divide(exact_pixels, compared_pixels),
  category_agreement = safe_divide(category_pixels, compared_pixels)
)


# 4. Confusion and error-direction tables ====================================

zone_matrix <- dcast(
  confusion,
  original_zone ~ assigned_zone,
  value.var = "pixels",
  fill = 0
)
category_matrix <- dcast(
  category_confusion,
  original_category ~ assigned_category,
  value.var = "pixels",
  fill = 0
)

errors_from <- confusion[
  original_zone != assigned_zone,
  .(pixels = sum(pixels)),
  by = .(method, original_zone, assigned_zone)
]
original_totals <- confusion[, .(
  original_pixels = sum(pixels)
), by = .(method, original_zone)]
errors_from <- merge(
  errors_from,
  original_totals,
  by = c("method", "original_zone")
)
errors_from[, percent_of_original := 100 * pixels / original_pixels]
setorder(errors_from, original_zone, -pixels)

errors_into <- confusion[
  original_zone != assigned_zone,
  .(pixels = sum(pixels)),
  by = .(method, assigned_zone, original_zone)
]
assigned_totals <- confusion[, .(
  assigned_pixels = sum(pixels)
), by = .(method, assigned_zone)]
errors_into <- merge(
  errors_into,
  assigned_totals,
  by = c("method", "assigned_zone")
)
errors_into[, percent_of_assigned := 100 * pixels / assigned_pixels]
setorder(errors_into, assigned_zone, -pixels)

fwrite(confusion, file.path(assessment_dir, "normal_map_confusion_long.csv"))
fwrite(
  zone_matrix,
  file.path(assessment_dir, "normal_map_confusion_matrix.csv")
)
fwrite(map_zone_metrics, file.path(assessment_dir, "normal_map_zone_metrics.csv"))
fwrite(
  category_confusion,
  file.path(assessment_dir, "normal_map_category_confusion_long.csv")
)
fwrite(
  category_matrix,
  file.path(assessment_dir, "normal_map_category_confusion_matrix.csv")
)
fwrite(
  map_category_metrics,
  file.path(assessment_dir, "normal_map_category_metrics.csv")
)
fwrite(map_overall, file.path(assessment_dir, "normal_map_overall_metrics.csv"))
fwrite(
  errors_from,
  file.path(assessment_dir, "normal_map_errors_from_original_zone.csv")
)
fwrite(
  errors_into,
  file.path(assessment_dir, "normal_map_errors_into_assigned_zone.csv")
)

cat("\nBINARY MODEL SUMMARY\n")
print(rf_summary)
cat("\nNORMAL MAP SUMMARY\n")
print(map_overall)
cat("\nCOMPLETE\nTables: ", assessment_dir, "\n", sep = "")
}
