# 7. Multiclass RF: assessment and projection
# =============================================================================
# Run from the repository root after scripts 2.15, 2.4, and 4.
#
# This robustness workflow:
#   1. joins reference-period climate and static topsoil predictors by cell;
#   2. performs multiclass-compatible backward selection to 30 predictors;
#   3. fits one class-balanced multiclass random forest;
#   4. assesses a held-out, stratified test sample;
#   5. writes class probabilities and assigned-zone maps for the reference
#      period and all six future climate combinations; and
#   6. compares its assigned maps with the binary Multi-Forest maps on their
#      common non-missing pixels.
#
# Multiclass probabilities are not binary dual-suitability values. Therefore
# no 0.4 novel-ecosystem threshold, and no Zone 99, is used here.
# =============================================================================

library(data.table)
library(randomForest)
library(terra)

rm(list = ls())
gc()


# 0. Settings ==================================================================

find_project_root <- function(path = getwd()) {
  configured <- Sys.getenv("ECOCHINA2_DIR", unset = "")
  if (nzchar(configured)) {
    return(normalizePath(configured, winslash = "/", mustWork = TRUE))
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "script", "7. multiclassrf.R"))) return(path)
    parent <- dirname(path)
    if (parent == path) stop("Run inside the repository or set ECOCHINA2_DIR.")
    path <- parent
  }
}

project_root <- find_project_root()

processed_dir <- file.path(project_root, "data", "processed")
binary_climate_dir <- file.path(project_root, "rf")
binary_soil_dir <- file.path(project_root, "rf_soil")
binary_map_dir <- file.path(project_root, "result maps", "mf_var")

# These are the established cache locations from the original workflow. Keep
# them stable: a clean code branch must not force a second multiclass fit or a
# second set of national rasters merely because directories were renamed.
model_dir <- file.path(project_root, "rf_multiclass")
map_dir <- file.path(project_root, "result maps", "multiclass_rf")
table_dir <- file.path(project_root, "assessment", "multiclass_rf")
probability_dir <- file.path(map_dir, "probability")
temp_dir <- file.path(project_root, "tmp_multiclass_rf")

for (d in c(model_dir, map_dir, table_dir, probability_dir, temp_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ECOCHINA2_CLIMATE_DIR and ECOCHINA2_SOIL_RASTER_DIR keep machine-specific
# raw-data locations out of the public code; defaults follow the README layout.
climate_root <- Sys.getenv(
  "ECOCHINA2_CLIMATE_DIR",
  unset = file.path(project_root, "data", "rasters", "climate")
)
soil_raster_dir <- Sys.getenv(
  "ECOCHINA2_SOIL_RASTER_DIR",
  unset = file.path(project_root, "data", "rasters", "soil")
)

zones <- c(1:7, 9:50, 52:55)
base_seed <- 49L
target_predictors <- 30L
selection_ntree <- 100L
selection_repeats <- 10L
final_ntree <- 500L
candidate_per_zone <- 10000L
selection_per_zone <- 1000L
max_train_per_zone <- 5000L
max_test_per_zone <- 2000L
train_fraction <- 0.70
detected_cores <- parallel::detectCores()
if (is.na(detected_cores)) detected_cores <- 1L
predict_cores <- min(4L, max(1L, detected_cores - 1L))
force_rebuild <- identical(
  tolower(Sys.getenv("ECOCHINA2_FORCE_MULTICLASS", "false")),
  "true"
)

terraOptions(tempdir = temp_dir, memfrac = 0.25)

scenario_table <- data.table(
  scenario = c(
    "normal",
    "2011-2040SSP245", "2041-2070SSP245", "2071-2100SSP245",
    "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585"
  ),
  climate_subdir = c(
    "Normal_1961_1990",
    "8GCMs_ensemble_ssp245_2011-2040",
    "8GCMs_ensemble_ssp245_2041-2070",
    "8GCMs_ensemble_ssp245_2071-2100",
    "8GCMs_ensemble_ssp585_2011-2040",
    "8GCMs_ensemble_ssp585_2041-2070",
    "8GCMs_ensemble_ssp585_2071-2100"
  )
)

model_file <- file.path(model_dir, "multiclass_climate_soil_rf.Rdata")
selection_file <- file.path(table_dir, "multiclass_mcRFop_variable_selection.csv")
selected_file <- file.path(table_dir, "multiclass_selected_variables.csv")
train_cache <- file.path(temp_dir, "multiclass_train_data.rds")
test_cache <- file.path(temp_dir, "multiclass_test_data.rds")


# 1. Small helpers =============================================================

first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths) | dir.exists(paths)]
  if (!length(hit)) {
    stop(label, " not found. Checked:\n", paste(paths, collapse = "\n"))
  }
  hit[[1]]
}

read_rf <- function(file) {
  if (tolower(tools::file_ext(file)) == "rds") {
    objects <- list(readRDS(file))
  } else {
    environment <- new.env(parent = emptyenv())
    object_names <- load(file, envir = environment)
    objects <- mget(object_names, envir = environment, inherits = FALSE)
  }
  candidates <- Filter(function(x) inherits(x, "randomForest"), objects)
  if (length(candidates) == 1L) return(candidates[[1]])
  candidates <- Filter(
    function(x) inherits(x, "randomForest"),
    unlist(lapply(objects, function(x) if (is.list(x)) x else list()), recursive = FALSE)
  )
  if (length(candidates) != 1L) {
    stop("Could not identify one randomForest object in: ", file)
  }
  candidates[[1]]
}

model_variables <- function(model) {
  variables <- model$varlist
  if (is.null(variables) || !length(variables)) {
    variables <- rownames(randomForest::importance(model, scale = FALSE))
  }
  if (is.null(variables) || !length(variables)) {
    stop("A fitted RF does not contain predictor names.")
  }
  as.character(variables)
}

read_table_variables <- function(files) {
  file <- files[file.exists(files)][1]
  if (is.na(file)) return(character())

  variables <- names(fread(file, nrows = 0))
  metadata <- c(
    "cell", "zoneid", "zone", "x", "y", "lon", "lat", "longitude",
    "latitude", "dem", "elevation", "china_90m", "split", "split_id"
  )
  keep <- !(tolower(variables) %in% metadata) &
    !grepl("^zone[0-9]+$", tolower(variables)) &
    !grepl("^v[0-9]+$", tolower(variables))
  variables <- variables[keep]

  # A leading row-name column is a common legacy CSV artifact.
  setdiff(unique(variables), c("row", "X", "V1"))
}

union_binary_variables <- function(niche) {
  files <- if (niche == "climate") {
    file.path(binary_climate_dir, paste0("clm_mfVar_zone", zones, ".Rdata"))
  } else {
    file.path(binary_soil_dir, paste0("soil_mf_zone", zones, ".Rdata"))
  }
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop("Missing final binary models:\n", paste(missing, collapse = "\n"))
  }
  unique(unlist(lapply(files, function(f) model_variables(read_rf(f)))))
}

read_stack <- function(variables, raster_dir, layer_names = variables) {
  files <- file.path(raster_dir, paste0(variables, ".tif"))
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop("Missing predictor rasters:\n", paste(missing, collapse = "\n"))
  }
  x <- rast(files)
  names(x) <- layer_names
  x
}

extract_values <- function(x, points) {
  p <- if (same.crs(points, x)) points else project(points, crs(x))
  answer <- as.data.table(terra::extract(x, p))
  if ("ID" %in% names(answer)) answer[, ID := NULL]
  answer
}

extract_predictors <- function(cells, reference, climate, soil) {
  xy <- xyFromCell(reference, cells$cell)
  points <- vect(
    data.frame(x = xy[, 1], y = xy[, 2]),
    geom = c("x", "y"),
    crs = crs(reference)
  )
  cbind(
    cells,
    extract_values(climate, points),
    extract_values(soil, points)
  )
}

safe_divide <- function(a, b) {
  ifelse(is.finite(b) & b > 0, a / b, NA_real_)
}

rank_auc <- function(y, probability) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(probability, ties.method = "average")[y == 1]) -
     n1 * (n1 + 1) / 2) / (n1 * n0)
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

files_current <- function(outputs, inputs) {
  outputs <- outputs[!is.na(outputs) & nzchar(outputs)]
  if (!length(outputs) || !all(file.exists(outputs))) return(FALSE)
  inputs <- inputs[!is.na(inputs) & nzchar(inputs) & file.exists(inputs)]
  if (!length(inputs)) return(TRUE)
  isTRUE(min(file.info(outputs)$mtime) >= max(file.info(inputs)$mtime))
}


# 2. Candidate predictors and reference cells =================================

reference_file <- first_existing(
  c(
    file.path(project_root, "raster", "ecosys_ori.tif"),
    file.path(project_root, "data", "ecosys_ori.tif"),
    file.path(processed_dir, "ecosys_ori.tif")
  ),
  "Reference vegetation raster"
)

cell_file <- first_existing(
  c(
    file.path(project_root, "data raw", "1. zoneID_Clm_800m_Normal_1961_1990SY.csv"),
    file.path(project_root, "results", "train_data.csv"),
    file.path(project_root, "data", "zone_cells.csv"),
    file.path(project_root, "data", "climate_reference.csv"),
    file.path(processed_dir, "reference_cells.csv"),
    file.path(processed_dir, "zone_cells.csv"),
    file.path(processed_dir, "climate_reference.csv"),
    file.path(processed_dir, "zone_climate_reference.csv"),
    file.path(project_root, "data", "reference_cells.csv"),
    file.path(project_root, "data", "climate_train.csv")
  ),
  "Reference cell table with cell and zoneID columns"
)

climate_variables <- read_table_variables(
  c(
    file.path(project_root, "results", "train_data.csv"),
    file.path(processed_dir, "climate_train.csv"),
    file.path(processed_dir, "train_data.csv"),
    file.path(project_root, "data", "climate_train.csv")
  )
)
soil_variables <- read_table_variables(
  c(
    file.path(project_root, "results", "soil_train_data.csv"),
    file.path(processed_dir, "soil_train.csv"),
    file.path(processed_dir, "soil_train_data.csv"),
    file.path(project_root, "data", "soil_train.csv")
  )
)

# The processed training tables are preferred because they retain the complete
# pre-selection candidate sets. Final binary-model unions are a reproducible
# fallback when only the fitted models and rasters are available.
if (!length(climate_variables)) {
  climate_variables <- union_binary_variables("climate")
}
if (!length(soil_variables)) {
  soil_variables <- union_binary_variables("soil")
}

soil_model_variables <- paste0("soil_", soil_variables)
all_variables <- c(climate_variables, soil_model_variables)
if (anyDuplicated(all_variables)) stop("Duplicated climate/soil predictor names.")
if (length(all_variables) < target_predictors) {
  stop("Only ", length(all_variables), " candidate predictors; 30 are required.")
}

normal_climate_dir <- file.path(climate_root, "Normal_1961_1990")
if (!dir.exists(normal_climate_dir) &&
    basename(normalizePath(climate_root, winslash = "/", mustWork = FALSE)) ==
      "Normal_1961_1990") {
  normal_climate_dir <- climate_root
}

reference <- rast(reference_file)


# 3. Stratified training and held-out test data ================================

cache_ok <- file.exists(train_cache) && file.exists(test_cache)
if (cache_ok) {
  train_data <- as.data.table(readRDS(train_cache))
  test_data <- as.data.table(readRDS(test_cache))
  cache_ok <- all(c("cell", "zoneID", all_variables) %in% names(train_data)) &&
    all(c("cell", "zoneID", all_variables) %in% names(test_data))
}

if (!cache_ok) {
  climate_normal <- read_stack(climate_variables, normal_climate_dir)
  soil_normal <- read_stack(
    soil_variables,
    soil_raster_dir,
    layer_names = soil_model_variables
  )
  cells <- fread(cell_file, select = c("cell", "zoneID"))
  cells[, `:=`(
    cell = as.integer(cell),
    zoneID = as.integer(as.character(zoneID))
  )]
  cells <- cells[zoneID %in% zones & !is.na(cell)]

  missing_zones <- setdiff(zones, unique(cells$zoneID))
  if (length(missing_zones)) {
    stop("Zones absent from the reference cell table: ",
         paste(missing_zones, collapse = ", "))
  }

  set.seed(base_seed)
  sampled_cells <- cells[
    , .SD[sample.int(.N, min(.N, candidate_per_zone))],
    by = zoneID
  ]
  rm(cells)
  gc()

  sample_data <- extract_predictors(
    sampled_cells[, .(cell, zoneID)],
    reference,
    climate_normal,
    soil_normal
  )
  sample_data <- sample_data[complete.cases(sample_data[, ..all_variables])]

  missing_zones <- setdiff(zones, unique(sample_data$zoneID))
  if (length(missing_zones)) {
    stop("Zones without complete climate and soil data: ",
         paste(missing_zones, collapse = ", "))
  }

  train_list <- vector("list", length(zones))
  test_list <- vector("list", length(zones))
  for (j in seq_along(zones)) {
    z <- zones[j]
    d <- sample_data[zoneID == z]
    set.seed(base_seed + z)
    order_index <- sample.int(nrow(d))
    n_train <- min(max_train_per_zone, floor(train_fraction * nrow(d)))
    n_test <- min(max_test_per_zone, nrow(d) - n_train)
    if (n_train < 2L || n_test < 1L) {
      stop("Insufficient complete observations for Zone ", z, ".")
    }
    train_list[[j]] <- d[order_index[seq_len(n_train)]]
    test_list[[j]] <- d[order_index[n_train + seq_len(n_test)]]
  }
  train_data <- rbindlist(train_list)
  test_data <- rbindlist(test_list)
  saveRDS(train_data, train_cache)
  saveRDS(test_data, test_cache)
}

train_data[, `:=`(
  cell = as.integer(cell),
  zoneID = as.integer(as.character(zoneID))
)]
test_data[, `:=`(
  cell = as.integer(cell),
  zoneID = as.integer(as.character(zoneID))
)]
missing_train_zones <- setdiff(zones, unique(train_data$zoneID))
missing_test_zones <- setdiff(zones, unique(test_data$zoneID))
if (length(missing_train_zones) || length(missing_test_zones)) {
  stop(
    "Cached/generated split is missing modeled zones. Train: ",
    paste(missing_train_zones, collapse = ", "),
    "; test: ", paste(missing_test_zones, collapse = ", ")
  )
}
if (any(!complete.cases(train_data[, ..all_variables])) ||
    any(!complete.cases(test_data[, ..all_variables]))) {
  stop("Cached/generated train or test predictors contain missing values.")
}

setorder(train_data, zoneID, cell)
setorder(test_data, zoneID, cell)
fwrite(
  merge(
    train_data[, .(n_train = .N), by = zoneID],
    test_data[, .(n_test = .N), by = zoneID],
    by = "zoneID",
    all = TRUE
  ),
  file.path(table_dir, "multiclass_sample_counts.csv")
)


# 4. Multiclass-compatible backward selection =================================

fit_balanced_rf <- function(x, y, ntree, importance = TRUE) {
  class_n <- table(y)
  sample_n <- min(as.integer(class_n))
  randomForest(
    x = x,
    y = y,
    ntree = ntree,
    mtry = max(1L, floor(sqrt(ncol(x)))),
    importance = importance,
    replace = TRUE,
    sampsize = rep(sample_n, length(class_n))
  )
}

backward_select <- function(x, y, target, ntree, repeats, seed) {
  current <- names(x)
  path <- list()
  step <- 0L

  while (length(current) >= target) {
    step <- step + 1L
    importance_sum <- setNames(numeric(length(current)), current)
    oob_accuracy <- numeric(repeats)

    for (repeat_id in seq_len(repeats)) {
      set.seed(seed + step * 100L + repeat_id)
      fit <- fit_balanced_rf(x[, current, drop = FALSE], y, ntree)
      imp <- randomForest::importance(fit, scale = FALSE)
      importance_sum <- importance_sum + imp[current, "MeanDecreaseAccuracy"]
      oob_accuracy[repeat_id] <- 1 - tail(fit$err.rate[, "OOB"], 1)
    }

    importance_mean <- importance_sum / repeats
    path[[step]] <- data.table(
      step = step,
      n_predictors = length(current),
      mean_oob_accuracy = mean(oob_accuracy),
      variables = paste(current, collapse = ",")
    )

    if (length(current) == target) break
    n_remove <- min(2L, length(current) - target)
    remove <- names(sort(importance_mean, decreasing = FALSE))[seq_len(n_remove)]
    current <- setdiff(current, remove)
  }

  list(variables = current, path = rbindlist(path))
}

cached_multiclass <- NULL
reuse_selection <- FALSE
if (!force_rebuild && file.exists(model_file)) {
  cached_multiclass <- tryCatch(read_rf(model_file), error = function(e) NULL)
  if (!is.null(cached_multiclass)) {
    cached_variables <- model_variables(cached_multiclass)
    reuse_selection <- length(cached_variables) == target_predictors &&
      all(cached_variables %in% all_variables)
    if (reuse_selection) {
      selected_variables <- cached_variables
      cat("[REUSE MULTICLASS MODEL]", model_file, "\n")
    }
  }
}
if (!reuse_selection && file.exists(selected_file)) {
  selected_table <- fread(selected_file)
  reuse_selection <- "variable" %in% names(selected_table) &&
    nrow(selected_table) == target_predictors &&
    all(selected_table$variable %in% all_variables)
}
if (!reuse_selection && file.exists(selection_file)) {
  selection_table <- fread(selection_file)
  variable_column <- intersect(c("variable", "variables"), names(selection_table))[1]
  if (!is.na(variable_column)) {
    paths <- lapply(selection_table[[variable_column]], function(x) {
      if (is.na(x) || !nzchar(x)) character() else trimws(strsplit(x, ",", fixed = TRUE)[[1]])
    })
    target_row <- which(lengths(paths) == target_predictors)[1]
    if (!is.na(target_row) && all(paths[[target_row]] %in% all_variables)) {
      selected_variables <- paths[[target_row]]
      reuse_selection <- TRUE
      cat("[REUSE MULTICLASS SELECTION]", selection_file, "\n")
    }
  }
}

if (reuse_selection && !exists("selected_variables")) {
  selected_variables <- selected_table$variable
} else if (!reuse_selection) {
  set.seed(base_seed + 1000L)
  selection_data <- train_data[
    , .SD[sample.int(.N, min(.N, selection_per_zone))],
    by = zoneID
  ]
  selection_result <- backward_select(
    x = as.data.frame(selection_data[, ..all_variables]),
    y = factor(selection_data$zoneID, levels = zones),
    target = target_predictors,
    ntree = selection_ntree,
    repeats = selection_repeats,
    seed = base_seed + 2000L
  )
  selected_variables <- selection_result$variables
  fwrite(selection_result$path, selection_file)
  fwrite(
    data.table(
      order = seq_along(selected_variables),
      variable = selected_variables,
      niche = ifelse(grepl("^soil_", selected_variables), "soil", "climate"),
      raster_variable = sub("^soil_", "", selected_variables)
    ),
    selected_file
  )
}

if (length(selected_variables) != target_predictors) {
  stop("Backward selection did not retain exactly 30 predictors.")
}
if (!any(selected_variables %in% climate_variables) ||
    !any(selected_variables %in% soil_model_variables)) {
  stop("The 30-predictor multiclass set must retain climate and soil variables.")
}


# 5. Final class-balanced model and held-out assessment ========================

if (!is.null(cached_multiclass) && reuse_selection && !force_rebuild) {
  multiclass_rf <- cached_multiclass
} else {
  x_train <- as.data.frame(train_data[, ..selected_variables])
  y_train <- factor(train_data$zoneID, levels = zones)
  set.seed(base_seed)
  multiclass_rf <- fit_balanced_rf(x_train, y_train, final_ntree)
  multiclass_rf$varlist <- selected_variables
  multiclass_rf$zones <- zones
  multiclass_rf$climate_vars <- intersect(selected_variables, climate_variables)
  multiclass_rf$soil_vars <- intersect(selected_variables, soil_model_variables)
  multiclass_rf$base_seed <- base_seed
  save(multiclass_rf, file = model_file)
  cat("[SAVED MULTICLASS MODEL]", model_file, "\n")
}

x_test <- as.data.frame(test_data[, ..selected_variables])
y_test <- factor(test_data$zoneID, levels = zones)
test_probability <- predict(multiclass_rf, x_test, type = "prob")
test_prediction <- factor(
  colnames(test_probability)[max.col(test_probability, ties.method = "first")],
  levels = as.character(zones)
)

zone_test_metrics <- rbindlist(lapply(zones, function(z) {
  observed <- as.integer(y_test == as.character(z))
  predicted <- as.integer(test_prediction == as.character(z))
  tp <- sum(observed == 1 & predicted == 1)
  tn <- sum(observed == 0 & predicted == 0)
  fp <- sum(observed == 0 & predicted == 1)
  fn <- sum(observed == 1 & predicted == 0)
  sensitivity <- safe_divide(tp, tp + fn)
  specificity <- safe_divide(tn, tn + fp)
  precision <- safe_divide(tp, tp + fp)
  data.table(
    model = "multiclass_rf",
    zoneID = z,
    n_test = tp + fn,
    TP = tp, TN = tn, FP = fp, FN = fn,
    accuracy = safe_divide(tp + tn, tp + tn + fp + fn),
    balanced_accuracy = (sensitivity + specificity) / 2,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1 = safe_divide(2 * tp, 2 * tp + fp + fn),
    tss = sensitivity + specificity - 1,
    auc = rank_auc(observed, test_probability[, as.character(z)])
  )
}))

overall_test_metrics <- data.table(
  model = "multiclass_rf",
  n_test = length(y_test),
  n_predictors = length(selected_variables),
  n_trees = final_ntree,
  oob_accuracy = 1 - tail(multiclass_rf$err.rate[, "OOB"], 1),
  overall_accuracy = mean(y_test == test_prediction),
  macro_balanced_accuracy = mean_or_na(zone_test_metrics$balanced_accuracy),
  macro_sensitivity = mean_or_na(zone_test_metrics$sensitivity),
  macro_specificity = mean_or_na(zone_test_metrics$specificity),
  macro_precision = mean_or_na(zone_test_metrics$precision),
  macro_f1 = mean_or_na(zone_test_metrics$f1),
  macro_tss = mean_or_na(zone_test_metrics$tss),
  macro_auc = mean_or_na(zone_test_metrics$auc)
)

fwrite(overall_test_metrics, file.path(table_dir, "multiclass_test_overall.csv"))
fwrite(zone_test_metrics, file.path(table_dir, "multiclass_test_by_zone.csv"))
fwrite(overall_test_metrics, file.path(table_dir, "multiclass_rf_test_overall.csv"))
fwrite(zone_test_metrics, file.path(table_dir, "multiclass_rf_test_zone_metrics.csv"))
fwrite(
  as.data.table(table(observed = y_test, predicted = test_prediction)),
  file.path(table_dir, "multiclass_test_confusion.csv")
)
fwrite(
  as.data.table(table(observed = y_test, predicted = test_prediction)),
  file.path(table_dir, "multiclass_rf_test_confusion_long.csv")
)
fwrite(
  data.table(
    cell = test_data$cell,
    observed_zone = as.integer(as.character(y_test)),
    predicted_zone = as.integer(as.character(test_prediction))
  ),
  file.path(table_dir, "multiclass_test_predictions.csv")
)
fwrite(
  data.table(
    cell = test_data$cell,
    observed_zone = as.integer(as.character(y_test)),
    predicted_zone = as.integer(as.character(test_prediction))
  ),
  file.path(table_dir, "multiclass_rf_test_predictions.csv")
)


# 6. Reference and future projections =========================================

rf_probability <- function(model, data) {
  probability <- predict(model, as.data.frame(data), type = "prob")
  colnames(probability) <- paste0("zone_", colnames(probability))
  as.data.frame(probability)
}

argmax_first <- function(values) {
  if (is.null(dim(values))) {
    if (all(is.na(values))) return(NA_integer_)
    return(which.max(values))
  }
  answer <- rep(NA_integer_, nrow(values))
  valid <- rowSums(is.finite(values)) > 0
  answer[valid] <- max.col(values[valid, , drop = FALSE], ties.method = "first")
  answer
}

projection_area <- function(map, scenario) {
  area <- cellSize(map, unit = "km")
  summary <- as.data.table(zonal(area, map, fun = "sum", na.rm = TRUE))
  if (ncol(summary) < 2L) return(data.table())
  setnames(summary, names(summary)[1:2], c("zoneID", "area_km2"))
  summary[, `:=`(scenario = scenario, zoneID = as.integer(zoneID))]
  summary[, .(scenario, zoneID, area_km2)]
}

map_candidates <- function(scenario) {
  established <- if (scenario == "normal") {
    file.path(map_dir, "assigned_zone_normal_multiclass_rf.tif")
  } else {
    file.path(map_dir, paste0("assigned_zone_", scenario, "_multiclass_rf.tif"))
  }
  unique(c(established, file.path(map_dir, paste0("assigned_zone_", scenario, ".tif"))))
}

valid_assigned_map <- function(file) {
  if (!file.exists(file)) return(FALSE)
  tryCatch({
    x <- rast(file)
    nlyr(x) == 1L && hasValues(x) &&
      compareGeom(x, reference, stopOnError = FALSE)
  }, error = function(e) FALSE)
}

valid_probability_stack <- function(file) {
  if (!file.exists(file)) return(FALSE)
  tryCatch({
    x <- rast(file)
    nlyr(x) == length(zones) && hasValues(x) && files_current(file, model_file)
  }, error = function(e) FALSE)
}

area_cache_file <- file.path(table_dir, "multiclass_area_change.csv")
cached_area <- if (!force_rebuild && file.exists(area_cache_file)) {
  tryCatch(fread(area_cache_file), error = function(e) data.table())
} else {
  data.table()
}

project_scenario <- function(scenario, climate_subdir) {
  cat("[MULTICLASS PROJECTION] ", scenario, "\n", sep = "")
  candidates <- map_candidates(scenario)
  reusable <- candidates[
    vapply(candidates, valid_assigned_map, logical(1)) &
      vapply(candidates, files_current, logical(1), inputs = model_file)
  ]
  probability_file <- file.path(
    probability_dir,
    paste0("class_probability_", scenario, ".tif")
  )

  if (!force_rebuild && length(reusable)) {
    map_file <- reusable[[1]]
    cat("[REUSE MULTICLASS MAP]", map_file, "\n")
    scenario_value <- scenario
    area <- if (nrow(cached_area) && files_current(area_cache_file, map_file) &&
      all(c("scenario", "zoneID", "area_km2") %in% names(cached_area)) &&
      scenario %in% cached_area$scenario) {
      cached_area[scenario == scenario_value, .(scenario, zoneID, area_km2)]
    } else {
      projection_area(rast(map_file), scenario)
    }
    return(list(
      scenario = scenario,
      map_file = map_file,
      probability_file = if (valid_probability_stack(probability_file)) {
        probability_file
      } else {
        NA_character_
      },
      area = area
    ))
  }

  scenario_climate_dir <- file.path(climate_root, climate_subdir)
  if (scenario == "normal" && identical(normal_climate_dir, climate_root)) {
    scenario_climate_dir <- normal_climate_dir
  }
  if (!dir.exists(scenario_climate_dir)) {
    stop("Climate directory not found: ", scenario_climate_dir)
  }

  selected_climate <- intersect(selected_variables, climate_variables)
  selected_soil_model <- intersect(selected_variables, soil_model_variables)
  selected_soil_raw <- sub("^soil_", "", selected_soil_model)

  climate <- read_stack(selected_climate, scenario_climate_dir)
  soil <- read_stack(
    selected_soil_raw,
    soil_raster_dir,
    layer_names = selected_soil_model
  )
  if (!compareGeom(soil, climate[[1]], stopOnError = FALSE)) {
    soil <- resample(soil, climate[[1]], method = "bilinear")
    names(soil) <- selected_soil_model
  }

  predictors <- c(climate, soil)[[selected_variables]]
  probability <- terra::predict(
    predictors,
    multiclass_rf,
    fun = rf_probability,
    na.rm = TRUE,
    cores = predict_cores,
    cpkgs = "randomForest",
    filename = probability_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = c("COMPRESS=LZW", "BIGTIFF=YES")
    )
  )
  names(probability) <- paste0("zone_", zones)

  index_map <- app(probability, argmax_first)
  climate_grid_map <- subst(
    index_map,
    from = seq_along(zones),
    to = zones,
    others = NA
  )
  names(climate_grid_map) <- "zoneID"

  if (!compareGeom(climate_grid_map, reference, stopOnError = FALSE)) {
    climate_grid_map <- resample(climate_grid_map, reference, method = "near")
  }
  reference_mask <- subst(reference, zones, zones, others = NA)
  map_file <- candidates[[1]]
  assigned <- mask(
    climate_grid_map,
    reference_mask,
    filename = map_file,
    overwrite = TRUE,
    wopt = list(datatype = "INT2S", gdal = "COMPRESS=LZW")
  )

  assigned_values <- as.integer(freq(assigned)$value)
  assigned_values <- assigned_values[!is.na(assigned_values)]
  bad_values <- setdiff(assigned_values, zones)
  if (length(bad_values)) {
    stop("Unexpected multiclass map values: ", paste(bad_values, collapse = ", "))
  }

  list(
    scenario = scenario,
    map_file = map_file,
    probability_file = probability_file,
    area = projection_area(assigned, scenario)
  )
}

projection_results <- lapply(seq_len(nrow(scenario_table)), function(i) {
  project_scenario(scenario_table$scenario[i], scenario_table$climate_subdir[i])
})

projection_inventory <- rbindlist(lapply(projection_results, function(x) {
  data.table(
    scenario = x$scenario,
    assigned_map = x$map_file,
    class_probability = x$probability_file
  )
}))
area_table <- rbindlist(lapply(projection_results, `[[`, "area"), fill = TRUE)
area_table <- merge(
  CJ(scenario = scenario_table$scenario, zoneID = zones, unique = TRUE),
  area_table,
  by = c("scenario", "zoneID"),
  all.x = TRUE
)
area_table[is.na(area_km2), area_km2 := 0]
normal_area <- area_table[scenario == "normal", .(zoneID, normal_area_km2 = area_km2)]
area_table <- merge(area_table, normal_area, by = "zoneID", all.x = TRUE)
area_table[, `:=`(
  change_km2 = area_km2 - normal_area_km2,
  change_percent = 100 * safe_divide(area_km2 - normal_area_km2, normal_area_km2)
)]
area_table[, scenario_order := match(scenario, scenario_table$scenario)]
setorder(area_table, scenario_order, zoneID)
area_table[, scenario_order := NULL]

fwrite(projection_inventory, file.path(table_dir, "multiclass_projection_inventory.csv"))
fwrite(area_table, file.path(table_dir, "multiclass_area_change.csv"))


# 7. Reference-map reconstruction on the modeled-zone common mask =============

assess_reference_map <- function(original, predicted) {
  if (!compareGeom(predicted, original, stopOnError = FALSE)) {
    predicted <- resample(predicted, original, method = "near")
  }
  original <- subst(original, zones, zones, others = NA)
  predicted <- subst(predicted, zones, zones, others = NA)
  names(original) <- "reference_zone"
  names(predicted) <- "predicted_zone"

  valid_reference <- global(!is.na(original), "sum", na.rm = TRUE)[1, 1]
  confusion <- as.data.table(
    crosstab(c(original, predicted), long = TRUE, useNA = FALSE)
  )
  setnames(confusion, names(confusion)[1:3], c(
    "reference_zone", "predicted_zone", "pixels"
  ))
  confusion[, `:=`(
    model = "multiclass_rf",
    reference_zone = as.integer(reference_zone),
    predicted_zone = as.integer(predicted_zone),
    pixels = as.numeric(pixels)
  )]
  n_common <- sum(confusion$pixels)

  by_zone <- rbindlist(lapply(zones, function(z) {
    tp <- confusion[reference_zone == z & predicted_zone == z, sum(pixels)]
    fn <- confusion[reference_zone == z & predicted_zone != z, sum(pixels)]
    fp <- confusion[reference_zone != z & predicted_zone == z, sum(pixels)]
    tn <- n_common - tp - fn - fp
    sensitivity <- safe_divide(tp, tp + fn)
    specificity <- safe_divide(tn, tn + fp)
    precision <- safe_divide(tp, tp + fp)
    data.table(
      model = "multiclass_rf", zoneID = z,
      TP = tp, TN = tn, FP = fp, FN = fn,
      reference_pixels = tp + fn,
      predicted_pixels = tp + fp,
      accuracy = safe_divide(tp + tn, n_common),
      balanced_accuracy = (sensitivity + specificity) / 2,
      sensitivity = sensitivity,
      specificity = specificity,
      precision = precision,
      f1 = safe_divide(2 * tp, 2 * tp + fp + fn),
      tss = sensitivity + specificity - 1
    )
  }))

  overall <- data.table(
    model = "multiclass_rf",
    valid_reference_pixels = valid_reference,
    common_pixels = n_common,
    coverage = safe_divide(n_common, valid_reference),
    exact_agreement = safe_divide(
      confusion[reference_zone == predicted_zone, sum(pixels)], n_common
    ),
    macro_balanced_accuracy = mean_or_na(by_zone$balanced_accuracy),
    macro_sensitivity = mean_or_na(by_zone$sensitivity),
    macro_specificity = mean_or_na(by_zone$specificity),
    macro_precision = mean_or_na(by_zone$precision),
    macro_f1 = mean_or_na(by_zone$f1),
    macro_tss = mean_or_na(by_zone$tss)
  )
  list(overall = overall, by_zone = by_zone, confusion = confusion)
}

normal_map_file <- projection_inventory[
  scenario == "normal", assigned_map
]
reference_files <- file.path(
  table_dir,
  c(
    "reference_map_overall.csv", "reference_map_by_zone.csv",
    "reference_zone_metrics.csv", "reference_map_confusion.csv"
  )
)
if (!force_rebuild && files_current(reference_files, c(reference_file, normal_map_file))) {
  cat("[REUSE MULTICLASS REFERENCE ASSESSMENT]", table_dir, "\n")
} else {
  reference_assessment <- assess_reference_map(reference, rast(normal_map_file))
  fwrite(reference_assessment$overall, reference_files[[1]])
  fwrite(reference_assessment$by_zone, reference_files[[2]])
  fwrite(reference_assessment$by_zone, reference_files[[3]])
  fwrite(reference_assessment$confusion, reference_files[[4]])
}


# 8. Common-mask comparison with the binary Multi-Forest ======================

compare_assigned_maps <- function(binary_file, multiclass_file, scenario) {
  binary <- rast(binary_file)
  multiclass <- rast(multiclass_file)
  if (!compareGeom(multiclass, binary, stopOnError = FALSE)) {
    multiclass <- resample(multiclass, binary, method = "near")
  }
  names(binary) <- "binary_zone"
  names(multiclass) <- "multiclass_zone"

  confusion <- as.data.table(crosstab(c(binary, multiclass), long = TRUE, useNA = FALSE))
  setnames(confusion, names(confusion)[1:3], c("binary_zone", "multiclass_zone", "pixels"))
  confusion[, `:=`(
    scenario = scenario,
    binary_zone = as.integer(binary_zone),
    multiclass_zone = as.integer(multiclass_zone),
    pixels = as.numeric(pixels)
  )]

  n_common <- sum(confusion$pixels)
  n_current <- confusion[binary_zone != 99, sum(pixels)]
  summary <- data.table(
    scenario = scenario,
    common_pixels = n_common,
    current_zone_common_pixels = n_current,
    binary_novel_pixels = confusion[binary_zone == 99, sum(pixels)],
    binary_novel_share = safe_divide(
      confusion[binary_zone == 99, sum(pixels)], n_common
    ),
    exact_agreement_all = safe_divide(
      confusion[binary_zone == multiclass_zone, sum(pixels)], n_common
    ),
    exact_agreement_current_zones = safe_divide(
      confusion[binary_zone != 99 & binary_zone == multiclass_zone, sum(pixels)],
      n_current
    )
  )
  list(summary = summary, confusion = confusion)
}

comparison_summary_file <- file.path(table_dir, "binary_mf_vs_multiclass_common_mask.csv")
comparison_confusion_file <- file.path(table_dir, "binary_mf_vs_multiclass_confusion.csv")
binary_comparison_inputs <- vapply(scenario_table$scenario, function(scenario) {
  candidates <- c(
    file.path(
      binary_map_dir,
      paste0(
        "assigned_zone_", scenario,
        "_threshold0.4_tol1e-04_novel99_maskNA8_noNovelNormal.tif"
      )
    ),
    file.path(project_root, "maps", paste0("assigned_zone_", scenario, ".tif"))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[[1]] else NA_character_
}, character(1))
comparison_inputs <- c(projection_inventory$assigned_map, binary_comparison_inputs)
comparison_cached <- FALSE
if (!force_rebuild && all(file.exists(c(comparison_summary_file, comparison_confusion_file)))) {
  cached_summary <- tryCatch(fread(comparison_summary_file), error = function(e) data.table())
  comparison_cached <- "scenario" %in% names(cached_summary) &&
    all(scenario_table$scenario %in% cached_summary$scenario) &&
    files_current(
      c(comparison_summary_file, comparison_confusion_file),
      comparison_inputs
    )
}

comparison_results <- if (comparison_cached) list() else lapply(seq_len(nrow(projection_inventory)), function(i) {
  scenario <- projection_inventory$scenario[i]
  binary_candidates <- c(
    file.path(
      binary_map_dir,
      paste0(
        "assigned_zone_", scenario,
        "_threshold0.4_tol1e-04_novel99_maskNA8_noNovelNormal.tif"
      )
    ),
    file.path(project_root, "maps", paste0("assigned_zone_", scenario, ".tif"))
  )
  binary_file <- first_existing(binary_candidates, paste("binary map", scenario))
  compare_assigned_maps(
    binary_file,
    projection_inventory$assigned_map[i],
    scenario
  )
})

if (comparison_cached) {
  cat("[REUSE BINARY–MULTICLASS COMPARISON]", table_dir, "\n")
} else {
  fwrite(rbindlist(lapply(comparison_results, `[[`, "summary")), comparison_summary_file)
  fwrite(rbindlist(lapply(comparison_results, `[[`, "confusion")), comparison_confusion_file)
}

cat(
  "\nSECTION 7 COMPLETE\n",
  "Model: ", model_file, "\n",
  "Tables: ", table_dir, "\n",
  "Maps: ", map_dir, "\n",
  "Probabilities: ", probability_dir, "\n",
  sep = ""
)
