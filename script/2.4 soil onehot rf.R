# Soil one-vs-rest Multi-Forest models.
# Backward purging selects exactly nine predictors before final training.
# Existing selected-variable Multi-Forest models are reused by default.

library(CEMT)
library(data.table)
library(doSNOW)
library(foreach)
library(randomForest)

rm(list = ls())
gc()

find_project_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "script", "2.4 soil onehot rf.R"))) return(path)
    parent <- dirname(path)
    if (parent == path) stop("Run inside the repository or set ECOCHINA2_DIR.")
    path <- parent
  }
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "true" else "false")
  tolower(trimws(value)) %in% c("1", "true", "yes", "y")
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[1] else NA_character_
}

env_root <- Sys.getenv("ECOCHINA2_DIR", unset = "")
base_dir <- if (nzchar(env_root)) {
  normalizePath(env_root, winslash = "/", mustWork = TRUE)
} else {
  find_project_root()
}
source(file.path(base_dir, "functions", "mcRFop_cls3.R"))
source(file.path(base_dir, "functions", "proportional sampling with print.R"))

zoneID <- c(1:7, 9:50, 52:55)
base_seed <- 49L
n_selected <- 9L
n_tree_selection <- 100L
n_tree_mf <- 100L
n_forest <- 10L
absence_ratio <- 1.3
soil_predictors <- c(
  "T_GRAVEL", "T_SAND", "T_SILT", "T_CLAY", "T_REF_BULK_DENSITY",
  "T_OC", "T_PH_H2O", "T_CEC_CLAY", "T_CEC_SOIL", "T_BS", "T_TEB",
  "T_CACO3", "T_CASO4", "T_ESP", "T_ECE"
)

results_dir <- file.path(base_dir, "results")
model_dir <- file.path(base_dir, "rf_soil")
accuracy_dir <- file.path(base_dir, "accuracy_soil")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(accuracy_dir, recursive = TRUE, showWarnings = FALSE)

force_split <- env_flag("ECOCHINA2_FORCE_SPLIT", FALSE)
force_selection <- env_flag("ECOCHINA2_FORCE_SELECTION", FALSE) || force_split
force_training <- env_flag("ECOCHINA2_FORCE_SOIL_TRAINING", FALSE) ||
  force_selection

accuracy_summary_file <- file.path(accuracy_dir, "soil_rf_accuracy_summary.csv")
if (file.exists(accuracy_summary_file)) {
  existing_accuracy <- fread(accuracy_summary_file)
  if ("model" %in% names(existing_accuracy)) {
    retained_accuracy <- existing_accuracy[model == "plain_mf_rf"]
    if (nrow(retained_accuracy) != nrow(existing_accuracy)) {
      fwrite(retained_accuracy, accuracy_summary_file)
    }
  }
}

model_file <- function(zone) {
  file.path(model_dir, paste0("soil_mf_zone", zone, ".Rdata"))
}

valid_model <- function(file, object_name = "soil_mf") {
  if (!file.exists(file) || is.na(file.info(file)$size) || file.info(file)$size <= 0) {
    return(FALSE)
  }
  tryCatch({
    e <- new.env(parent = emptyenv())
    load(file, envir = e)
    if (!exists(object_name, envir = e, inherits = FALSE)) return(FALSE)
    model <- get(object_name, envir = e, inherits = FALSE)
    length(model$varlist) == n_selected
  }, error = function(e) FALSE)
}

read_n_variables <- function(path, n, available) {
  selection <- fread(path)
  variable_col <- intersect(c("variable", "Variable"), names(selection))[1]
  if (is.na(variable_col)) stop("Invalid selection file: ", path)

  variable_sets <- lapply(selection[[variable_col]], function(value) {
    if (is.na(value) || !nzchar(value) || value == "0") return(character())
    trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  })
  hit <- which(lengths(variable_sets) == n)
  if (!length(hit)) stop("No ", n, "-predictor set in ", path)

  vars <- unique(variable_sets[[hit[1]]])
  if (length(vars) != n || length(setdiff(vars, available))) {
    stop("Invalid ", n, "-predictor set in ", path)
  }
  vars
}

sample_soil <- function(x, zone_col, seed) {
  n_presence <- sum(x[[zone_col]] == 1L, na.rm = TRUE)
  if (n_presence < 2L) stop("Too few presences for ", zone_col)

  n_rounds <- if (n_presence >= 6000L) {
    3L
  } else if (n_presence >= 3000L) {
    2L
  } else {
    1L
  }
  set.seed(seed)
  sampled <- x
  for (round_index in seq_len(n_rounds)) {
    sampled <- proSysSmpl(
      sampled,
      byCol = match(zone_col, names(sampled)),
      minSz = 2000L
    )
  }
  sampled
}

multi_forest <- function(xy_y, xy_n, varlist, y_col, seed) {
  detected_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(detected_cores)) detected_cores <- 1L
  n_core <- min(max(1L, detected_cores - 1L), n_tree_mf)
  ntree_vec <- rep(n_tree_mf %/% n_core, n_core)
  if (n_tree_mf %% n_core) {
    ntree_vec[seq_len(n_tree_mf %% n_core)] <-
      ntree_vec[seq_len(n_tree_mf %% n_core)] + 1L
  }
  ntree_vec <- ntree_vec[ntree_vec > 0L]

  cl <- parallel::makeCluster(length(ntree_vec), type = "SOCK")
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterSetRNGStream(cl, iseed = seed)
  set.seed(seed)

  combined <- NULL
  n_absence <- min(nrow(xy_n), floor(nrow(xy_y) * absence_ratio))
  for (forest_index in seq_len(n_forest)) {
    train <- rbind(
      xy_y,
      xy_n[sample.int(nrow(xy_n), n_absence), , drop = FALSE]
    )
    x <- train[, varlist, drop = FALSE]
    y <- factor(train[[y_col]], levels = c(0, 1))
    forest <- foreach::foreach(
      ntree = ntree_vec,
      .combine = randomForest::combine,
      .packages = "randomForest"
    ) %dopar% {
      randomForest::randomForest(x, y, ntree = ntree, importance = TRUE)
    }
    combined <- if (is.null(combined)) forest else randomForest::combine(combined, forest)
  }
  combined
}

save_training_accuracy <- function(model, x, y, zone) {
  prediction <- factor(predict(model, x, type = "response"), levels = c(0, 1))
  confusion <- table(observed = factor(y, levels = c(0, 1)), predicted = prediction)
  fwrite(
    as.data.table(confusion),
    file.path(accuracy_dir, paste0("plain_mf_rf_confusion_train_zone", zone, ".csv"))
  )

  oob_accuracy <- NA_real_
  if (!is.null(model$confusion)) {
    oob_raw <- as.matrix(model$confusion)
    oob <- oob_raw[, setdiff(colnames(oob_raw), "class.error"), drop = FALSE]
    oob_accuracy <- sum(diag(oob)) / sum(oob)
    fwrite(
      as.data.table(oob, keep.rownames = "observed"),
      file.path(accuracy_dir, paste0("plain_mf_rf_confusion_oob_zone", zone, ".csv"))
    )
  }

  data.table(
    zone = zone,
    model = "plain_mf_rf",
    n = length(y),
    n_presence = sum(y == 1),
    n_absence = sum(y == 0),
    n_predictors = ncol(x),
    train_accuracy = sum(diag(confusion)) / sum(confusion),
    oob_accuracy = oob_accuracy
  )
}

# `plain_mf` is retained only as the historical soil-output label. These models
# use the nine backward-purged variables and are the final Multi-Forest models.
all_models_valid <- !force_split && !force_selection && !force_training && all(vapply(
  zoneID,
  function(zone) valid_model(model_file(zone)),
  logical(1)
))

if (all_models_valid) {
  cat(
    "[USE EXISTING] All 53 soil Multi-Forest models are valid.\n",
    "Model pattern: rf_soil/soil_mf_zone<ID>.Rdata\n",
    "Set ECOCHINA2_FORCE_SOIL_TRAINING=true to rebuild them.\n",
    sep = ""
  )
} else {
  env_soil <- Sys.getenv("ECOCHINA2_SOIL_TABLE", unset = "")
  soil_file <- first_existing(c(
    env_soil,
    file.path(base_dir, "data raw", "new_soil_raster.csv"),
    file.path(base_dir, "data", "raw", "new_soil_raster.csv")
  ))
  if (is.na(soil_file)) {
    stop("Missing soil table. Set ECOCHINA2_SOIL_TABLE.")
  }

  train_file <- file.path(results_dir, "soil_train_data.csv")
  coord_file <- file.path(results_dir, "soil_train_coords.csv")
  test_file <- file.path(results_dir, "soil_test_data.csv")
  onehot_file <- file.path(results_dir, "soil_hot_encoded.csv")

  split_complete <- all(file.exists(c(train_file, coord_file, test_file, onehot_file)))
  if (force_split || !split_complete) {
    soil <- fread(soil_file)
    required <- c("zoneID", soil_predictors)
    missing <- setdiff(required, names(soil))
    if (length(missing)) stop("Missing soil columns: ", paste(missing, collapse = ", "))

    coordinate_names <- intersect(c("x", "y"), names(soil))
    if (length(coordinate_names) != 2L) stop("Soil table must include x and y.")
    soil <- soil[complete.cases(soil[, ..required])]
    set.seed(base_seed)
    train_id <- sample.int(nrow(soil), floor(0.70 * nrow(soil)))

    soil_train <- soil[train_id, ..required]
    soil_test <- soil[-train_id, ..required]
    soil_coords <- soil[train_id, ..coordinate_names]
    fwrite(soil_train, train_file)
    fwrite(soil_coords, coord_file)
    fwrite(soil_test, test_file)

    soil_hot <- cbind(soil_train, soil_coords)
    for (zone in sort(unique(soil_train$zoneID))) {
      soil_hot[, (paste0("zone", zone)) := as.integer(zoneID == zone)]
    }
    fwrite(soil_hot, onehot_file)
    rm(soil, soil_train, soil_test, soil_coords)
    gc()
  } else {
    cat("[USE EXISTING] Soil train/test split and one-hot table.\n")
    soil_hot <- fread(onehot_file)
  }

  for (zone in zoneID) {
    zone_col <- paste0("zone", zone)
    if (!zone_col %in% names(soil_hot)) {
      soil_hot[, (zone_col) := as.integer(zoneID == zone)]
    }

    selection_file <- file.path(results_dir, paste0("soilopList_zone", zone, ".csv"))
    if (!force_selection && file.exists(selection_file)) {
      read_n_variables(selection_file, n_selected, soil_predictors)
      cat("[USE SELECTION] zone", zone, "\n")
    } else {
      cat("[SELECT] zone", zone, "\n")
      sampled <- as.data.frame(sample_soil(soil_hot, zone_col, base_seed + zone))
      selection <- mcRFop_cls(
        sampled[, soil_predictors, drop = FALSE],
        factor(sampled[[zone_col]], levels = c(0, 1)),
        nTree = n_tree_selection,
        seed = base_seed + zone
      )
      fwrite(selection, selection_file)
      read_n_variables(selection_file, n_selected, soil_predictors)
      rm(sampled, selection)
      gc()
    }
  }

  accuracy_rows <- list()
  for (zone in zoneID) {
    output_file <- model_file(zone)
    if (!force_training && valid_model(output_file)) {
      cat("[USE MODEL] zone", zone, "\n")
      next
    }

    zone_col <- paste0("zone", zone)
    selection_file <- file.path(results_dir, paste0("soilopList_zone", zone, ".csv"))
    varlist <- read_n_variables(selection_file, n_selected, soil_predictors)
    cat("[TRAIN] zone", zone, "|", paste(varlist, collapse = ", "), "\n")

    sampled <- as.data.frame(sample_soil(soil_hot, zone_col, base_seed + zone))
    sampled <- sampled[complete.cases(sampled[, c(zone_col, varlist)]), , drop = FALSE]
    n_positive <- sum(sampled[[zone_col]] == 1L)
    negative <- which(sampled[[zone_col]] == 0L)
    if (length(negative) > n_positive * absence_ratio) {
      set.seed(base_seed + zone)
      keep_negative <- sample(
        negative,
        min(length(negative), round(n_positive * absence_ratio))
      )
      sampled <- sampled[c(which(sampled[[zone_col]] == 1L), keep_negative), , drop = FALSE]
    }

    xy_y <- sampled[sampled[[zone_col]] == 1L, , drop = FALSE]
    xy_n <- sampled[sampled[[zone_col]] == 0L, , drop = FALSE]
    if (nrow(xy_y) < 2L || nrow(xy_n) < 2L) stop("Too few rows for zone ", zone)

    soil_mf <- multi_forest(
      xy_y,
      xy_n,
      varlist = varlist,
      y_col = zone_col,
      seed = base_seed + zone
    )
    soil_mf$varlist <- varlist
    soil_mf$zoneID <- zone
    soil_mf$model <- "plain_mf"
    soil_mf$nForest <- n_forest
    soil_mf$nTree_per_forest <- n_tree_mf
    save(soil_mf, file = output_file)

    accuracy_rows[[length(accuracy_rows) + 1L]] <- save_training_accuracy(
      soil_mf,
      sampled[, varlist, drop = FALSE],
      factor(sampled[[zone_col]], levels = c(0, 1)),
      zone
    )
    rm(sampled, xy_y, xy_n, soil_mf)
    gc()
  }

  if (length(accuracy_rows)) {
    updated <- rbindlist(accuracy_rows)
    if (file.exists(accuracy_summary_file) && !force_training) {
      previous <- fread(accuracy_summary_file)
      if ("model" %in% names(previous) && "zone" %in% names(previous)) {
        previous <- previous[
          model == "plain_mf_rf" & !(zone %in% updated$zone)
        ]
      } else {
        previous <- data.table()
      }
      updated <- rbind(previous, updated, fill = TRUE)
    }
    setorder(updated, model, zone)
    fwrite(updated, accuracy_summary_file)
  }

  cat("COMPLETE: selected-variable soil Multi-Forest workflow.\n")
}
