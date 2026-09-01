# Climate one-vs-rest Multi-Forest models.
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
    if (file.exists(file.path(path, "script", "2.15 climate onehot rf.R"))) return(path)
    parent <- dirname(path)
    if (parent == path) stop("Run inside the repository or set ECOCHINA2_DIR.")
    path <- parent
  }
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "true" else "false")
  tolower(trimws(value)) %in% c("1", "true", "yes", "y")
}

env_root <- Sys.getenv("ECOCHINA2_DIR", unset = "")
base_dir <- if (nzchar(env_root)) {
  normalizePath(env_root, winslash = "/", mustWork = TRUE)
} else {
  find_project_root()
}
source(file.path(base_dir, "functions", "mcRFop_cls3.R"))

zoneID <- c(1:7, 9:50, 52:55)
base_seed <- 49L
n_selected <- 9L
n_tree_selection <- 100L
n_tree_mf <- 100L
n_forest <- 10L
absence_ratio <- 1.3

results_dir <- file.path(base_dir, "results")
model_dir <- file.path(base_dir, "rf")
accuracy_dir <- file.path(base_dir, "accuracy_climate_var")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(accuracy_dir, recursive = TRUE, showWarnings = FALSE)

force_split <- env_flag("ECOCHINA2_FORCE_SPLIT", FALSE)
force_selection <- env_flag("ECOCHINA2_FORCE_SELECTION", FALSE) || force_split
force_training <- env_flag("ECOCHINA2_FORCE_CLIMATE_TRAINING", FALSE) ||
  force_selection

accuracy_summary_file <- file.path(
  accuracy_dir, "climate_rf_var_accuracy_summary.csv"
)
if (file.exists(accuracy_summary_file)) {
  existing_accuracy <- fread(accuracy_summary_file)
  if ("model" %in% names(existing_accuracy)) {
    retained_accuracy <- existing_accuracy[model == "mf_var"]
    if (nrow(retained_accuracy) != nrow(existing_accuracy)) {
      fwrite(retained_accuracy, accuracy_summary_file)
    }
  }
}

model_file <- function(zone) {
  file.path(model_dir, paste0("clm_mfVar_zone", zone, ".Rdata"))
}

valid_model <- function(file, object_name = "clm_mfVar") {
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

sample_climate <- function(x, zone_col, predictors, max_presence, seed) {
  n_presence <- sum(x[[zone_col]] == 1L, na.rm = TRUE)
  if (n_presence < 2L) stop("Too few presences for ", zone_col)

  pos <- min(max_presence, n_presence)
  noise <- max(1L, min(as.integer(round(0.10 * pos)), pos - 1L))
  smpl_pa(
    x,
    zone_col,
    cols = c(zone_col, predictors),
    pos = pos,
    noise = noise,
    pa = 1 / absence_ratio,
    max_n = 20000L,
    seed = seed
  )
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
    file.path(accuracy_dir, paste0("mf_var_zone", zone, "_confusion_train.csv"))
  )

  oob_accuracy <- NA_real_
  if (!is.null(model$confusion)) {
    oob <- as.matrix(model$confusion[, c("0", "1"), drop = FALSE])
    oob_accuracy <- sum(diag(oob)) / sum(oob)
    fwrite(
      as.data.table(oob, keep.rownames = "observed"),
      file.path(accuracy_dir, paste0("mf_var_zone", zone, "_confusion_oob.csv"))
    )
  }

  data.table(
    zone = zone,
    model = "mf_var",
    n = length(y),
    oob_accuracy = oob_accuracy,
    train_accuracy = sum(diag(confusion)) / sum(confusion)
  )
}

# Avoid reading the multi-million-row climate table when every final model exists.
all_models_valid <- !force_split && !force_selection && !force_training && all(vapply(
  zoneID,
  function(zone) valid_model(model_file(zone)),
  logical(1)
))

if (all_models_valid) {
  cat(
    "[USE EXISTING] All 53 climate Multi-Forest models are valid.\n",
    "Model pattern: rf/clm_mfVar_zone<ID>.Rdata\n",
    "Set ECOCHINA2_FORCE_CLIMATE_TRAINING=true to rebuild them.\n",
    sep = ""
  )
} else {
  climate_file <- file.path(
    base_dir,
    "data raw",
    "1. zoneID_Clm_800m_Normal_1961_1990SY.csv"
  )
  train_file <- file.path(results_dir, "train_data.csv")
  test_file <- file.path(results_dir, "test_data.csv")
  onehot_file <- file.path(results_dir, "train_combined_data_onehot.csv")

  split_complete <- all(file.exists(c(train_file, test_file, onehot_file)))
  if (force_split || !split_complete) {
    if (!file.exists(climate_file)) stop("Run script/1. data.R first.")

    climate <- fread(climate_file)
    climate[climate == -9999] <- NA
    if (!"zoneID" %in% names(climate)) stop("Climate table lacks zoneID.")

    predictors <- names(climate)[6:ncol(climate)]
    climate <- climate[complete.cases(climate[, c("zoneID", predictors), with = FALSE])]
    set.seed(base_seed)
    train_id <- sample.int(nrow(climate), floor(0.70 * nrow(climate)))
    climate_train <- climate[train_id]
    climate_test <- climate[-train_id]
    fwrite(climate_train, train_file)
    fwrite(climate_test, test_file)

    for (zone in sort(unique(climate_train$zoneID))) {
      climate_train[, (paste0("zone", zone)) := as.integer(zoneID == zone)]
    }
    fwrite(climate_train, onehot_file)
    rm(climate, climate_test)
    gc()
  } else {
    cat("[USE EXISTING] Climate train/test split and one-hot table.\n")
    climate_train <- fread(onehot_file)
    train_header <- fread(train_file, nrows = 0)
    predictors <- names(train_header)[6:ncol(train_header)]
  }

  if (!exists("predictors")) {
    train_header <- fread(train_file, nrows = 0)
    predictors <- names(train_header)[6:ncol(train_header)]
  }

  for (zone in zoneID) {
    zone_col <- paste0("zone", zone)
    if (!zone_col %in% names(climate_train)) {
      climate_train[, (zone_col) := as.integer(zoneID == zone)]
    }

    selection_file <- file.path(results_dir, paste0("clmhot_opList_zone", zone, ".csv"))
    if (!force_selection && file.exists(selection_file)) {
      read_n_variables(selection_file, n_selected, predictors)
      cat("[USE SELECTION] zone", zone, "\n")
    } else {
      cat("[SELECT] zone", zone, "\n")
      sampled <- as.data.frame(sample_climate(
        climate_train,
        zone_col,
        predictors,
        max_presence = 5000L,
        seed = base_seed + zone
      ))
      selection <- mcRFop_cls(
        sampled[, predictors, drop = FALSE],
        factor(sampled[[zone_col]], levels = c(0, 1)),
        nTree = n_tree_selection,
        seed = base_seed + zone
      )
      fwrite(selection, selection_file)
      read_n_variables(selection_file, n_selected, predictors)
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
    selection_file <- file.path(results_dir, paste0("clmhot_opList_zone", zone, ".csv"))
    varlist <- read_n_variables(selection_file, n_selected, predictors)
    cat("[TRAIN] zone", zone, "|", paste(varlist, collapse = ", "), "\n")

    sampled <- as.data.frame(sample_climate(
      climate_train,
      zone_col,
      predictors,
      max_presence = 8000L,
      seed = base_seed + zone
    ))
    sampled <- sampled[complete.cases(sampled[, c(zone_col, varlist)]), , drop = FALSE]
    xy_y <- sampled[sampled[[zone_col]] == 1L, , drop = FALSE]
    xy_n <- sampled[sampled[[zone_col]] == 0L, , drop = FALSE]
    if (nrow(xy_y) < 2L || nrow(xy_n) < 2L) stop("Too few rows for zone ", zone)

    clm_mfVar <- multi_forest(
      xy_y,
      xy_n,
      varlist = varlist,
      y_col = zone_col,
      seed = base_seed + zone
    )
    clm_mfVar$varlist <- varlist
    clm_mfVar$zoneID <- zone
    clm_mfVar$model <- "mf_var"
    clm_mfVar$nForest <- n_forest
    clm_mfVar$nTree_per_forest <- n_tree_mf
    save(clm_mfVar, file = output_file)

    accuracy_rows[[length(accuracy_rows) + 1L]] <- save_training_accuracy(
      clm_mfVar,
      sampled[, varlist, drop = FALSE],
      factor(sampled[[zone_col]], levels = c(0, 1)),
      zone
    )
    rm(sampled, xy_y, xy_n, clm_mfVar)
    gc()
  }

  if (length(accuracy_rows)) {
    updated <- rbindlist(accuracy_rows)
    if (file.exists(accuracy_summary_file) && !force_training) {
      previous <- fread(accuracy_summary_file)
      if ("model" %in% names(previous) && "zone" %in% names(previous)) {
        previous <- previous[
          model == "mf_var" & !(zone %in% updated$zone)
        ]
      } else {
        previous <- data.table()
      }
      updated <- rbind(previous, updated, fill = TRUE)
    }
    setorder(updated, model, zone)
    fwrite(updated, accuracy_summary_file)
  }

  cat("COMPLETE: selected-variable climate Multi-Forest workflow.\n")
}
