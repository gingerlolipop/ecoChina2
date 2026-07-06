# Section 7. Multiclass RF: climate + soil
# ============================================================
# Purpose:
#   1. Use mcRFop backward elimination to select variables.
#   2. Fit one ordinary multiclass random forest for all 53 zones.
#   3. Use selected climate and soil predictors in the same model.
#   4. Evaluate the model on a stratified independent test set.
#   5. Predict the normal-period vegetation-zone map.
#   6. Compare the multiclass map with the four existing overlay maps.
#
# Notes:
#   - This is one standard randomForest model, not the four binary-RF variants.
#   - All predictors used by the existing plain climate and plain soil RFs are used.
#   - Class-balanced sampling is applied within each tree because zone sizes differ greatly.
#   - Normal-map agreement is reconstruction agreement, not independent validation.
# ============================================================

library(terra)
library(randomForest)
library(data.table)
library(foreach)
library(doSNOW)

rm(list = ls())
gc()

# 0. Parameters ================================================================

base_dir <- "H:/Jing/ecoChina2"

source(file.path(base_dir, "functions", "mcRFop_cls3.R"))

r_file <- file.path(base_dir, "raster/ecosys_ori.tif")
clm_data_file <- file.path(
  base_dir,
  "data raw/1. zoneID_Clm_800m_Normal_1961_1990SY.csv"
)

climdir_normal <- paste0(
  "H:/Jing/ecoChina/play/China/ClimateData/CN/",
  "800m/Normal_1961_1990"
)

soil_dir <- "H:/Jing/soil rasters/tif2"

result_dir <- file.path(base_dir, "results")
model_dir <- file.path(base_dir, "rf_multiclass")
map_dir <- file.path(base_dir, "result maps", "multiclass_rf")
assess_dir <- file.path(base_dir, "assessment", "multiclass_rf")
tmp_dir <- file.path(base_dir, "tmp_multiclass_rf")

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(assess_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

terraOptions(tempdir = tmp_dir, memfrac = 0.25)

zones <- c(1:7, 9:50, 52:55)

base_seed <- 49L
ntree <- 500L

# Backward variable selection before the final RF.
op_ntree <- 100L
op_nrep <- 10L
op_per_zone <- 1000L
op_accuracy_tolerance <- 0.001

# Maximum observations retained after complete-case filtering.
# Smaller zones use all available observations and keep a 70:30 split.
candidate_per_zone <- 10000L
max_train_per_zone <- 5000L
max_test_per_zone <- 2000L
train_fraction <- 0.70

predict_cores <- min(
  4L,
  max(1L, parallel::detectCores() - 1L)
)

model_file <- file.path(model_dir, "multiclass_climate_soil_rf.Rdata")
map_file <- file.path(map_dir, "assigned_zone_normal_multiclass_rf.tif")

set.seed(base_seed)


# 1. Helper functions ==========================================================

# mcRFop_cls3.R was written for binary RFs. In a binary randomForest,
# importance columns 3:4 are MeanDecreaseAccuracy and MeanDecreaseGini.
# In a 53-class RF, columns 3:4 are class-specific importance values.
# This multiclass version keeps the same backward-elimination logic but
# explicitly uses the two global importance columns.
mcRFop_cls_multiclass <- function(
    x,
    y,
    nTree = op_ntree,
    nRep = op_nrep) {
  
  x <- as.data.frame(x)
  y <- droplevels(as.factor(y))
  
  if (ncol(x) < 3) {
    stop("At least three predictors are required for variable selection.")
  }
  
  n_core <- min(
    max(1L, parallel::detectCores() - 1L),
    nTree
  )
  
  ntree_vec <- rep(floor(nTree / n_core), n_core)
  
  if (nTree %% n_core > 0) {
    ntree_vec[seq_len(nTree %% n_core)] <-
      ntree_vec[seq_len(nTree %% n_core)] + 1L
  }
  
  ntree_vec <- ntree_vec[ntree_vec > 0]
  
  cl <- makeCluster(length(ntree_vec), type = "SOCK")
  registerDoSNOW(cl)
  on.exit(stopCluster(cl), add = TRUE)
  
  op_list <- matrix(
    NA_character_,
    nrow = ncol(x),
    ncol = 2,
    dimnames = list(
      NULL,
      c("Accy", "variable")
    )
  )
  
  x_current <- x
  
  while (ncol(x_current) >= 3) {
    cat(
      "[OPLIST] predictors:", ncol(x_current),
      "| trees:", nTree,
      "| repeats:", nRep,
      "\n"
    )
    
    imp_accum <- NULL
    acc_rep <- numeric(nRep)
    
    for (rr in seq_len(nRep)) {
      set.seed(base_seed + 10000L + ncol(x_current) * 100L + rr)
      
      rf_rep <- foreach(
        ntree_i = ntree_vec,
        .combine = combine,
        .packages = "randomForest"
      ) %dopar% {
        randomForest(
          x_current,
          y,
          ntree = ntree_i,
          importance = TRUE
        )
      }
      
      imp_rep <- importance(rf_rep)
      
      required_cols <- c(
        "MeanDecreaseAccuracy",
        "MeanDecreaseGini"
      )
      
      if (!all(required_cols %in% colnames(imp_rep))) {
        stop(
          "Global multiclass importance columns were not found."
        )
      }
      
      imp_rep <- imp_rep[, required_cols, drop = FALSE]
      
      if (is.null(imp_accum)) {
        imp_accum <- imp_rep
      } else {
        imp_accum <- imp_accum + imp_rep
      }
      
      pred_rep <- predict(rf_rep, x_current)
      acc_rep[rr] <- mean(pred_rep == y)
    }
    
    imp_mean <- imp_accum / nRep
    imp_mean <- imp_mean[
      order(imp_mean[, "MeanDecreaseAccuracy"]),
      ,
      drop = FALSE
    ]
    
    n_var <- nrow(imp_mean)
    
    op_list[n_var, "Accy"] <- sprintf(
      "%.10f",
      mean(acc_rep)
    )
    
    op_list[n_var, "variable"] <- paste(
      rownames(imp_mean),
      collapse = ", "
    )
    
    if (n_var <= 4) {
      break
    }
    
    remove_vars <- rownames(imp_mean)[1:2]
    
    x_current <- x_current[
      ,
      setdiff(names(x_current), remove_vars),
      drop = FALSE
    ]
  }
  
  as.data.frame(
    op_list,
    stringsAsFactors = FALSE
  )
}

select_op_variables <- function(
    op_list,
    tolerance = op_accuracy_tolerance) {
  
  op <- as.data.table(op_list)
  
  op[, row_id := seq_len(.N)]
  op[, accuracy := suppressWarnings(as.numeric(Accy))]
  op[, variable := trimws(variable)]
  op <- op[
    is.finite(accuracy) &
      !is.na(variable) &
      nzchar(variable)
  ]
  
  if (!nrow(op)) {
    stop("mcRFop did not return a valid variable set.")
  }
  
  op[, n_variables := lengths(strsplit(variable, ","))]
  
  best_accuracy <- max(op$accuracy)
  
  candidates <- op[
    accuracy >= best_accuracy - tolerance
  ]
  
  # Prefer the smallest set whose accuracy is effectively tied with the best.
  setorder(
    candidates,
    n_variables,
    -accuracy
  )
  
  selected <- trimws(
    unlist(
      strsplit(candidates$variable[1], ",")
    )
  )
  
  list(
    variables = selected,
    selected_row = candidates[1],
    best_accuracy = best_accuracy,
    all_candidates = op
  )
}

div <- function(a, b) {
  ifelse(is.finite(b) & b > 0, a / b, NA_real_)
}

mean_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

auc_rank <- function(y, p) {
  n1 <- as.numeric(sum(y == 1))
  n0 <- as.numeric(sum(y == 0))
  
  if (n1 == 0 || n0 == 0) return(NA_real_)
  
  (
    sum(rank(p, ties.method = "average")[y == 1]) -
      n1 * (n1 + 1) / 2
  ) / (n1 * n0)
}

remove_raster_files <- function(filepath) {
  base <- tools::file_path_sans_ext(filepath)
  
  files <- unique(c(
    filepath,
    paste0(filepath, ".aux.xml"),
    paste0(filepath, ".ovr"),
    paste0(filepath, ".msk"),
    paste0(base, ".aux.xml"),
    paste0(base, ".ovr"),
    paste0(base, ".tfw")
  ))
  
  suppressWarnings(unlink(files[file.exists(files)], force = TRUE))
  invisible(TRUE)
}

read_predictor_names <- function() {
  clm_ref_file <- file.path(result_dir, "train_data.csv")
  soil_ref_file <- file.path(result_dir, "soil_train_data.csv")
  
  if (!file.exists(clm_ref_file)) {
    stop("Missing climate reference table: ", clm_ref_file)
  }
  
  if (!file.exists(soil_ref_file)) {
    stop("Missing soil reference table: ", soil_ref_file)
  }
  
  clm_names <- names(fread(clm_ref_file, nrows = 0))
  soil_names <- names(fread(soil_ref_file, nrows = 0))
  
  if (length(clm_names) < 6) {
    stop("Unexpected climate reference-table structure.")
  }
  
  clm_vars <- clm_names[6:length(clm_names)]
  soil_vars <- setdiff(soil_names, "zoneID")
  
  list(
    climate = clm_vars,
    soil_raw = soil_vars,
    soil_model = paste0("soil_", soil_vars)
  )
}

get_stack <- function(varlist, raster_dir, layer_names = varlist) {
  files <- file.path(raster_dir, paste0(varlist, ".tif"))
  missing_files <- files[!file.exists(files)]
  
  if (length(missing_files) > 0) {
    stop(
      "Missing predictor rasters:\n",
      paste(missing_files, collapse = "\n")
    )
  }
  
  s <- rast(files)
  names(s) <- layer_names
  s
}

extract_without_id <- function(x, points) {
  out <- as.data.table(terra::extract(x, points))
  
  if ("ID" %in% names(out)) {
    out[, ID := NULL]
  }
  
  out
}

rf_class <- function(model, data) {
  pred <- predict(
    model,
    as.data.frame(data),
    type = "response"
  )
  
  as.integer(as.character(pred))
}

assess_multiclass <- function(observed, predicted, probability, zone_values) {
  observed <- factor(observed, levels = as.character(zone_values))
  predicted <- factor(predicted, levels = as.character(zone_values))
  
  zone_metrics <- rbindlist(
    lapply(zone_values, function(z) {
      z_chr <- as.character(z)
      
      y <- as.integer(observed == z_chr)
      p <- as.integer(predicted == z_chr)
      
      TP <- sum(y == 1 & p == 1)
      TN <- sum(y == 0 & p == 0)
      FP <- sum(y == 0 & p == 1)
      FN <- sum(y == 1 & p == 0)
      
      recall <- div(TP, TP + FN)
      specificity <- div(TN, TN + FP)
      precision <- div(TP, TP + FP)
      
      prob_z <- probability[, z_chr]
      
      data.table(
        zone = z,
        n_test = TP + FN,
        TP, TN, FP, FN,
        accuracy = div(TP + TN, TP + TN + FP + FN),
        balanced_accuracy = div(recall + specificity, 2),
        recall,
        specificity,
        precision,
        f1 = div(2 * precision * recall, precision + recall),
        tss = recall + specificity - 1,
        auc = auc_rank(y, prob_z)
      )
    })
  )
  
  overall <- data.table(
    method = "multiclass_rf",
    n_test = length(observed),
    overall_accuracy = mean(observed == predicted),
    macro_balanced_accuracy = mean_na(zone_metrics$balanced_accuracy),
    macro_recall = mean_na(zone_metrics$recall),
    macro_specificity = mean_na(zone_metrics$specificity),
    macro_precision = mean_na(zone_metrics$precision),
    macro_f1 = mean_na(zone_metrics$f1),
    macro_tss = mean_na(zone_metrics$tss),
    macro_auc = mean_na(zone_metrics$auc)
  )
  
  list(
    overall = overall,
    zone = zone_metrics,
    confusion = as.data.table(
      table(
        observed = observed,
        predicted = predicted
      )
    )
  )
}

assess_map <- function(
    original,
    predicted,
    method,
    zone_values,
    category_lut = NULL) {
  
  pred <- subst(
    predicted,
    from = zone_values,
    to = zone_values,
    others = NA
  )
  names(pred) <- "pred"
  
  valid_original <- global(
    !is.na(original),
    "sum",
    na.rm = TRUE
  )[1, 1]
  
  compared <- global(
    !is.na(original) & !is.na(pred),
    "sum",
    na.rm = TRUE
  )[1, 1]
  
  ct <- as.data.table(
    crosstab(
      c(original, pred),
      long = TRUE,
      useNA = FALSE
    )
  )
  
  setnames(ct, c("ori", "pred", "n"))
  
  ct[, `:=`(
    method = method,
    ori = as.integer(ori),
    pred = as.integer(pred),
    n = as.numeric(n)
  )]
  
  total <- sum(ct$n, na.rm = TRUE)
  correct <- ct[ori == pred, sum(n, na.rm = TRUE)]
  
  zone_metrics <- rbindlist(
    lapply(zone_values, function(z) {
      TP <- ct[ori == z & pred == z, sum(n, na.rm = TRUE)]
      FN <- ct[ori == z & pred != z, sum(n, na.rm = TRUE)]
      FP <- ct[ori != z & pred == z, sum(n, na.rm = TRUE)]
      TN <- total - TP - FN - FP
      
      recall <- div(TP, TP + FN)
      specificity <- div(TN, TN + FP)
      precision <- div(TP, TP + FP)
      
      data.table(
        method,
        zone = z,
        original_pixels = TP + FN,
        predicted_pixels = TP + FP,
        predicted_to_original = div(TP + FP, TP + FN),
        TP, TN, FP, FN,
        accuracy = div(TP + TN, total),
        balanced_accuracy = div(recall + specificity, 2),
        recall,
        specificity,
        precision,
        f1 = div(2 * precision * recall, precision + recall),
        tss = recall + specificity - 1
      )
    })
  )
  
  broad_accuracy <- NA_real_
  
  if (!is.null(category_lut)) {
    ct[, ori_category := category_lut[as.character(ori)]]
    ct[, pred_category := category_lut[as.character(pred)]]
    
    broad_accuracy <- div(
      ct[
        !is.na(ori_category) &
          !is.na(pred_category) &
          ori_category == pred_category,
        sum(n)
      ],
      total
    )
  }
  
  overall <- data.table(
    method,
    valid_original_pixels = valid_original,
    compared_pixels = compared,
    missing_predictions = valid_original - compared,
    coverage = div(compared, valid_original),
    exact_zone_accuracy = div(correct, total),
    broad_category_accuracy = broad_accuracy,
    macro_balanced_accuracy = mean_na(zone_metrics$balanced_accuracy),
    macro_recall = mean_na(zone_metrics$recall),
    macro_specificity = mean_na(zone_metrics$specificity),
    macro_precision = mean_na(zone_metrics$precision),
    macro_f1 = mean_na(zone_metrics$f1),
    macro_tss = mean_na(zone_metrics$tss)
  )
  
  list(
    overall = overall,
    zone = zone_metrics,
    confusion = ct
  )
}


# 2. Predictor names and raster stacks =========================================

predictor_info <- read_predictor_names()

clm_vars <- predictor_info$climate
soil_vars <- predictor_info$soil_raw
soil_model_vars <- predictor_info$soil_model
predictor_vars <- c(clm_vars, soil_model_vars)

if (anyDuplicated(predictor_vars)) {
  stop("Duplicated predictor names after adding the soil_ prefix.")
}

cat("Climate predictors:", length(clm_vars), "\n")
cat("Soil predictors:", length(soil_vars), "\n")
cat("Total predictors:", length(predictor_vars), "\n")

r <- rast(r_file)

clm_stack <- get_stack(
  clm_vars,
  climdir_normal,
  layer_names = clm_vars
)

soil_stack <- get_stack(
  soil_vars,
  soil_dir,
  layer_names = soil_model_vars
)


# 3. Stratified training and testing samples ===================================

if (!file.exists(clm_data_file)) {
  stop("Missing cell-zone table: ", clm_data_file)
}

cat("\n[READ CELL INDEX]\n")

cell_dt <- fread(
  clm_data_file,
  select = c("cell", "zoneID")
)

cell_dt[, `:=`(
  cell = as.integer(cell),
  zoneID = as.integer(as.character(zoneID))
)]

cell_dt <- cell_dt[
  zoneID %in% zones &
    !is.na(cell)
]

missing_zone <- setdiff(zones, unique(cell_dt$zoneID))

if (length(missing_zone) > 0) {
  stop(
    "Modeled zones missing from the cell table: ",
    paste(missing_zone, collapse = ", ")
  )
}

set.seed(base_seed)

sample_cells <- cell_dt[
  ,
  {
    n_take <- min(.N, candidate_per_zone)
    .SD[sample.int(.N, n_take)]
  },
  by = zoneID
]

rm(cell_dt)
gc()

xy <- xyFromCell(r, sample_cells$cell)

points <- vect(
  data.frame(
    x = xy[, 1],
    y = xy[, 2]
  ),
  geom = c("x", "y"),
  crs = crs(r)
)

points_clm <- if (
  same.crs(points, clm_stack)
) {
  points
} else {
  project(points, crs(clm_stack))
}

points_soil <- if (
  same.crs(points, soil_stack)
) {
  points
} else {
  project(points, crs(soil_stack))
}

cat("\n[EXTRACT CLIMATE]\n")
clm_values <- extract_without_id(clm_stack, points_clm)

cat("\n[EXTRACT SOIL]\n")
soil_values <- extract_without_id(soil_stack, points_soil)

sample_data <- cbind(
  sample_cells[, .(cell, zoneID)],
  clm_values,
  soil_values
)

rm(points, points_clm, points_soil, xy, clm_values, soil_values, sample_cells)
gc()

sample_data <- sample_data[
  complete.cases(sample_data[, ..predictor_vars])
]

available_count <- sample_data[
  ,
  .(complete_observations = .N),
  by = zoneID
][order(zoneID)]

fwrite(
  available_count,
  file.path(assess_dir, "multiclass_complete_samples_by_zone.csv")
)

missing_zone <- setdiff(zones, unique(sample_data$zoneID))

if (length(missing_zone) > 0) {
  stop(
    "No complete climate-soil observations for zones: ",
    paste(missing_zone, collapse = ", ")
  )
}

train_list <- vector("list", length(zones))
test_list <- vector("list", length(zones))

for (k in seq_along(zones)) {
  z <- zones[k]
  d <- sample_data[zoneID == z]
  
  set.seed(base_seed + z)
  idx <- sample.int(nrow(d))
  
  n_train <- min(
    max_train_per_zone,
    floor(train_fraction * nrow(d))
  )
  
  n_test <- min(
    max_test_per_zone,
    nrow(d) - n_train
  )
  
  if (n_train < 2 || n_test < 1) {
    stop(
      "Insufficient complete observations for zone ",
      z,
      ": train = ", n_train,
      ", test = ", n_test
    )
  }
  
  train_idx <- idx[seq_len(n_train)]
  test_idx <- idx[n_train + seq_len(n_test)]
  
  train_list[[k]] <- d[train_idx]
  test_list[[k]] <- d[test_idx]
}

train_data <- rbindlist(train_list)
test_data <- rbindlist(test_list)

rm(train_list, test_list, sample_data)
gc()

setorder(train_data, zoneID, cell)
setorder(test_data, zoneID, cell)

sample_summary <- merge(
  train_data[, .(train_n = .N), by = zoneID],
  test_data[, .(test_n = .N), by = zoneID],
  by = "zoneID",
  all = TRUE
)

fwrite(
  sample_summary,
  file.path(assess_dir, "multiclass_train_test_counts.csv")
)

fwrite(
  train_data[, .(cell, zoneID)],
  file.path(assess_dir, "multiclass_train_cells.csv")
)

fwrite(
  test_data[, .(cell, zoneID)],
  file.path(assess_dir, "multiclass_test_cells.csv")
)

cat("\nTraining observations:", nrow(train_data), "\n")
cat("Testing observations:", nrow(test_data), "\n")


# 4. Pre-training variable selection with mcRFop ==============================

op_train <- train_data[
  ,
  {
    n_take <- min(.N, op_per_zone)
    .SD[sample.int(.N, n_take)]
  },
  by = zoneID
]

x_op <- as.data.frame(
  op_train[, ..predictor_vars]
)

y_op <- factor(
  op_train$zoneID,
  levels = as.character(zones)
)

cat("\n[PRE-TRAIN VARIABLE SELECTION]\n")
cat("Observations:", nrow(op_train), "\n")
cat("Candidate predictors:", length(predictor_vars), "\n")

set.seed(base_seed)

multiclass_opList <- mcRFop_cls_multiclass(
  x = x_op,
  y = y_op,
  nTree = op_ntree,
  nRep = op_nrep
)

fwrite(
  as.data.table(
    multiclass_opList,
    keep.rownames = "row"
  ),
  file.path(
    assess_dir,
    "multiclass_mcRFop_variable_selection.csv"
  )
)

op_selected <- select_op_variables(
  multiclass_opList,
  tolerance = op_accuracy_tolerance
)

selected_vars <- op_selected$variables
selected_clm_vars <- intersect(selected_vars, clm_vars)
selected_soil_vars <- intersect(selected_vars, soil_model_vars)

if (
  !length(selected_vars) ||
  !all(selected_vars %in% predictor_vars)
) {
  stop("Invalid predictor set returned by mcRFop.")
}

selected_summary <- data.table(
  selected_accuracy = op_selected$selected_row$accuracy,
  best_accuracy = op_selected$best_accuracy,
  accuracy_tolerance = op_accuracy_tolerance,
  n_selected = length(selected_vars),
  n_climate_selected = length(selected_clm_vars),
  n_soil_selected = length(selected_soil_vars),
  variables = paste(selected_vars, collapse = ", ")
)

fwrite(
  selected_summary,
  file.path(
    assess_dir,
    "multiclass_selected_variables_summary.csv"
  )
)

fwrite(
  data.table(
    variable = selected_vars,
    niche = ifelse(
      selected_vars %in% selected_clm_vars,
      "climate",
      "soil"
    )
  ),
  file.path(
    assess_dir,
    "multiclass_selected_variables.csv"
  )
)

cat("Selected predictors:", length(selected_vars), "\n")
cat("Selected climate predictors:", length(selected_clm_vars), "\n")
cat("Selected soil predictors:", length(selected_soil_vars), "\n")
cat("Selected accuracy:", op_selected$selected_row$accuracy, "\n")
cat("Best accuracy:", op_selected$best_accuracy, "\n")
cat("Variables:", paste(selected_vars, collapse = ", "), "\n")

rm(
  x_op,
  y_op,
  op_train,
  multiclass_opList,
  op_selected
)
gc()


# 5. Fit one ordinary multiclass random forest ================================

x_train <- as.data.frame(train_data[, ..selected_vars])
y_train <- factor(
  train_data$zoneID,
  levels = as.character(zones)
)

class_n <- table(y_train)

if (any(class_n == 0)) {
  stop(
    "Training data are missing zones: ",
    paste(names(class_n)[class_n == 0], collapse = ", ")
  )
}

# Equal class samples are drawn for every tree.
samp_per_class <- min(as.integer(class_n))
sampsize <- rep(samp_per_class, length(class_n))
names(sampsize) <- names(class_n)

mtry_value <- floor(sqrt(length(selected_vars)))

cat("\n[FIT FINAL MULTICLASS RF]\n")
cat("ntree:", ntree, "\n")
cat("mtry:", mtry_value, "\n")
cat("selected predictors:", length(selected_vars), "\n")
cat("samples per class per tree:", samp_per_class, "\n")

set.seed(base_seed)

multiclass_rf <- randomForest(
  x = x_train,
  y = y_train,
  ntree = ntree,
  mtry = mtry_value,
  importance = TRUE,
  replace = TRUE,
  sampsize = sampsize,
  keep.forest = TRUE
)

multiclass_rf$varlist <- selected_vars
multiclass_rf$all_candidate_vars <- predictor_vars
multiclass_rf$climate_vars <- selected_clm_vars
multiclass_rf$soil_vars <- selected_soil_vars
multiclass_rf$zones <- zones
multiclass_rf$base_seed <- base_seed
multiclass_rf$train_cells <- train_data$cell
multiclass_rf$test_cells <- test_data$cell

save(
  multiclass_rf,
  file = model_file
)

importance_dt <- as.data.table(
  importance(multiclass_rf),
  keep.rownames = "variable"
)

fwrite(
  importance_dt,
  file.path(assess_dir, "multiclass_rf_variable_importance.csv")
)

oob_accuracy <- 1 - tail(
  multiclass_rf$err.rate[, "OOB"],
  1
)

cat("Final OOB accuracy:", oob_accuracy, "\n")


# 6. Independent multiclass test assessment ===================================

x_test <- as.data.frame(test_data[, ..selected_vars])
y_test <- factor(
  test_data$zoneID,
  levels = as.character(zones)
)

test_pred <- predict(
  multiclass_rf,
  x_test,
  type = "response"
)

test_prob <- predict(
  multiclass_rf,
  x_test,
  type = "prob"
)

test_result <- assess_multiclass(
  observed = y_test,
  predicted = test_pred,
  probability = test_prob,
  zone_values = zones
)

test_result$overall[, `:=`(
  oob_accuracy = oob_accuracy,
  ntree = ntree,
  mtry = mtry_value,
  samples_per_class_per_tree = samp_per_class,
  n_predictors = length(selected_vars),
  n_climate_predictors = length(selected_clm_vars),
  n_soil_predictors = length(selected_soil_vars),
  candidate_predictors = length(predictor_vars)
)]

fwrite(
  test_result$overall,
  file.path(assess_dir, "multiclass_rf_test_overall.csv")
)

fwrite(
  test_result$zone,
  file.path(assess_dir, "multiclass_rf_test_zone_metrics.csv")
)

fwrite(
  test_result$confusion,
  file.path(assess_dir, "multiclass_rf_test_confusion_long.csv")
)

test_predictions <- data.table(
  cell = test_data$cell,
  observed_zone = as.integer(as.character(y_test)),
  predicted_zone = as.integer(as.character(test_pred))
)

fwrite(
  test_predictions,
  file.path(assess_dir, "multiclass_rf_test_predictions.csv")
)

cat("\n[INDEPENDENT TEST COMPLETE]\n")
print(test_result$overall)


# 7. Predict the normal-period multiclass map ==================================

soil_aligned_file <- file.path(
  tmp_dir,
  "soil_predictors_aligned_to_normal_climate.tif"
)

if (
  compareGeom(
    soil_stack,
    clm_stack[[1]],
    stopOnError = FALSE
  )
) {
  soil_aligned <- soil_stack
  
} else {
  use_existing_aligned <- FALSE
  
  if (file.exists(soil_aligned_file)) {
    soil_check <- rast(soil_aligned_file)
    
    use_existing_aligned <- (
      nlyr(soil_check) == length(soil_model_vars) &&
        compareGeom(
          soil_check,
          clm_stack[[1]],
          stopOnError = FALSE
        )
    )
    
    if (use_existing_aligned) {
      names(soil_check) <- soil_model_vars
      soil_aligned <- soil_check
    }
  }
  
  if (!use_existing_aligned) {
    if (exists("soil_check")) {
      rm(soil_check)
      gc()
    }
    
    cat("\n[RESAMPLE SOIL TO CLIMATE GRID]\n")
    
    remove_raster_files(soil_aligned_file)
    
    soil_aligned <- resample(
      soil_stack,
      clm_stack[[1]],
      method = "bilinear",
      filename = soil_aligned_file,
      overwrite = TRUE,
      wopt = list(
        datatype = "FLT4S",
        gdal = "COMPRESS=LZW"
      )
    )
    
    names(soil_aligned) <- soil_model_vars
  }
}

normal_predictors_all <- c(
  clm_stack,
  soil_aligned
)

names(normal_predictors_all) <- predictor_vars

normal_predictors <- normal_predictors_all[[selected_vars]]

if (!identical(names(normal_predictors), multiclass_rf$varlist)) {
  stop("Prediction-layer names do not match the selected RF predictors.")
}

normal_climate_grid_file <- file.path(
  tmp_dir,
  "multiclass_normal_on_climate_grid.tif"
)

remove_raster_files(normal_climate_grid_file)

cat("\n[PREDICT NORMAL MULTICLASS MAP]\n")

pred_climate_grid <- terra::predict(
  normal_predictors,
  multiclass_rf,
  fun = rf_class,
  na.rm = TRUE,
  cores = predict_cores,
  cpkgs = "randomForest",
  filename = normal_climate_grid_file,
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2S",
    gdal = "COMPRESS=LZW"
  )
)

names(pred_climate_grid) <- "multiclass_zone"

normal_on_original_grid_file <- file.path(
  tmp_dir,
  "multiclass_normal_on_original_grid.tif"
)

remove_raster_files(normal_on_original_grid_file)

if (
  compareGeom(
    pred_climate_grid,
    r,
    stopOnError = FALSE
  )
) {
  pred_original_grid <- pred_climate_grid
  
} else {
  pred_original_grid <- resample(
    pred_climate_grid,
    r,
    method = "near",
    filename = normal_on_original_grid_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "INT2S",
      gdal = "COMPRESS=LZW"
    )
  )
}

ori <- subst(
  r,
  from = zones,
  to = zones,
  others = NA
)
names(ori) <- "ori"

remove_raster_files(map_file)

assigned_multiclass <- mask(
  pred_original_grid,
  ori,
  filename = map_file,
  overwrite = TRUE,
  wopt = list(
    datatype = "INT2S",
    gdal = "COMPRESS=LZW"
  )
)

assigned_freq <- freq(assigned_multiclass)
bad_value <- setdiff(assigned_freq$value, zones)

if (length(bad_value) > 0) {
  stop(
    "Multiclass map contains unexpected zone values: ",
    paste(bad_value, collapse = ", ")
  )
}

cat("[SAVED]", map_file, "\n")


# 8. Compare with the existing overlay maps ====================================

existing_threshold <- 0.2
existing_tie_tol <- 1e-4

existing_methods <- c(
  "optimized_mf",
  "optimized_rf",
  "plain_mf",
  "plain_rf"
)

existing_files <- file.path(
  base_dir,
  "result maps",
  existing_methods,
  paste0(
    "assigned_zone_normal",
    "_threshold", existing_threshold,
    "_tol", existing_tie_tol,
    "_novel99_maskNA8_noNovelNormal.tif"
  )
)

names(existing_files) <- existing_methods

map_files <- c(
  existing_files[file.exists(existing_files)],
  multiclass_rf = map_file
)

palette_file <- file.path(base_dir, "color_palette_China.csv")
category_lut <- NULL
zone_colors <- NULL

if (file.exists(palette_file)) {
  palette <- fread(palette_file)
  
  if (all(c("zoneID", "category2") %in% names(palette))) {
    category_lut <- setNames(
      as.character(palette$category2),
      as.character(palette$zoneID)
    )
  }
  
  if (all(c("zoneID", "COLOR") %in% names(palette))) {
    zone_colors <- palette$COLOR[
      match(zones, palette$zoneID)
    ]
    
    if (anyNA(zone_colors)) {
      zone_colors <- NULL
    }
  }
}

comparison_list <- list()

for (method in names(map_files)) {
  cat("[COMPARE MAP]", method, "\n")
  
  p <- rast(map_files[[method]])
  
  if (!compareGeom(p, ori, stopOnError = FALSE)) {
    cat("[SKIP GEOMETRY]", method, "\n")
    next
  }
  
  comparison_list[[method]] <- assess_map(
    original = ori,
    predicted = p,
    method = method,
    zone_values = zones,
    category_lut = category_lut
  )
}

if (!length(comparison_list)) {
  stop("No aligned maps were available for comparison.")
}

comparison_overall <- rbindlist(
  lapply(comparison_list, `[[`, "overall"),
  fill = TRUE
)

comparison_zone <- rbindlist(
  lapply(comparison_list, `[[`, "zone"),
  fill = TRUE
)

comparison_ct <- rbindlist(
  lapply(comparison_list, `[[`, "confusion"),
  fill = TRUE
)

setorder(comparison_overall, -exact_zone_accuracy)
setorder(comparison_zone, method, zone)
setorder(comparison_ct, method, ori, pred)

fwrite(
  comparison_overall,
  file.path(
    assess_dir,
    "multiclass_vs_overlay_map_overall_metrics.csv"
  )
)

fwrite(
  comparison_zone,
  file.path(
    assess_dir,
    "multiclass_vs_overlay_map_zone_metrics.csv"
  )
)

fwrite(
  comparison_ct,
  file.path(
    assess_dir,
    "multiclass_vs_overlay_map_confusion_long.csv"
  )
)

# Confusion matrix for every available method.
for (m in unique(comparison_ct$method)) {
  mat <- dcast(
    comparison_ct[method == m],
    ori ~ pred,
    value.var = "n",
    fill = 0
  )
  
  fwrite(
    mat,
    file.path(
      assess_dir,
      paste0(
        "multiclass_comparison_confusion_matrix_",
        m,
        ".csv"
      )
    )
  )
}

# Where pixels from each original zone were assigned.
error_out <- comparison_ct[
  ori != pred,
  .(pixels = sum(n)),
  by = .(method, ori, pred)
]

ori_total <- comparison_ct[
  ,
  .(original_pixels = sum(n)),
  by = .(method, ori)
]

error_out <- merge(
  error_out,
  ori_total,
  by = c("method", "ori")
)

error_out[
  ,
  percent_of_original := 100 * pixels / original_pixels
]

setorder(error_out, method, ori, -pixels)

setnames(
  error_out,
  c("ori", "pred"),
  c("original_zone", "assigned_zone")
)

fwrite(
  error_out,
  file.path(
    assess_dir,
    "multiclass_vs_overlay_errors_from_original_zone.csv"
  )
)

# Where false-positive pixels in each assigned zone came from.
error_in <- comparison_ct[
  ori != pred,
  .(pixels = sum(n)),
  by = .(method, pred, ori)
]

pred_total <- comparison_ct[
  ,
  .(predicted_pixels = sum(n)),
  by = .(method, pred)
]

error_in <- merge(
  error_in,
  pred_total,
  by = c("method", "pred")
)

error_in[
  ,
  percent_of_predicted := 100 * pixels / predicted_pixels
]

setorder(error_in, method, pred, -pixels)

setnames(
  error_in,
  c("pred", "ori"),
  c("assigned_zone", "original_source_zone")
)

fwrite(
  error_in,
  file.path(
    assess_dir,
    "multiclass_vs_overlay_errors_into_assigned_zone.csv"
  )
)

cat("\n[MAP COMPARISON COMPLETE]\n")
print(comparison_overall)


# 9. Optional six-panel map comparison =========================================

if (!is.null(zone_colors)) {
  plot_order <- c(
    "original",
    "multiclass_rf",
    "optimized_rf",
    "plain_rf",
    "plain_mf",
    "optimized_mf"
  )
  
  plot_titles <- c(
    original = "Original",
    multiclass_rf = "Multiclass RF",
    optimized_rf = "Binary overlay: optimized RF",
    plain_rf = "Binary overlay: plain RF",
    plain_mf = "Binary overlay: plain multi-Forest",
    optimized_mf = "Binary overlay: optimized multi-Forest"
  )
  
  plot_rasters <- list(original = ori)
  
  for (method in setdiff(plot_order, "original")) {
    if (method %in% names(map_files)) {
      plot_rasters[[method]] <- rast(map_files[[method]])
    }
  }
  
  plot_indexed <- function(x, title) {
    x_index <- subst(
      x,
      from = zones,
      to = seq_along(zones),
      others = NA
    )
    
    plot(
      x_index,
      col = zone_colors,
      breaks = seq(
        0.5,
        length(zones) + 0.5,
        by = 1
      ),
      legend = FALSE,
      axes = FALSE,
      main = title
    )
  }
  
  comparison_figure <- file.path(
    assess_dir,
    "normal_map_comparison_multiclass_vs_overlay.png"
  )
  
  png(
    comparison_figure,
    width = 18,
    height = 10,
    units = "in",
    res = 250
  )
  
  par(
    mfrow = c(2, 3),
    mar = c(1, 1, 2.2, 1)
  )
  
  for (method in plot_order) {
    if (method %in% names(plot_rasters)) {
      plot_indexed(
        plot_rasters[[method]],
        plot_titles[[method]]
      )
    } else {
      plot.new()
      title(main = paste(plot_titles[[method]], "(missing)"))
    }
  }
  
  dev.off()
  
  cat("[SAVED]", comparison_figure, "\n")
}


# 10. Clean up ==================================================================

rm(
  x_train,
  y_train,
  x_test,
  y_test,
  test_prob,
  test_pred,
  train_data,
  test_data,
  clm_stack,
  soil_stack,
  soil_aligned,
  normal_predictors,
  normal_predictors_all,
  pred_climate_grid,
  pred_original_grid,
  assigned_multiclass
)

gc()

cat(
  "\nSECTION 7 COMPLETE\n",
  "Model: ", model_file, "\n",
  "Map: ", map_file, "\n",
  "Assessment: ", assess_dir, "\n",
  sep = ""
)
