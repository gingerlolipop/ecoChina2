# 5.1 Assessment for corrected climate models
# ==============================================================================
# Compares the two original plain climate models with the two corrected models:
#
#   plain_rf  : original ordinary RF using all climate variables
#   rf_var    : corrected ordinary RF using mcRFop variables
#   plain_mf  : original plain multi-Forest using all climate variables
#   mf_var    : corrected plain multi-Forest using mcRFop variables
#
# Assessment levels:
#   1. Independent balanced climate test-set assessment.
#   2. Normal-period assigned-map reconstruction assessment.
#
# Soil models are not reassessed because rf_var and mf_var reuse the unchanged
# plain_rf and plain_mf soil suitability rasters.
#
# New outputs are written to:
#   H:/Jing/ecoChina2/assessment_var
#
# Existing assessment files are not overwritten.
# ==============================================================================

library(terra)
library(data.table)
library(randomForest)

rm(list = ls())
gc()


# 0. Paths and settings ==========================================================

base_dir <- "H:/Jing/ecoChina2"

result_dir <- file.path(
  base_dir,
  "results"
)

model_dir <- file.path(
  base_dir,
  "rf"
)

result_map_root <- file.path(
  base_dir,
  "result maps"
)

reference_file <- file.path(
  base_dir,
  "raster/ecosys_ori.tif"
)

assessment_dir <- file.path(
  base_dir,
  "assessment_var"
)

dir.create(
  assessment_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

zoneID <- c(
  1:7,
  9:50,
  52:55
)

prob_threshold <- 0.5
tie_tol <- 1e-4
base_seed <- 49L

method_set <- data.table(
  method = c(
    "plain_rf",
    "rf_var",
    "plain_mf",
    "mf_var"
  ),
  model_family = c(
    "single RF",
    "single RF",
    "multi-Forest",
    "multi-Forest"
  ),
  model_version = c(
    "original",
    "corrected",
    "original",
    "corrected"
  ),
  clm_prefix = c(
    "clm_plain_zone",
    "clm_rfVar_zone",
    "clm_mf_zone",
    "clm_mfVar_zone"
  ),
  clm_object = c(
    "clm_plain",
    "clm_rfVar",
    "clm_mf",
    "clm_mfVar"
  ),
  map_threshold = c(
    0.2,
    0.4,
    0.2,
    0.4
  )
)

method_order <- method_set$method


# 1. Helpers ====================================================================

div <- function(a, b) {
  ifelse(
    is.finite(b) & b > 0,
    a / b,
    NA_real_
  )
}


mean_na <- function(x) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    mean(
      x,
      na.rm = TRUE
    )
  }
}


auc_rank <- function(y, probability) {
  
  n1 <- as.numeric(
    sum(y == 1)
  )
  
  n0 <- as.numeric(
    sum(y == 0)
  )
  
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }
  
  (
    sum(
      rank(
        probability,
        ties.method = "average"
      )[y == 1]
    ) -
      n1 * (n1 + 1) / 2
  ) / (
    n1 * n0
  )
}


load_rf <- function(
    model_file,
    object_name) {
  
  if (!file.exists(model_file)) {
    return(NULL)
  }
  
  e <- new.env()
  
  load(
    model_file,
    envir = e
  )
  
  if (!exists(
    object_name,
    envir = e
  )) {
    return(NULL)
  }
  
  get(
    object_name,
    envir = e
  )
}


get_vars <- function(model) {
  
  if (
    !is.null(model$varlist) &&
    length(model$varlist) > 0
  ) {
    return(model$varlist)
  }
  
  if (
    !is.null(model$importance) &&
    nrow(model$importance) > 0
  ) {
    return(
      rownames(model$importance)
    )
  }
  
  NULL
}


balance_test <- function(
    test,
    zone,
    seed) {
  
  positive <- which(
    test$zoneID == zone
  )
  
  negative <- which(
    test$zoneID %in% zoneID &
      test$zoneID != zone
  )
  
  if (
    !length(positive) ||
    !length(negative)
  ) {
    return(integer())
  }
  
  set.seed(seed)
  
  if (length(negative) > length(positive)) {
    
    negative <- negative[
      sample.int(
        length(negative),
        length(positive)
      )
    ]
  }
  
  c(
    positive,
    negative
  )
}


normal_map_file <- function(
    method,
    threshold) {
  
  file.path(
    result_map_root,
    method,
    paste0(
      "assigned_zone_normal",
      "_threshold",
      threshold,
      "_tol",
      tie_tol,
      "_novel99_maskNA8_noNovelNormal.tif"
    )
  )
}


# 2. Independent climate test-set assessment ====================================

climate_test_file <- file.path(
  result_dir,
  "test_data.csv"
)

if (!file.exists(climate_test_file)) {
  stop(
    "Missing climate test data: ",
    climate_test_file
  )
}

climate_test <- as.data.frame(
  fread(climate_test_file)
)

climate_test$zoneID <- as.numeric(
  as.character(
    climate_test$zoneID
  )
)

# The same balanced observations are used for all four climate methods.
test_index <- setNames(
  lapply(
    zoneID,
    function(zone) {
      balance_test(
        climate_test,
        zone,
        base_seed + zone
      )
    }
  ),
  zoneID
)


assess_climate_model <- function(
    method_key,
    zone) {
  
  # Use match() rather than a data.table expression such as
  # method_set$method == method. Inside data.table, the argument name
  # "method" can be masked by the column with the same name, selecting
  # every row and producing a model_file vector of length > 1.
  config_row <- match(
    method_key,
    method_set$method
  )
  
  if (is.na(config_row)) {
    stop(
      "No model configuration found for method: ",
      method_key
    )
  }
  
  if (sum(method_set$method == method_key) != 1L) {
    stop(
      "Method configuration is not unique for: ",
      method_key
    )
  }
  
  config <- method_set[
    config_row
  ]
  
  model_file <- file.path(
    model_dir,
    paste0(
      config$clm_prefix[[1]],
      zone,
      ".Rdata"
    )
  )
  
  object_name <- config$clm_object[[1]]
  
  model <- load_rf(
    model_file,
    object_name
  )
  
  if (is.null(model)) {
    
    cat(
      "[SKIP MODEL]",
      method_key,
      "| zone",
      zone,
      "\n"
    )
    
    return(NULL)
  }
  
  variables <- get_vars(model)
  
  if (
    is.null(variables) ||
    !all(
      variables %in% names(climate_test)
    )
  ) {
    
    cat(
      "[SKIP VARIABLES]",
      method_key,
      "| zone",
      zone,
      "\n"
    )
    
    return(NULL)
  }
  
  index <- test_index[[as.character(zone)]]
  
  if (!length(index)) {
    
    cat(
      "[SKIP TEST]",
      method_key,
      "| zone",
      zone,
      "\n"
    )
    
    return(NULL)
  }
  
  x <- climate_test[
    index,
    variables,
    drop = FALSE
  ]
  
  y <- as.integer(
    climate_test$zoneID[index] == zone
  )
  
  keep <- complete.cases(x)
  
  x <- x[
    keep,
    ,
    drop = FALSE
  ]
  
  y <- y[keep]
  
  if (
    !nrow(x) ||
    length(unique(y)) < 2
  ) {
    return(NULL)
  }
  
  probability <- predict(
    model,
    x,
    type = "prob"
  )
  
  if (!(
    "1" %in% colnames(probability)
  )) {
    
    cat(
      "[SKIP PROBABILITY]",
      method_key,
      "| zone",
      zone,
      "\n"
    )
    
    return(NULL)
  }
  
  probability <- as.numeric(
    probability[
      ,
      "1"
    ]
  )
  
  keep <- is.finite(probability)
  
  probability <- probability[keep]
  y <- y[keep]
  
  prediction <- as.integer(
    probability >= prob_threshold
  )
  
  TP <- sum(
    y == 1 &
      prediction == 1
  )
  
  TN <- sum(
    y == 0 &
      prediction == 0
  )
  
  FP <- sum(
    y == 0 &
      prediction == 1
  )
  
  FN <- sum(
    y == 1 &
      prediction == 0
  )
  
  recall <- div(
    TP,
    TP + FN
  )
  
  specificity <- div(
    TN,
    TN + FP
  )
  
  precision <- div(
    TP,
    TP + FP
  )
  
  balanced_accuracy <- div(
    recall + specificity,
    2
  )
  
  data.table(
    method = method_key,
    model_family = config$model_family[[1]],
    model_version = config$model_version[[1]],
    zone = zone,
    probability_threshold = prob_threshold,
    sampling = "all presence + equal modeled-zone absence",
    n_test = length(y),
    presence = sum(y == 1),
    absence = sum(y == 0),
    n_variables = length(variables),
    variables = paste(
      variables,
      collapse = ","
    ),
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    accuracy = div(
      TP + TN,
      length(y)
    ),
    balanced_accuracy = balanced_accuracy,
    recall = recall,
    specificity = specificity,
    precision = precision,
    f1 = div(
      2 * precision * recall,
      precision + recall
    ),
    tss = recall + specificity - 1,
    auc = auc_rank(
      y,
      probability
    )
  )
}


climate_zone_results <- list()

for (method in method_order) {
  for (zone in zoneID) {
    
    result <- assess_climate_model(
      method,
      zone
    )
    
    if (!is.null(result)) {
      climate_zone_results[[length(climate_zone_results) + 1L]] <- result
    }
  }
}

if (!length(climate_zone_results)) {
  stop(
    "No climate models could be assessed."
  )
}

climate_zone_metrics <- rbindlist(
  climate_zone_results,
  fill = TRUE
)

climate_model_summary <- climate_zone_metrics[
  ,
  .(
    zones_assessed = .N,
    mean_n_variables = mean_na(n_variables),
    mean_accuracy = mean_na(accuracy),
    mean_balanced_accuracy = mean_na(
      balanced_accuracy
    ),
    mean_recall = mean_na(recall),
    mean_specificity = mean_na(specificity),
    mean_precision = mean_na(precision),
    mean_f1 = mean_na(f1),
    mean_tss = mean_na(tss),
    mean_auc = mean_na(auc)
  ),
  by = .(
    method,
    model_family,
    model_version
  )
]

setorder(
  climate_model_summary,
  model_family,
  model_version
)

fwrite(
  climate_zone_metrics,
  file.path(
    assessment_dir,
    "climate_test_zone_metrics_var.csv"
  )
)

fwrite(
  climate_model_summary,
  file.path(
    assessment_dir,
    "climate_test_model_summary_var.csv"
  )
)

cat(
  "\n[CLIMATE TEST ASSESSMENT COMPLETE]\n"
)

print(
  climate_model_summary[
    order(
      model_family,
      -mean_auc
    )
  ]
)


# 3. Normal-period map reconstruction assessment ===============================

reference_map <- rast(
  reference_file
)

# Exclude unmodeled original Zones 8 and 51.
original_modeled <- subst(
  reference_map,
  from = zoneID,
  to = zoneID,
  others = NA
)

names(original_modeled) <- "original_zone"

valid_original_pixels <- global(
  !is.na(original_modeled),
  fun = "sum",
  na.rm = TRUE
)[1, 1]

map_confusion_results <- list()
map_zone_results <- list()
map_overall_results <- list()

for (row_index in seq_len(nrow(method_set))) {
  
  method_key <- method_set$method[[row_index]]
  
  prediction_file <- normal_map_file(
    method_key,
    method_set$map_threshold[[row_index]]
  )
  
  if (!file.exists(prediction_file)) {
    
    cat(
      "[SKIP NORMAL MAP]",
      method_key,
      "| missing:",
      prediction_file,
      "\n"
    )
    
    next
  }
  
  cat(
    "[ASSESS NORMAL MAP]",
    method_key,
    "\n"
  )
  
  predicted_map <- rast(
    prediction_file
  )
  
  if (!compareGeom(
    original_modeled,
    predicted_map,
    stopOnError = FALSE
  )) {
    
    cat(
      "[SKIP GEOMETRY]",
      method_key,
      "\n"
    )
    
    next
  }
  
  predicted_modeled <- subst(
    predicted_map,
    from = zoneID,
    to = zoneID,
    others = NA
  )
  
  names(predicted_modeled) <- "predicted_zone"
  
  compared_pixels <- global(
    !is.na(original_modeled) &
      !is.na(predicted_modeled),
    fun = "sum",
    na.rm = TRUE
  )[1, 1]
  
  confusion <- as.data.table(
    crosstab(
      c(
        original_modeled,
        predicted_modeled
      ),
      long = TRUE,
      useNA = FALSE
    )
  )
  
  setnames(
    confusion,
    c(
      "original_zone",
      "predicted_zone",
      "n"
    )
  )
  
  confusion[
    ,
    `:=`(
      method = method_key,
      model_family = method_set$model_family[[row_index]],
      model_version = method_set$model_version[[row_index]],
      original_zone = as.integer(
        original_zone
      ),
      predicted_zone = as.integer(
        predicted_zone
      ),
      n = as.numeric(n)
    )
  ]
  
  total_compared <- sum(
    confusion$n,
    na.rm = TRUE
  )
  
  correct <- confusion[
    original_zone == predicted_zone,
    sum(
      n,
      na.rm = TRUE
    )
  ]
  
  map_overall_results[[method_key]] <- data.table(
    method = method_key,
    model_family = method_set$model_family[[row_index]],
    model_version = method_set$model_version[[row_index]],
    map_threshold_in_filename =
      method_set$map_threshold[[row_index]],
    valid_original_pixels =
      valid_original_pixels,
    compared_pixels =
      compared_pixels,
    missing_predictions =
      valid_original_pixels -
      compared_pixels,
    coverage = div(
      compared_pixels,
      valid_original_pixels
    ),
    exact_accuracy = div(
      correct,
      compared_pixels
    ),
    prediction_file =
      prediction_file
  )
  
  map_zone_results[[method_key]] <- rbindlist(
    lapply(
      zoneID,
      function(zone) {
        
        TP <- confusion[
          original_zone == zone &
            predicted_zone == zone,
          sum(
            n,
            na.rm = TRUE
          )
        ]
        
        FN <- confusion[
          original_zone == zone &
            predicted_zone != zone,
          sum(
            n,
            na.rm = TRUE
          )
        ]
        
        FP <- confusion[
          original_zone != zone &
            predicted_zone == zone,
          sum(
            n,
            na.rm = TRUE
          )
        ]
        
        TN <- total_compared -
          TP -
          FN -
          FP
        
        recall <- div(
          TP,
          TP + FN
        )
        
        specificity <- div(
          TN,
          TN + FP
        )
        
        precision <- div(
          TP,
          TP + FP
        )
        
        data.table(
          method = method_key,
          model_family =
            method_set$model_family[[row_index]],
          model_version =
            method_set$model_version[[row_index]],
          zone = zone,
          original_pixels = TP + FN,
          predicted_pixels = TP + FP,
          TP = TP,
          TN = TN,
          FP = FP,
          FN = FN,
          accuracy = div(
            TP + TN,
            total_compared
          ),
          balanced_accuracy = div(
            recall + specificity,
            2
          ),
          recall = recall,
          specificity = specificity,
          precision = precision,
          f1 = div(
            2 * precision * recall,
            precision + recall
          ),
          tss = recall +
            specificity -
            1
        )
      }
    ),
    fill = TRUE
  )
  
  map_confusion_results[[method_key]] <- confusion
}

if (!length(map_overall_results)) {
  stop(
    "No completed aligned normal maps could be assessed."
  )
}

normal_map_confusion <- rbindlist(
  map_confusion_results,
  fill = TRUE
)

normal_map_zone_metrics <- rbindlist(
  map_zone_results,
  fill = TRUE
)

normal_map_overall_metrics <- rbindlist(
  map_overall_results,
  fill = TRUE
)

fwrite(
  normal_map_confusion,
  file.path(
    assessment_dir,
    "normal_map_confusion_long_var.csv"
  )
)

fwrite(
  normal_map_zone_metrics,
  file.path(
    assessment_dir,
    "normal_map_zone_metrics_var.csv"
  )
)

fwrite(
  normal_map_overall_metrics,
  file.path(
    assessment_dir,
    "normal_map_overall_metrics_var.csv"
  )
)


# Confusion matrices: rows = original zones; columns = predicted zones.
for (method_key in unique(
  normal_map_confusion[["method"]]
)) {
  
  # Explicit vector comparison avoids data.table's .. pronoun and avoids
  # masking the loop variable with the column named "method".
  confusion_method <- normal_map_confusion[
    normal_map_confusion[["method"]] == method_key
  ]
  
  confusion_matrix <- dcast(
    confusion_method,
    original_zone ~ predicted_zone,
    value.var = "n",
    fill = 0
  )
  
  fwrite(
    confusion_matrix,
    file.path(
      assessment_dir,
      paste0(
        "normal_map_confusion_matrix_",
        method_key,
        "_var.csv"
      )
    )
  )
}


# Error directions from original zones.
errors_from_original <- normal_map_confusion[
  original_zone != predicted_zone,
  .(
    pixels = sum(
      n,
      na.rm = TRUE
    )
  ),
  by = .(
    method,
    model_family,
    model_version,
    original_zone,
    predicted_zone
  )
]

original_totals <- normal_map_confusion[
  ,
  .(
    original_total = sum(
      n,
      na.rm = TRUE
    )
  ),
  by = .(
    method,
    original_zone
  )
]

errors_from_original <- merge(
  errors_from_original,
  original_totals,
  by = c(
    "method",
    "original_zone"
  ),
  all.x = TRUE
)

errors_from_original[
  ,
  share_of_original := div(
    pixels,
    original_total
  )
]

setorder(
  errors_from_original,
  method,
  original_zone,
  -pixels
)

fwrite(
  errors_from_original,
  file.path(
    assessment_dir,
    "normal_map_errors_from_original_zone_var.csv"
  )
)


# Error directions into assigned zones.
errors_into_assigned <- normal_map_confusion[
  original_zone != predicted_zone,
  .(
    pixels = sum(
      n,
      na.rm = TRUE
    )
  ),
  by = .(
    method,
    model_family,
    model_version,
    original_zone,
    predicted_zone
  )
]

predicted_totals <- normal_map_confusion[
  ,
  .(
    predicted_total = sum(
      n,
      na.rm = TRUE
    )
  ),
  by = .(
    method,
    predicted_zone
  )
]

errors_into_assigned <- merge(
  errors_into_assigned,
  predicted_totals,
  by = c(
    "method",
    "predicted_zone"
  ),
  all.x = TRUE
)

errors_into_assigned[
  ,
  share_of_prediction := div(
    pixels,
    predicted_total
  )
]

setorder(
  errors_into_assigned,
  method,
  predicted_zone,
  -pixels
)

fwrite(
  errors_into_assigned,
  file.path(
    assessment_dir,
    "normal_map_errors_into_assigned_zone_var.csv"
  )
)

cat(
  "\n[NORMAL MAP ASSESSMENT COMPLETE]\n"
)

print(
  normal_map_overall_metrics[
    order(
      model_family,
      -exact_accuracy
    )
  ]
)


# 4. Corrected-minus-original comparison ========================================

comparison_pairs <- data.table(
  model_family = c(
    "single RF",
    "multi-Forest"
  ),
  original_method = c(
    "plain_rf",
    "plain_mf"
  ),
  corrected_method = c(
    "rf_var",
    "mf_var"
  )
)

comparison_results <- list()

for (pair_index in seq_len(
  nrow(comparison_pairs)
)) {
  
  family <- comparison_pairs$model_family[
    pair_index
  ]
  
  original_method <- comparison_pairs$original_method[
    pair_index
  ]
  
  corrected_method <- comparison_pairs$corrected_method[
    pair_index
  ]
  
  climate_original <- climate_model_summary[
    method == original_method
  ]
  
  climate_corrected <- climate_model_summary[
    method == corrected_method
  ]
  
  map_original <- normal_map_overall_metrics[
    method == original_method
  ]
  
  map_corrected <- normal_map_overall_metrics[
    method == corrected_method
  ]
  
  comparison_results[[length(comparison_results) + 1L]] <- data.table(
    model_family = family,
    original_method = original_method,
    corrected_method = corrected_method,
    delta_mean_n_variables =
      climate_corrected$mean_n_variables -
      climate_original$mean_n_variables,
    delta_mean_balanced_accuracy =
      climate_corrected$mean_balanced_accuracy -
      climate_original$mean_balanced_accuracy,
    delta_mean_f1 =
      climate_corrected$mean_f1 -
      climate_original$mean_f1,
    delta_mean_auc =
      climate_corrected$mean_auc -
      climate_original$mean_auc,
    delta_normal_exact_accuracy =
      map_corrected$exact_accuracy -
      map_original$exact_accuracy,
    delta_normal_coverage =
      map_corrected$coverage -
      map_original$coverage
  )
}

comparison_delta <- rbindlist(
  comparison_results,
  fill = TRUE
)

fwrite(
  comparison_delta,
  file.path(
    assessment_dir,
    "corrected_minus_original_summary_var.csv"
  )
)

cat(
  "\n[CORRECTED MINUS ORIGINAL]\n"
)

print(comparison_delta)

cat(
  "\nCOMPLETE\n",
  "Assessment directory: ",
  assessment_dir,
  "\n",
  "Independent assessment: climate test set only\n",
  "Map assessment: normal-period reconstruction\n",
  sep = ""
)
