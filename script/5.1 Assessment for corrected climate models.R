# 5.1 Final assessment: Plain RF and Plain MF RF
# ==============================================================================
# Final manuscript assessment uses ONLY the two correctly specified workflows:
#   rf_var : Plain RF
#   mf_var : Plain MF RF
#
# Climate models:
#   rf/clm_rfVar_zone*.Rdata
#   rf/clm_mfVar_zone*.Rdata
# Soil models:
#   rf_soil/soil_plain_zone*.Rdata
#   rf_soil/soil_mf_zone*.Rdata
#
# Both climate and soil models use zone-specific variables selected by mcRFop.
# Old incorrectly coded climate models are NOT evaluated here. Their previous
# assessment files may be kept separately for reviewer responses.
#
# Evaluation metrics include:
#   accuracy, balanced accuracy, recall, specificity, precision, F1, TSS, AUC
#
# Normal-map assessment excludes:
#   * original Zones 8 and 51 (not modeled)
#   * predicted NA pixels
# ============================================================================== 

library(terra)
library(data.table)
library(randomForest)

rm(list = ls())
gc()


# 0. Paths and settings ==========================================================

base_dir <- "H:/Jing/ecoChina2"
result_dir <- file.path(base_dir, "results")
climate_model_dir <- file.path(base_dir, "rf")
soil_model_dir <- file.path(base_dir, "rf_soil")
result_map_root <- file.path(base_dir, "result maps")
reference_file <- file.path(base_dir, "raster/ecosys_ori.tif")
assessment_dir <- file.path(base_dir, "assessment_var")

dir.create(assessment_dir, recursive = TRUE, showWarnings = FALSE)

zoneID <- c(1:7, 9:50, 52:55)
prob_threshold <- 0.5
tie_tol <- 1e-4
map_threshold <- 0.4
base_seed <- 49L

method_order <- c("rf_var", "mf_var")

method_labels <- c(
  rf_var = "Plain RF",
  mf_var = "Plain MF RF"
)

model_set <- data.table(
  method = method_order,
  climate_prefix = c("clm_rfVar_zone", "clm_mfVar_zone"),
  climate_object = c("clm_rfVar", "clm_mfVar"),
  soil_prefix = c("soil_plain_zone", "soil_mf_zone"),
  soil_object = c("soil_plain", "soil_mf")
)


# 1. Helpers ====================================================================

div <- function(a, b) {
  ifelse(is.finite(b) & b > 0, a / b, NA_real_)
}

mean_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

auc_rank <- function(y, probability) {
  # Force floating-point arithmetic.  On large test sets, integer n1 * n0
  # can exceed R's 32-bit integer range and silently turn the AUC into NA.
  n1 <- as.numeric(sum(y == 1))
  n0 <- as.numeric(sum(y == 0))
  
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }
  
  (
    sum(rank(probability, ties.method = "average")[y == 1]) -
      n1 * (n1 + 1) / 2
  ) / (n1 * n0)
}

load_rf <- function(model_file, object_name) {
  if (!file.exists(model_file)) {
    return(NULL)
  }
  
  e <- new.env()
  load(model_file, envir = e)
  
  if (!exists(object_name, envir = e)) {
    return(NULL)
  }
  
  get(object_name, envir = e)
}

get_vars <- function(model) {
  if (!is.null(model$varlist) && length(model$varlist) > 0) {
    return(model$varlist)
  }
  
  if (!is.null(model$importance) && nrow(model$importance) > 0) {
    return(rownames(model$importance))
  }
  
  NULL
}

balance_test <- function(test, zone, seed) {
  positive <- which(test$zoneID == zone)
  negative <- which(test$zoneID %in% zoneID & test$zoneID != zone)
  
  if (!length(positive) || !length(negative)) {
    return(integer())
  }
  
  set.seed(seed)
  
  if (length(negative) > length(positive)) {
    negative <- negative[
      sample.int(length(negative), length(positive))
    ]
  }
  
  c(positive, negative)
}

normal_map_file <- function(method) {
  file.path(
    result_map_root,
    method,
    paste0(
      "assigned_zone_normal_threshold",
      map_threshold,
      "_tol",
      tie_tol,
      "_novel99_maskNA8_noNovelNormal.tif"
    )
  )
}

assess_binary_model <- function(
    method_key,
    niche,
    zone,
    test,
    test_index) {
  
  config_row <- match(method_key, model_set$method)
  
  if (is.na(config_row)) {
    stop("No configuration for method: ", method_key)
  }
  
  if (niche == "climate") {
    model_file <- file.path(
      climate_model_dir,
      paste0(model_set$climate_prefix[[config_row]], zone, ".Rdata")
    )
    object_name <- model_set$climate_object[[config_row]]
  } else {
    model_file <- file.path(
      soil_model_dir,
      paste0(model_set$soil_prefix[[config_row]], zone, ".Rdata")
    )
    object_name <- model_set$soil_object[[config_row]]
  }
  
  model <- load_rf(model_file, object_name)
  
  if (is.null(model)) {
    cat("[SKIP MODEL]", method_key, "|", niche, "| zone", zone, "\n")
    return(NULL)
  }
  
  variables <- get_vars(model)
  
  if (is.null(variables) || !all(variables %in% names(test))) {
    cat("[SKIP VARIABLES]", method_key, "|", niche, "| zone", zone, "\n")
    return(NULL)
  }
  
  index <- test_index[[as.character(zone)]]
  
  if (!length(index)) {
    return(NULL)
  }
  
  x <- test[index, variables, drop = FALSE]
  y <- as.integer(test$zoneID[index] == zone)
  
  keep <- complete.cases(x)
  x <- x[keep, , drop = FALSE]
  y <- y[keep]
  
  if (!nrow(x) || length(unique(y)) < 2) {
    return(NULL)
  }
  
  probability <- predict(model, x, type = "prob")
  
  if (!("1" %in% colnames(probability))) {
    return(NULL)
  }
  
  probability <- as.numeric(probability[, "1"])
  keep <- is.finite(probability)
  probability <- probability[keep]
  y <- y[keep]
  
  prediction <- as.integer(probability >= prob_threshold)
  
  TP <- sum(y == 1 & prediction == 1)
  TN <- sum(y == 0 & prediction == 0)
  FP <- sum(y == 0 & prediction == 1)
  FN <- sum(y == 1 & prediction == 0)
  
  recall <- div(TP, TP + FN)
  specificity <- div(TN, TN + FP)
  precision <- div(TP, TP + FP)
  balanced_accuracy <- div(recall + specificity, 2)
  
  data.table(
    method = method_key,
    method_label = method_labels[[method_key]],
    niche = niche,
    zone = zone,
    probability_threshold = prob_threshold,
    sampling = "all presence + equal modeled-zone absence",
    n_test = length(y),
    presence = sum(y == 1),
    absence = sum(y == 0),
    n_variables = length(variables),
    variables = paste(variables, collapse = ","),
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    accuracy = div(TP + TN, length(y)),
    balanced_accuracy = balanced_accuracy,
    recall = recall,
    specificity = specificity,
    precision = precision,
    f1 = div(2 * precision * recall, precision + recall),
    tss = recall + specificity - 1,
    auc = auc_rank(y, probability)
  )
}


# 2. Independent climate + soil test assessment =================================

climate_test_file <- file.path(result_dir, "test_data.csv")
soil_test_file <- file.path(result_dir, "soil_test_data.csv")

if (!file.exists(climate_test_file)) {
  stop("Missing climate test data: ", climate_test_file)
}

if (!file.exists(soil_test_file)) {
  stop("Missing soil test data: ", soil_test_file)
}

climate_test <- as.data.frame(fread(climate_test_file))
soil_test <- as.data.frame(fread(soil_test_file))

climate_test$zoneID <- as.numeric(as.character(climate_test$zoneID))
soil_test$zoneID <- as.numeric(as.character(soil_test$zoneID))

climate_test_index <- setNames(
  lapply(
    zoneID,
    function(zone) balance_test(climate_test, zone, base_seed + zone)
  ),
  zoneID
)

soil_test_index <- setNames(
  lapply(
    zoneID,
    function(zone) balance_test(soil_test, zone, base_seed + 1000L + zone)
  ),
  zoneID
)

model_results <- list()

for (method_key in method_order) {
  for (zone in zoneID) {
    out <- assess_binary_model(
      method_key,
      "climate",
      zone,
      climate_test,
      climate_test_index
    )
    
    if (!is.null(out)) {
      model_results[[length(model_results) + 1L]] <- out
    }
    
    out <- assess_binary_model(
      method_key,
      "soil",
      zone,
      soil_test,
      soil_test_index
    )
    
    if (!is.null(out)) {
      model_results[[length(model_results) + 1L]] <- out
    }
  }
}

if (!length(model_results)) {
  stop("No final RF models could be assessed.")
}

model_zone_metrics <- rbindlist(model_results, fill = TRUE)

model_summary <- model_zone_metrics[
  ,
  .(
    zones_assessed = .N,
    zones_with_auc = sum(is.finite(auc)),
    mean_n_variables = mean_na(n_variables),
    mean_accuracy = mean_na(accuracy),
    mean_balanced_accuracy = mean_na(balanced_accuracy),
    mean_recall = mean_na(recall),
    mean_specificity = mean_na(specificity),
    mean_precision = mean_na(precision),
    mean_f1 = mean_na(f1),
    mean_tss = mean_na(tss),
    mean_auc = mean_na(auc)
  ),
  by = .(
    niche,
    method,
    method_label
  )
]

# The final workflow must contain one finite AUC for every modeled ecotype,
# model and niche.  This catches missing model files and any future regression
# of the integer-overflow bug before manuscript tables are regenerated.
assessment_audit <- model_zone_metrics[
  ,
  .(
    zones_assessed = .N,
    zones_with_finite_auc = sum(is.finite(auc)),
    expected_zones = length(zoneID)
  ),
  by = .(niche, method, method_label)
]

fwrite(
  assessment_audit,
  file.path(assessment_dir, "model_test_completeness_audit_var.csv")
)

incomplete <- assessment_audit[
  zones_assessed != expected_zones |
    zones_with_finite_auc != expected_zones
]

if (nrow(incomplete) > 0L) {
  print(incomplete)
  stop(
    "Final model assessment is incomplete: every niche-method combination ",
    "must have one finite AUC for each of the ",
    length(zoneID),
    " modeled zones."
  )
}

setorder(model_zone_metrics, niche, method, zone)
setorder(model_summary, niche, method)

fwrite(
  model_zone_metrics,
  file.path(assessment_dir, "model_test_zone_metrics_var.csv")
)

fwrite(
  model_summary,
  file.path(assessment_dir, "model_test_summary_var.csv")
)

fwrite(
  model_zone_metrics[niche == "climate"],
  file.path(assessment_dir, "climate_test_zone_metrics_var.csv")
)

fwrite(
  model_zone_metrics[niche == "soil"],
  file.path(assessment_dir, "soil_test_zone_metrics_var.csv")
)

fwrite(
  model_summary[niche == "climate"],
  file.path(assessment_dir, "climate_test_model_summary_var.csv")
)

fwrite(
  model_summary[niche == "soil"],
  file.path(assessment_dir, "soil_test_model_summary_var.csv")
)

cat("\n[FINAL CLIMATE + SOIL TEST ASSESSMENT COMPLETE]\n")
print(model_summary)


# 3. Normal-period assigned-map assessment ======================================

reference_map <- rast(reference_file)

original_modeled <- subst(
  reference_map,
  from = zoneID,
  to = zoneID,
  others = NA
)

names(original_modeled) <- "original_zone"

valid_original_pixels <- as.numeric(
  global(!is.na(original_modeled), "sum", na.rm = TRUE)[1, 1]
)

map_confusion_results <- list()
map_zone_results <- list()
map_overall_results <- list()

for (method_key in method_order) {
  
  prediction_file <- normal_map_file(method_key)
  
  if (!file.exists(prediction_file)) {
    cat("[SKIP NORMAL MAP]", method_key, "|", prediction_file, "\n")
    next
  }
  
  predicted_map <- rast(prediction_file)
  
  if (!compareGeom(original_modeled, predicted_map, stopOnError = FALSE)) {
    stop("Geometry mismatch for normal map: ", method_key)
  }
  
  predicted_modeled <- subst(
    predicted_map,
    from = zoneID,
    to = zoneID,
    others = NA
  )
  
  names(predicted_modeled) <- "predicted_zone"
  
  common_valid <- !is.na(original_modeled) & !is.na(predicted_modeled)
  
  original_common <- ifel(common_valid, original_modeled, NA)
  predicted_common <- ifel(common_valid, predicted_modeled, NA)
  
  names(original_common) <- "original_zone"
  names(predicted_common) <- "predicted_zone"
  
  compared_pixels <- as.numeric(
    global(common_valid, "sum", na.rm = TRUE)[1, 1]
  )
  
  confusion <- as.data.table(
    crosstab(
      c(original_common, predicted_common),
      long = TRUE,
      useNA = FALSE
    )
  )
  
  setnames(
    confusion,
    names(confusion),
    c("original_zone", "predicted_zone", "n")
  )
  
  confusion[
    ,
    `:=`(
      original_zone = as.integer(original_zone),
      predicted_zone = as.integer(predicted_zone),
      n = as.numeric(n),
      method = method_key,
      method_label = method_labels[[method_key]]
    )
  ]
  
  confusion <- confusion[
    original_zone %in% zoneID &
      predicted_zone %in% zoneID &
      n > 0
  ]
  
  total <- sum(confusion$n)
  
  zone_metrics <- rbindlist(
    lapply(
      zoneID,
      function(zone_value) {
        TP <- confusion[
          original_zone == zone_value & predicted_zone == zone_value,
          sum(n)
        ]
        
        FN <- confusion[
          original_zone == zone_value & predicted_zone != zone_value,
          sum(n)
        ]
        
        FP <- confusion[
          original_zone != zone_value & predicted_zone == zone_value,
          sum(n)
        ]
        
        TN <- total - TP - FN - FP
        
        recall <- div(TP, TP + FN)
        specificity <- div(TN, TN + FP)
        precision <- div(TP, TP + FP)
        
        data.table(
          method = method_key,
          method_label = method_labels[[method_key]],
          zone = zone_value,
          original_pixels = TP + FN,
          predicted_pixels = TP + FP,
          TP = TP,
          TN = TN,
          FP = FP,
          FN = FN,
          accuracy = div(TP + TN, total),
          balanced_accuracy = div(recall + specificity, 2),
          recall = recall,
          specificity = specificity,
          precision = precision,
          f1 = div(2 * precision * recall, precision + recall),
          tss = recall + specificity - 1
        )
      }
    ),
    fill = TRUE
  )
  
  exact_accuracy <- confusion[
    original_zone == predicted_zone,
    sum(n)
  ] / total
  
  overall <- data.table(
    method = method_key,
    method_label = method_labels[[method_key]],
    valid_original_pixels = valid_original_pixels,
    compared_pixels = compared_pixels,
    missing_predictions = valid_original_pixels - compared_pixels,
    coverage = div(compared_pixels, valid_original_pixels),
    exact_accuracy = exact_accuracy,
    macro_balanced_accuracy = mean_na(zone_metrics$balanced_accuracy),
    macro_recall = mean_na(zone_metrics$recall),
    macro_specificity = mean_na(zone_metrics$specificity),
    macro_precision = mean_na(zone_metrics$precision),
    macro_f1 = mean_na(zone_metrics$f1),
    macro_tss = mean_na(zone_metrics$tss),
    source_map = prediction_file
  )
  
  map_confusion_results[[method_key]] <- confusion
  map_zone_results[[method_key]] <- zone_metrics
  map_overall_results[[method_key]] <- overall
  
  confusion_matrix <- dcast(
    confusion,
    original_zone ~ predicted_zone,
    value.var = "n",
    fill = 0
  )
  
  fwrite(
    confusion_matrix,
    file.path(
      assessment_dir,
      paste0("normal_map_confusion_matrix_", method_key, ".csv")
    )
  )
}

if (!length(map_confusion_results)) {
  stop("No final normal maps could be assessed.")
}

normal_map_confusion <- rbindlist(map_confusion_results, fill = TRUE)
normal_map_zone_metrics <- rbindlist(map_zone_results, fill = TRUE)
normal_map_overall_metrics <- rbindlist(map_overall_results, fill = TRUE)

setorder(normal_map_confusion, method, original_zone, predicted_zone)
setorder(normal_map_zone_metrics, method, zone)
setorder(normal_map_overall_metrics, method)

fwrite(
  normal_map_confusion,
  file.path(assessment_dir, "normal_map_confusion_long_var.csv")
)

fwrite(
  normal_map_zone_metrics,
  file.path(assessment_dir, "normal_map_zone_metrics_var.csv")
)

fwrite(
  normal_map_overall_metrics,
  file.path(assessment_dir, "normal_map_overall_metrics_var.csv")
)

errors <- normal_map_confusion[original_zone != predicted_zone]

errors_from <- errors[
  ,
  .(error_pixels = sum(n)),
  by = .(method, method_label, original_zone)
][order(method, -error_pixels)]

errors_into <- errors[
  ,
  .(error_pixels = sum(n)),
  by = .(method, method_label, predicted_zone)
][order(method, -error_pixels)]

fwrite(
  errors_from,
  file.path(assessment_dir, "normal_map_errors_from_original_zone_var.csv")
)

fwrite(
  errors_into,
  file.path(assessment_dir, "normal_map_errors_into_assigned_zone_var.csv")
)

cat("\n[NORMAL MAP ASSESSMENT COMPLETE]\n")
print(normal_map_overall_metrics)

cat("\nCOMPLETE\n")
