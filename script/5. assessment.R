library(CEMT)
library(terra)
library(data.table)
library(randomForest)

# Outputs:
#   1. rf_test_zone_metrics.csv: balanced test metrics by method, niche and zone
#   2. rf_test_model_summary.csv: mean RF metrics by method and niche
#   3. normal_map_confusion_long.csv: original-predicted pixel counts
#   4. normal_map_zone_metrics.csv: map metrics by method and zone
#   5. normal_map_overall_metrics.csv: overall map accuracy and coverage
# The last three files include only completed normal maps.

base_dir <- "H:/Jing/ecoChina2"
result_dir <- file.path(base_dir, "results")
result_root <- file.path(base_dir, "result maps")
assess_dir <- file.path(base_dir, "assessment")
dir.create(assess_dir, recursive = TRUE, showWarnings = FALSE)

zoneID <- c(1:7, 9:50, 52:55)
prob_threshold <- 0.5
map_threshold <- 0.2
tie_tol <- 1e-4
base_seed <- 49L

method_order <- c(
  "optimized_mf",
  "optimized_rf",
  "plain_mf",
  "plain_rf"
)

model_set <- data.frame(
  method = method_order,
  clm_prefix = c(
    "clm_mfOp_zone", "clm_zOp_zone",
    "clm_mf_zone", "clm_plain_zone"
  ),
  clm_object = c(
    "clim_mfOp", "clim_zOp",
    "clm_mf", "clm_plain"
  ),
  soil_prefix = c(
    "soil_mfOp_zone", "soil_zOp_zone",
    "soil_mf_zone", "soil_plain_zone"
  ),
  soil_object = c(
    "soil_mfOp", "soil_zOp",
    "soil_mf", "soil_plain"
  )
)

div <- function(a, b) {
  ifelse(is.finite(b) & b > 0, a / b, NA_real_)
}

mean_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

auc_rank <- function(y, p) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  
  if (n1 == 0 || n0 == 0) return(NA_real_)
  
  (sum(rank(p, ties.method = "average")[y == 1]) -
      n1 * (n1 + 1) / 2) / (n1 * n0)
}

load_rf <- function(file, object) {
  if (!file.exists(file)) return(NULL)
  
  e <- new.env()
  load(file, envir = e)
  
  if (!exists(object, envir = e)) return(NULL)
  get(object, envir = e)
}

get_vars <- function(m) {
  if (!is.null(m$varlist)) return(m$varlist)
  rownames(m$importance)
}

balance_test <- function(test, zone, seed) {
  pos <- which(test$zoneID == zone)
  neg <- which(test$zoneID != zone & !is.na(test$zoneID))
  
  if (!length(pos) || !length(neg)) return(integer())
  
  set.seed(seed)
  
  if (length(neg) > length(pos)) {
    neg <- neg[sample.int(length(neg), length(pos))]
  }
  
  c(pos, neg)
}


# 1. Independent climate and soil RF assessment ===============================

clm_test_file <- file.path(result_dir, "test_data.csv")
soil_test_file <- file.path(result_dir, "soil_test_data.csv")

if (!file.exists(clm_test_file)) {
  stop("Missing climate test data: ", clm_test_file)
}

if (!file.exists(soil_test_file)) {
  stop("Missing soil test data: ", soil_test_file)
}

clm_test <- as.data.frame(fread(clm_test_file))
soil_test <- as.data.frame(fread(soil_test_file))

clm_test$zoneID <- as.numeric(as.character(clm_test$zoneID))
soil_test$zoneID <- as.numeric(as.character(soil_test$zoneID))

# The same balanced test observations are used by all four models.
test_index <- list(
  climate = setNames(
    lapply(
      zoneID,
      function(z) balance_test(clm_test, z, base_seed + z)
    ),
    zoneID
  ),
  soil = setNames(
    lapply(
      zoneID,
      function(z) balance_test(soil_test, z, base_seed + 1000L + z)
    ),
    zoneID
  )
)

assess_rf <- function(method, niche, zone) {
  cfg <- model_set[model_set$method == method, , drop = FALSE]
  is_clm <- niche == "climate"
  
  prefix <- if (is_clm) cfg$clm_prefix else cfg$soil_prefix
  object <- if (is_clm) cfg$clm_object else cfg$soil_object
  test <- if (is_clm) clm_test else soil_test
  model_dir <- if (is_clm) "rf" else "rf_soil"
  
  file <- file.path(
    base_dir,
    model_dir,
    paste0(prefix, zone, ".Rdata")
  )
  
  m <- load_rf(file, object)
  
  if (is.null(m)) {
    cat("[SKIP MODEL]", method, "|", niche, "| zone", zone, "\n")
    return(NULL)
  }
  
  vars <- get_vars(m)
  
  if (is.null(vars) || !all(vars %in% names(test))) {
    cat("[SKIP VARS]", method, "|", niche, "| zone", zone, "\n")
    return(NULL)
  }
  
  idx <- test_index[[niche]][[as.character(zone)]]
  
  if (!length(idx)) {
    cat("[SKIP TEST]", method, "|", niche, "| zone", zone, "\n")
    return(NULL)
  }
  
  x <- test[idx, vars, drop = FALSE]
  y <- as.integer(test$zoneID[idx] == zone)
  
  keep <- complete.cases(x)
  x <- x[keep, , drop = FALSE]
  y <- y[keep]
  
  if (!nrow(x) || length(unique(y)) < 2) return(NULL)
  
  prob <- predict(m, x, type = "prob")
  
  if (!("1" %in% colnames(prob))) {
    cat("[SKIP PROB]", method, "|", niche, "| zone", zone, "\n")
    return(NULL)
  }
  
  prob <- as.numeric(prob[, "1"])
  keep <- is.finite(prob)
  prob <- prob[keep]
  y <- y[keep]
  
  pred <- as.integer(prob >= prob_threshold)
  
  TP <- sum(y == 1 & pred == 1)
  TN <- sum(y == 0 & pred == 0)
  FP <- sum(y == 0 & pred == 1)
  FN <- sum(y == 1 & pred == 0)
  
  recall <- div(TP, TP + FN)
  specificity <- div(TN, TN + FP)
  precision <- div(TP, TP + FP)
  balanced_accuracy <- div(recall + specificity, 2)
  tss <- recall + specificity - 1
  
  data.table(
    method,
    niche,
    zone,
    threshold = prob_threshold,
    sampling = "all presence + equal absence",
    n_test = length(y),
    presence = sum(y == 1),
    absence = sum(y == 0),
    TP, TN, FP, FN,
    accuracy = div(TP + TN, length(y)),
    balanced_accuracy,
    recall,
    specificity,
    precision,
    f1 = div(2 * precision * recall, precision + recall),
    tss,
    auc = auc_rank(y, prob)
  )
}

rf_list <- list()

for (method in method_order) {
  for (niche in c("climate", "soil")) {
    for (zone in zoneID) {
      out <- assess_rf(method, niche, zone)
      
      if (!is.null(out)) {
        rf_list[[length(rf_list) + 1L]] <- out
      }
    }
  }
}

if (!length(rf_list)) {
  stop("No RF models could be assessed.")
}

rf_test <- rbindlist(rf_list)

rf_summary <- rf_test[, .(
  zones_assessed = .N,
  zones_with_presence = sum(presence > 0),
  mean_accuracy = mean_na(accuracy),
  mean_balanced_accuracy = mean_na(balanced_accuracy),
  mean_recall = mean_na(recall),
  mean_specificity = mean_na(specificity),
  mean_precision = mean_na(precision),
  mean_f1 = mean_na(f1),
  mean_tss = mean_na(tss),
  mean_auc = mean_na(auc)
), by = .(method, niche)]

fwrite(
  rf_test,
  file.path(assess_dir, "rf_test_zone_metrics.csv")
)

fwrite(
  rf_summary,
  file.path(assess_dir, "rf_test_model_summary.csv")
)

cat("\n[RF ASSESSMENT COMPLETE]\n")
print(rf_summary[order(niche, -mean_auc)])


# 2. Completed normal-map assessment ==========================================

r <- rast(file.path(base_dir, "raster/ecosys_ori.tif"))
ori <- ifel(is.na(r) | r == 8, NA, r)
names(ori) <- "ori"

normal_files <- file.path(
  result_root,
  method_order,
  paste0(
    "assigned_zone_normal",
    "_threshold", map_threshold,
    "_tol", tie_tol,
    "_novel99_maskNA8_noNovelNormal.tif"
  )
)

names(normal_files) <- method_order
normal_files <- normal_files[file.exists(normal_files)]

if (!length(normal_files)) {
  cat("\n[SKIP MAP ASSESSMENT] No completed normal maps found.\n")
  
} else {
  map_ct <- list()
  map_zone <- list()
  map_overall <- list()
  
  for (method in names(normal_files)) {
    cat("[ASSESS MAP]", method, "\n")
    
    p <- rast(normal_files[[method]])
    
    if (!compareGeom(ori, p, stopOnError = FALSE)) {
      cat("[SKIP GEOMETRY]", method, "\n")
      next
    }
    
    pred <- ifel(is.na(ori), NA, ifel(is.na(p), -999, p))
    names(pred) <- "pred"
    
    ct <- as.data.table(
      crosstab(c(ori, pred), long = TRUE, useNA = FALSE)
    )
    
    setnames(ct, c("ori", "pred", "n"))
    
    ct[, `:=`(
      method = method,
      ori = as.integer(ori),
      pred = as.integer(pred),
      n = as.numeric(n)
    )]
    
    total <- sum(ct$n, na.rm = TRUE)
    predicted <- ct[pred != -999, sum(n, na.rm = TRUE)]
    correct <- ct[ori == pred, sum(n, na.rm = TRUE)]
    
    map_overall[[method]] <- data.table(
      method,
      valid_pixels = total,
      predicted_pixels = predicted,
      missing_predictions = total - predicted,
      coverage = div(predicted, total),
      accuracy_predicted = div(correct, predicted),
      accuracy_all = div(correct, total)
    )
    
    map_zone[[method]] <- rbindlist(
      lapply(sort(unique(ct$ori)), function(z) {
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
          TP, TN, FP, FN,
          recall,
          specificity,
          precision,
          f1 = div(2 * precision * recall, precision + recall),
          tss = recall + specificity - 1
        )
      })
    )
    
    ct[pred == -999, pred := NA_integer_]
    map_ct[[method]] <- ct
  }
  
  if (length(map_overall)) {
    map_overall <- rbindlist(map_overall)
    
    fwrite(
      rbindlist(map_ct, fill = TRUE),
      file.path(assess_dir, "normal_map_confusion_long.csv")
    )
    
    fwrite(
      rbindlist(map_zone, fill = TRUE),
      file.path(assess_dir, "normal_map_zone_metrics.csv")
    )
    
    fwrite(
      map_overall,
      file.path(assess_dir, "normal_map_overall_metrics.csv")
    )
    
    cat("\n[MAP ASSESSMENT COMPLETE]\n")
    print(map_overall[order(-accuracy_all)])
    
  } else {
    cat("\n[SKIP MAP ASSESSMENT] No aligned normal maps found.\n")
  }
}