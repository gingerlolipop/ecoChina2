library(terra)
library(data.table)

# ------------------------------------------------------------
# Pixel-by-pixel assessment: original map vs normal predictions
# for all four RF model versions.
# ------------------------------------------------------------

base_dir <- "H:/Jing/ecoChina2"
result_root <- file.path(base_dir, "result maps")
assess_root <- file.path(base_dir, "assessment")

r <- rast(file.path(base_dir, "raster/ecosys_ori.tif"))
names(r) <- "ori"

method_order <- c("optimized_mf", "optimized_rf", "plain_mf", "plain_rf")
methods <- method_order[dir.exists(file.path(result_root, method_order))]

assess_one <- function(method) {
  method_dir <- file.path(result_root, method)
  pred_files <- list.files(
    method_dir,
    pattern = "^assigned_zone_normal_.*\\.tif$",
    full.names = TRUE
  )
  
  if (length(pred_files) == 0) {
    warning("No normal assigned map found for method: ", method)
    return(NULL)
  }
  
  pred_file <- pred_files[1]
  p <- rast(pred_file)
  names(p) <- "pred"
  
  if (!compareGeom(r, p, stopOnError = FALSE)) {
    stop("Original and prediction rasters do not align: ", method)
  }
  
  s <- c(r, mask(p, r))
  
  ct <- as.data.table(crosstab(s, long = TRUE, useNA = FALSE))
  setnames(ct, c("ori", "pred", "n"))
  ct[, ori := as.integer(ori)]
  ct[, pred := as.integer(pred)]
  
  mat <- dcast(ct, ori ~ pred, value.var = "n", fill = 0)
  
  zones <- sort(unique(ct$ori))
  assess <- rbindlist(lapply(zones, function(z) {
    TP <- ct[ori == z & pred == z, sum(n)]
    FN <- ct[ori == z & pred != z, sum(n)]
    FP <- ct[ori != z & pred == z, sum(n)]
    
    TP <- ifelse(is.na(TP), 0, TP)
    FN <- ifelse(is.na(FN), 0, FN)
    FP <- ifelse(is.na(FP), 0, FP)
    
    recall <- TP / (TP + FN)
    precision <- TP / (TP + FP)
    f1 <- ifelse(
      is.finite(precision + recall) && precision + recall > 0,
      2 * precision * recall / (precision + recall),
      NA_real_
    )
    
    data.table(
      method = method,
      zone = z,
      original_pixels = TP + FN,
      predicted_pixels = TP + FP,
      true_positive = TP,
      false_negative = FN,
      false_positive = FP,
      recall = recall,
      precision = precision,
      f1 = f1
    )
  }))
  
  overall <- ct[, .(
    method = method,
    total_pixels = sum(n),
    correct_pixels = sum(n[ori == pred]),
    accuracy = sum(n[ori == pred]) / sum(n)
  )]
  
  out_dir <- file.path(assess_root, method)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  fwrite(mat, file.path(out_dir, paste0("normal_confusion_matrix_", method, ".csv")))
  fwrite(assess, file.path(out_dir, paste0("normal_zone_assessment_", method, ".csv")))
  fwrite(overall, file.path(out_dir, paste0("normal_overall_accuracy_", method, ".csv")))
  
  list(overall = overall, assess = assess)
}

res <- lapply(methods, assess_one)
res <- res[!vapply(res, is.null, logical(1))]

overall_all <- rbindlist(lapply(res, `[[`, "overall"), fill = TRUE)
assess_all <- rbindlist(lapply(res, `[[`, "assess"), fill = TRUE)

dir.create(assess_root, recursive = TRUE, showWarnings = FALSE)
fwrite(overall_all, file.path(assess_root, "normal_overall_accuracy_all_methods.csv"))
fwrite(assess_all, file.path(assess_root, "normal_zone_assessment_all_methods.csv"))

print(overall_all[order(-accuracy)])
print(assess_all[order(method, recall)])
