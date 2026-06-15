library(terra)
library(data.table)

# ------------------------------------------------------------
# Pixel-by-pixel assessment: original map vs normal prediction
#
# Outputs:
#   1. Confusion matrix:
#        rows    = original zone
#        columns = predicted zone
#
#   2. Zone-level table:
#        TP = original == zone and prediction == zone
#        FN = original == zone and prediction != zone
#        FP = original != zone and prediction == zone
# ------------------------------------------------------------

base_dir <- "H:/Jing/ecoChina2"

r <- rast(file.path(base_dir, "raster/ecosys_ori.tif"))

p <- rast(file.path(
  base_dir,
  "result maps",
  "assigned_zone_normal_threshold0.1_tol1e-04_novel99_maskNA_noNovelNormal.tif"
))

names(r) <- "ori"
names(p) <- "pred"

if (!compareGeom(r, p, stopOnError = FALSE)) {
  stop("Original and prediction rasters do not align.")
}

# Only assess pixels inside the original vegetation map.
s <- c(r, mask(p, r))

# Confusion matrix in long format.
ct <- as.data.table(crosstab(s, long = TRUE, useNA = FALSE))
setnames(ct, c("ori", "pred", "n"))

ct[, ori := as.integer(ori)]
ct[, pred := as.integer(pred)]

# Wide matrix: rows = original zone, columns = predicted zone.
mat <- dcast(ct, ori ~ pred, value.var = "n", fill = 0)

# Zone-level TP / FN / FP.
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
  total_pixels = sum(n),
  correct_pixels = sum(n[ori == pred]),
  accuracy = sum(n[ori == pred]) / sum(n)
)]

out_dir <- file.path(base_dir, "assessment")
dir.create(out_dir, showWarnings = FALSE)

fwrite(mat, file.path(out_dir, "normal_confusion_matrix_tol1e-04.csv"))
fwrite(assess, file.path(out_dir, "normal_zone_assessment_tol1e-04.csv"))
fwrite(overall, file.path(out_dir, "normal_overall_accuracy_tol1e-04.csv"))

print(overall)
print(assess[order(recall)])