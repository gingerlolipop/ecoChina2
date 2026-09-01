library(data.table)

base_dir <- "H:/Jing/ecoChina2"
acc_dir <- file.path(base_dir, "accuracy_soil")
mod_dir <- file.path(base_dir, "rf_soil")
out_dir <- file.path(acc_dir, "PA_check")

zoneID <- c(1:7, 9:55)
MF_NR <- 1.3

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read class counts stored inside an existing RF object.
# This does not rerun sampling, prediction, or model training.
read_model_counts <- function(file, object_name) {
  if (!file.exists(file)) return(NULL)
  
  e <- new.env()
  nm <- load(file, envir = e)
  if (!(object_name %in% nm)) return(NULL)
  
  m <- e[[object_name]]
  
  if (!is.null(m$y)) {
    y <- as.character(m$y)
    return(data.table(
      n_presence = sum(y == "1", na.rm = TRUE),
      n_absence = sum(y == "0", na.rm = TRUE),
      source = paste0(basename(file), ": model$y")
    ))
  }
  
  if (!is.null(m$confusion)) {
    cm <- as.matrix(m$confusion)
    cm <- cm[, setdiff(colnames(cm), "class.error"), drop = FALSE]
    n <- rowSums(cm, na.rm = TRUE)
    
    if (all(c("0", "1") %in% names(n))) {
      return(data.table(
        n_presence = as.integer(n["1"]),
        n_absence = as.integer(n["0"]),
        source = paste0(basename(file), ": model$confusion")
      ))
    }
  }
  
  NULL
}

# The training script already saved the sampled-data counts here.
summary_file <- file.path(acc_dir, "soil_rf_accuracy_summary.csv")
acc <- if (file.exists(summary_file)) fread(summary_file) else data.table()

# Counts in the sampled dataset entering the four model workflows.
# Prefer the small summary CSV; use the saved plain RF only for missing zones.
base_count <- rbindlist(lapply(zoneID, function(i) {
  x <- data.table()
  
  if (nrow(acc) > 0 &&
      all(c("zone", "n_presence", "n_absence") %in% names(acc))) {
    x <- acc[zone == i]
    if ("model" %in% names(acc) && any(x$model == "plain_rf")) {
      x <- x[model == "plain_rf"]
    }
  }
  
  if (nrow(x) > 0) {
    x <- x[.N]
    return(data.table(
      zone = i,
      n_presence = as.integer(x$n_presence),
      n_absence = as.integer(x$n_absence),
      source = basename(summary_file)
    ))
  }
  
  x <- read_model_counts(
    file.path(mod_dir, paste0("soil_plain_zone", i, ".Rdata")),
    "soil_plain"
  )
  
  if (!is.null(x)) return(cbind(data.table(zone = i), x))
  
  data.table(
    zone = i,
    n_presence = NA_integer_,
    n_absence = NA_integer_,
    source = "not found"
  )
}), fill = TRUE)

# Final class counts retained by optimized single RF after classOP().
# These may differ from the sampled input data if classOP removed observations.
opt_count <- rbindlist(lapply(zoneID, function(i) {
  x <- read_model_counts(
    file.path(mod_dir, paste0("soil_zOp_zone", i, ".Rdata")),
    "soil_zOp"
  )
  
  if (is.null(x)) {
    x <- data.table(
      n_presence = NA_integer_,
      n_absence = NA_integer_,
      source = "not recoverable from saved optimized RF"
    )
  }
  
  cbind(data.table(zone = i), x)
}), fill = TRUE)

# One row per zone and model workflow.
pa <- rbindlist(lapply(zoneID, function(i) {
  b <- base_count[zone == i]
  o <- opt_count[zone == i]
  
  p <- b$n_presence
  a <- b$n_absence
  
  # Each multi-Forest subforest uses all presences and at most 1.3x absences.
  mf_a <- if (is.na(p) || is.na(a)) NA_integer_ else min(a, floor(p * MF_NR))
  
  rbindlist(list(
    data.table(
      zone = i,
      model = "plain_rf",
      stage = "actual training data",
      n_presence = p,
      n_absence = a,
      source = b$source
    ),
    data.table(
      zone = i,
      model = "optimized_rf",
      stage = "final model after classOP",
      n_presence = o$n_presence,
      n_absence = o$n_absence,
      source = o$source
    ),
    data.table(
      zone = i,
      model = "plain_mf",
      stage = "input to each subforest",
      n_presence = p,
      n_absence = mf_a,
      source = paste0("saved base counts + MF_NR=", MF_NR)
    ),
    data.table(
      zone = i,
      model = "optimized_mf",
      stage = "each subforest before classOP",
      n_presence = p,
      n_absence = mf_a,
      source = paste0("saved base counts + MF_NR=", MF_NR)
    )
  ))
}), fill = TRUE)

pa[, `:=`(
  absence_per_presence = n_absence / n_presence,
  presence_percent = 100 * n_presence / (n_presence + n_absence),
  balanced = !is.na(n_presence) & n_presence > 0 & n_absence <= 1.5 * n_presence
)]

setorder(pa, zone, model)

fwrite(pa, file.path(out_dir, "soil_RF_presence_absence_ratio.csv"))
fwrite(
  pa[is.na(n_presence) | is.na(n_absence) | !balanced],
  file.path(out_dir, "soil_RF_presence_absence_problem.csv")
)

cat("\nPresence/absence ratios by soil RF model:\n")
print(pa)

cat("\nUnbalanced or missing:\n")
print(pa[is.na(n_presence) | is.na(n_absence) | !balanced])

cat("\nSaved to:\n", out_dir, "\n")
