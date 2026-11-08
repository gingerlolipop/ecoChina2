# Refresh binary Multi-Forest versus existing multiclass assignment comparisons.
# Reads only saved maps; never loads, fits or projects an RF model.
local({
  library(terra)
  library(data.table)
  root <- normalizePath(Sys.getenv("ECOCHINA2_DIR", getwd()), winslash = "/", mustWork = TRUE)
  rule_file <- file.path(root, "functions", "future_topk_rule.R")
  source(rule_file, local = environment())
  table_dir <- file.path(root, "assessment", "multiclass_rf")
  cache_dir <- file.path(table_dir, "comparison_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  scenarios <- c("normal", "2011-2040SSP245", "2041-2070SSP245", "2071-2100SSP245",
                 "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585")
  first_existing <- function(paths) {
    hit <- paths[file.exists(paths)]
    if (length(hit)) hit[[1]] else NA_character_
  }
  signature <- function(files) {
    if (!length(files) || !all(file.exists(files))) return(NA_character_)
    z <- file.info(files)
    paste(paste(normalizePath(files, winslash = "/", mustWork = TRUE), z$size,
                sprintf("%.6f", as.numeric(z$mtime)), sep = "|"), collapse = "\n")
  }
  safe_divide <- function(a, b) {
    ifelse(is.finite(b) & b > 0, a / b, NA_real_)
  }
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


  results <- list()
  status <- list()
  for (scenario in scenarios) {
    binary_file <- first_existing(c(
      file.path(root, "result maps", "mf_var", paste0("assigned_zone_", scenario,
        "_threshold0.4_tol1e-04_novel99_maskNA8_noNovelNormal.tif")),
      file.path(root, "result maps", "mf_var", paste0("assigned_zone_", scenario, ".tif")),
      file.path(root, "maps", paste0("assigned_zone_", scenario, ".tif"))
    ))
    established <- paste0("assigned_zone_", scenario, "_multiclass_rf.tif")
    multiclass_file <- first_existing(c(
      file.path(root, "result maps", "multiclass_rf", established),
      file.path(root, "result maps", "multiclass_rf", paste0("assigned_zone_", scenario, ".tif")),
      file.path(root, "multiclass", "maps", paste0("assigned_zone_", scenario, ".tif"))
    ))
    if (is.na(binary_file) || is.na(multiclass_file)) {
      status[[scenario]] <- data.table(scenario = scenario, status = "missing_inputs",
        binary_file = binary_file, multiclass_file = multiclass_file)
      eco_topk_milestone("SKIP MULTICLASS COMPARISON", paste(scenario, "| saved map unavailable"))
      next
    }
    inputs <- c(binary_file, multiclass_file,
                file.path(root, "script", "7.1 Refresh assigned-map multiclass comparison.R"))
    before <- signature(inputs)
    cache_file <- file.path(cache_dir, paste0(scenario, ".rds"))
    cache <- if (file.exists(cache_file)) tryCatch(readRDS(cache_file), error = function(e) NULL) else NULL
    if (!is.null(cache) && identical(cache$signature, before)) {
      result <- cache$result
      eco_topk_milestone("REUSE MULTICLASS COMPARISON", scenario)
    } else {
      eco_topk_milestone("MULTICLASS COMPARISON", scenario)
      result <- compare_assigned_maps(binary_file, multiclass_file, scenario)
      if (!identical(before, signature(inputs))) stop("Map inputs changed during comparison: ", scenario)
      saveRDS(list(signature = before, result = result), cache_file)
    }
    results[[scenario]] <- result
    status[[scenario]] <- data.table(scenario = scenario, status = "current",
      binary_file = binary_file, multiclass_file = multiclass_file)
  }
  empty_summary <- data.table(scenario = character(), common_pixels = numeric(),
    current_zone_common_pixels = numeric(), binary_novel_pixels = numeric(),
    binary_novel_share = numeric(), exact_agreement_all = numeric(),
    exact_agreement_current_zones = numeric())
  empty_confusion <- data.table(binary_zone = integer(), multiclass_zone = integer(),
    pixels = numeric(), scenario = character())
  summary <- if (length(results)) rbindlist(lapply(results, `[[`, "summary")) else empty_summary
  confusion <- if (length(results)) rbindlist(lapply(results, `[[`, "confusion")) else empty_confusion
  # Missing scenarios are omitted from rebuilt tables and recorded explicitly;
  # an older comparison is never left in place as if it matched the new maps.
  fwrite(summary, file.path(table_dir, "binary_mf_vs_multiclass_common_mask.csv"))
  fwrite(confusion, file.path(table_dir, "binary_mf_vs_multiclass_confusion.csv"))
  refresh_status <- rbindlist(status)
  fwrite(refresh_status, file.path(table_dir, "comparison_refresh_status.csv"))
  eco_topk_milestone("MULTICLASS COMPARISON COMPLETE", sprintf("%d/%d scenarios current",
    sum(refresh_status$status == "current"), length(scenarios)))
})
