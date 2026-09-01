# Assign vegetation zones from Multi-Forest dual suitability
# ==============================================================================
# Run after script 3.2.
#
# Rules:
#   - cells outside the 53 modeled reference zones remain NA;
#   - the normal-period map cannot contain novel Zone 99;
#   - a future cell is Zone 99 when every dual suitability is below 0.40;
#   - values within 1e-4 of the maximum are tied;
#   - a tied original zone is retained; other ties are resolved reproducibly.

library(terra)

rm(list = ls())
gc()


# 0. Paths and parameters =====================================================

find_project_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "script", "1. data.R"))) return(path)
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
project_dir <- if (nzchar(env_root)) {
  normalizePath(env_root, winslash = "/", mustWork = TRUE)
} else {
  find_project_root()
}

reference_file <- file.path(project_dir, "raster", "ecosys_ori.tif")
dual_root <- file.path(project_dir, "dual suit", "mf_var")
map_dir <- file.path(project_dir, "result maps", "mf_var")
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)

zoneID <- c(1:7, 9:50, 52:55)
scenario_order <- c(
  "normal",
  "2011-2040SSP245", "2041-2070SSP245", "2071-2100SSP245",
  "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585"
)

novel_threshold <- 0.40
novel_value <- 99L
tie_tol <- 1e-4
base_seed <- 49L
reuse_existing <- env_flag("ECOCHINA2_REUSE_MAPS", TRUE) &&
  !env_flag("ECOCHINA2_FORCE_MAPS", FALSE)

if (!file.exists(reference_file)) {
  stop("Missing reference raster: ", reference_file)
}

reference_map <- rast(reference_file)


# 1. Block assignment =========================================================

assign_zone_block <- function(values, is_normal) {
  if (is.null(dim(values))) {
    values <- matrix(values, nrow = 1L)
  }

  n <- nrow(values)
  nz <- length(zoneID)

  if (ncol(values) != nz + 1L) {
    stop("Expected ", nz + 1L, " input columns; received ", ncol(values), ".")
  }

  suit <- values[, seq_len(nz), drop = FALSE]
  original_zone <- values[, nz + 1L]
  assigned <- rep(NA_real_, n)

  valid <- !is.na(original_zone) & original_zone %in% zoneID

  if (!any(valid)) {
    attr(assigned, "n_draw") <- 0L
    return(assigned)
  }

  # -Inf cannot win. Rows that were entirely NA remain invalid below.
  suit[is.na(suit)] <- -Inf
  first_max <- max.col(suit, ties.method = "first")
  max_suit <- suit[cbind(seq_len(n), first_max)]
  valid <- valid & is.finite(max_suit)

  if (!any(valid)) {
    attr(assigned, "n_draw") <- 0L
    return(assigned)
  }

  if (!is_normal) {
    novel <- valid & max_suit < novel_threshold
    assigned[novel] <- novel_value
    valid <- valid & !novel
  }

  if (!any(valid)) {
    attr(assigned, "n_draw") <- 0L
    return(assigned)
  }

  tied <- max_suit - suit <= tie_tol
  tied[is.na(tied)] <- FALSE

  original_col <- match(original_zone, zoneID)
  retain_original <- rep(FALSE, n)
  rows_with_original <- which(valid & !is.na(original_col))

  if (length(rows_with_original)) {
    retain_original[rows_with_original] <- tied[
      cbind(rows_with_original, original_col[rows_with_original])
    ]
  }

  retained_rows <- which(valid & retain_original)
  assigned[retained_rows] <- original_zone[retained_rows]

  # Rows are processed in raster-cell order. Resetting the seed for every map
  # makes results independent of the order in which scenarios are run.
  random_rows <- which(valid & !retain_original)

  for (row in random_rows) {
    tied_zones <- zoneID[tied[row, ]]
    assigned[row] <- tied_zones[sample.int(length(tied_zones), 1L)]
  }

  attr(assigned, "n_draw") <- length(random_rows)
  assigned
}

valid_existing_map <- function(file, is_normal, inputs = character()) {
  if (!file.exists(file) || is.na(file.info(file)$size) || file.info(file)$size <= 0) {
    return(FALSE)
  }

  tryCatch({
    x <- rast(file)
    # Metadata-only validation avoids scanning every pixel on routine reruns.
    geometry_ok <- nlyr(x) == 1L &&
      compareGeom(x, reference_map, stopOnError = FALSE)
    inputs <- inputs[file.exists(inputs)]
    current <- !length(inputs) ||
      file.info(file)$mtime >= max(file.info(inputs)$mtime)
    geometry_ok && isTRUE(current)
  }, error = function(e) FALSE)
}


# 2. Assign each scenario =====================================================

for (scenario in scenario_order) {
  is_normal <- scenario == "normal"
  input_dir <- file.path(dual_root, scenario)
  dual_files <- file.path(
    input_dir,
    paste0("dual_suitability_zone", zoneID, ".tif")
  )
  output_file <- file.path(
    map_dir,
    paste0(
      "assigned_zone_", scenario,
      "_threshold", novel_threshold,
      "_tol", tie_tol,
      "_novel99_maskNA8_noNovelNormal.tif"
    )
  )

  if (reuse_existing && valid_existing_map(
    output_file, is_normal, c(reference_file, dual_files)
  )) {
    cat("[USE EXISTING]", output_file, "\n")
    next
  }

  missing <- zoneID[!file.exists(dual_files)]

  if (length(missing)) {
    stop(
      "Incomplete dual-suitability stack for ", scenario,
      ". Missing zone(s): ", paste(missing, collapse = ", ")
    )
  }

  cat("\n[ASSIGN]", scenario, "\n")
  set.seed(base_seed)

  dual_stack <- rast(dual_files)
  names(dual_stack) <- as.character(zoneID)

  if (!compareGeom(dual_stack, reference_map, stopOnError = FALSE)) {
    dual_stack <- resample(dual_stack, reference_map, method = "bilinear")
    names(dual_stack) <- as.character(zoneID)
  }

  assignment_input <- c(dual_stack, reference_map)
  names(assignment_input)[nlyr(assignment_input)] <- "original_zone"

  output_raster <- rast(reference_map)
  names(output_raster) <- "assigned_zone"

  blocks <- writeStart(
    output_raster,
    output_file,
    overwrite = TRUE,
    wopt = list(datatype = "INT2S", gdal = "COMPRESS=LZW")
  )
  readStart(assignment_input)

  draws <- 0L

  for (block in seq_len(blocks$n)) {
    input_values <- readValues(
      assignment_input,
      row = blocks$row[block],
      nrows = blocks$nrows[block],
      mat = TRUE
    )

    assigned_values <- assign_zone_block(input_values, is_normal)
    draws <- draws + attr(assigned_values, "n_draw")

    writeValues(
      output_raster,
      assigned_values,
      blocks$row[block],
      blocks$nrows[block]
    )

    cat("  block", block, "/", blocks$n, "\n")
  }

  readStop(assignment_input)
  writeStop(output_raster)

  saved_values <- freq(rast(output_file))$value
  saved_values <- saved_values[!is.na(saved_values)]
  allowed <- if (is_normal) zoneID else c(zoneID, novel_value)
  unexpected <- setdiff(saved_values, allowed)

  if (length(unexpected)) {
    stop("Unexpected assigned value(s): ", paste(unexpected, collapse = ", "))
  }

  if (is_normal && novel_value %in% saved_values) {
    stop("The normal-period map contains novel Zone 99.")
  }

  cat("[SAVED]", output_file, "| random tie breaks:", draws, "\n")

  rm(dual_stack, assignment_input, output_raster)
  gc()
}

cat(
  "\nCOMPLETE\n",
  "Novel threshold (future only): ", novel_threshold, "\n",
  "Tie tolerance: ", tie_tol, "\n",
  "Maps: ", map_dir, "\n",
  sep = ""
)
