# 4.1 Assign zones for rf_var and mf_var
# Vectorised block-matrix version with the original assignment rules
# ==============================================================================
#
# Inputs:
#   dual suit/rf_var/{scenario}/dual_suitability_zone*.tif
#   dual suit/mf_var/{scenario}/dual_suitability_zone*.tif
#
# Outputs:
#   result maps/rf_var/
#   result maps/mf_var/
#
# Assignment rules (identical to the original cell-wise script):
#   1. Cells outside the original map, and original Zone 8, remain NA.
#   2. If all zone suitability values are NA, the result is NA.
#   3. Normal period:
#        - no novel Zone 99;
#        - assign the zone with the maximum dual suitability.
#   4. Future periods:
#        - if maximum dual suitability < novel_threshold, assign novel Zone 99;
#        - otherwise assign the zone with the maximum dual suitability.
#   5. Zones within tie_tol of the maximum are treated as tied.
#   6. If the original zone is among the tied zones, retain the original zone.
#   7. Otherwise, randomly select one of the tied zones with equal probability.
#
# EXECUTION NOTES (behaviour differs from the app()-based draft, rules do not):
#   - Blocks are read and written explicitly with readValues / writeValues.
#     This does not rely on whether terra's app() passes a whole-block matrix or
#     one cell at a time; the matrix handed to assign_zone_block is cut here, so
#     the vectorised path is guaranteed and cannot silently fall back to a
#     per-cell call.
#   - Tie-breaking replays the original RNG stream exactly. Each cell that would
#     have reached sample.int() in the original is drawn for, one at a time, in
#     ascending cell order (the row-major order in which the original traversed
#     cells). sample.int(1, 1) is called even for a single tied zone, because
#     the original consumed a draw there too; skipping it would shift the whole
#     stream. Block boundaries do not affect the result.
#   - The seed is reset at the start of every job, so each map depends only on
#     its own inputs and not on how many draws earlier jobs consumed. Without
#     this, a normal-period map could change depending on the novel_threshold
#     used for the future jobs that ran before it, even though the normal rule
#     never uses the threshold.
#
# Existing plain_rf, plain_mf, optimized_rf and optimized_mf outputs are untouched.
# Re-running skips valid rf_var/mf_var maps and regenerates only missing or invalid outputs.
# ==============================================================================

library(terra)

rm(list = ls())
gc()


# 0. Paths and parameters ========================================================

base_dir <- "H:/Jing/ecoChina2"

dual_dir <- file.path(
  base_dir,
  "dual suit"
)

output_root <- file.path(
  base_dir,
  "result maps"
)

reference_file <- file.path(
  base_dir,
  "raster/ecosys_ori.tif"
)

zoneID <- c(
  1:7,
  9:30,
  31:50,
  52:55
)

# This threshold is used only to assign future novel Zone 99.
novel_threshold <- 0.4

novel_value <- 99
tie_tol <- 1e-4

# Retain the original random tie-breaking rule. The seed is re-applied per job
# (see the job loop below), not only here.
base_seed <- 49

# Resume mode: skip complete and valid assigned-zone maps.
reuse_existing_outputs <- TRUE

method_order <- c(
  "rf_var",
  "mf_var"
)

future_order <- c(
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

scenario_order <- c(
  "normal",
  future_order
)


# 1. Helpers ====================================================================

remove_raster_files <- function(
    filepath,
    retries = 5L,
    wait = 0.5) {
  
  base <- tools::file_path_sans_ext(
    filepath
  )
  
  files <- unique(c(
    filepath,
    paste0(filepath, ".aux.xml"),
    paste0(filepath, ".ovr"),
    paste0(filepath, ".msk"),
    paste0(base, ".aux.xml"),
    paste0(base, ".ovr"),
    paste0(base, ".tfw")
  ))
  
  files <- files[
    file.exists(files)
  ]
  
  if (length(files) == 0) {
    return(invisible(TRUE))
  }
  
  remaining <- files
  
  for (attempt in seq_len(retries)) {
    
    suppressWarnings(
      unlink(
        remaining,
        force = TRUE
      )
    )
    
    remaining <- files[
      file.exists(files)
    ]
    
    if (length(remaining) == 0) {
      return(invisible(TRUE))
    }
    
    if (attempt < retries) {
      gc()
      Sys.sleep(wait)
    }
  }
  
  stop(
    "Cannot remove existing raster file. It may still be open or locked:\n",
    paste(
      remaining,
      collapse = "\n"
    ),
    call. = FALSE
  )
}


# Matrix-based assignment.
#
# m:
#   one row per raster cell;
#   first length(zoneID) columns = zone suitability;
#   final column = original zone.
#
# Returns a numeric vector with one assigned value per row.
assign_zone_block <- function(
    m,
    is_normal,
    zoneID,
    novel_threshold,
    novel_value,
    tie_tol) {
  
  if (is.null(dim(m))) {
    m <- matrix(
      m,
      nrow = 1L
    )
  }
  
  n <- nrow(m)
  nz <- length(zoneID)
  
  if (ncol(m) != nz + 1L) {
    stop(
      "assign_zone_block expected ",
      nz + 1L,
      " columns but received ",
      ncol(m),
      "."
    )
  }
  
  suit <- m[
    ,
    seq_len(nz),
    drop = FALSE
  ]
  
  original_zone <- m[
    ,
    nz + 1L
  ]
  
  assigned <- rep(
    NA_real_,
    n
  )
  
  # Rule 1: cells outside the original map, and original Zone 8, remain NA.
  valid <- (
    !is.na(original_zone) &
      original_zone != 8
  )
  
  n_draw <- 0L
  
  if (!any(valid)) {
    attr(assigned, "n_draw") <- n_draw
    return(assigned)
  }
  
  # Replace NA suitability with -Inf so it cannot become a winning zone.
  # Rule 2: rows that were entirely NA then have a non-finite maximum.
  suit[
    is.na(suit)
  ] <- -Inf
  
  first_max_column <- max.col(
    suit,
    ties.method = "first"
  )
  
  max_suitability <- suit[
    cbind(
      seq_len(n),
      first_max_column
    )
  ]
  
  valid <- (
    valid &
      is.finite(max_suitability)
  )
  
  if (!any(valid)) {
    attr(assigned, "n_draw") <- n_draw
    return(assigned)
  }
  
  # Rule 4: future-only novel-zone assignment.
  if (!is_normal) {
    
    novel <- (
      valid &
        max_suitability < novel_threshold
    )
    
    assigned[novel] <- novel_value
    
    valid <- (
      valid &
        !novel
    )
    
    if (!any(valid)) {
      attr(assigned, "n_draw") <- n_draw
      return(assigned)
    }
  }
  
  # Rule 5: zones within tie_tol of the row maximum are tied.
  # -Inf entries give max - (-Inf) = Inf, which fails the test, so NA
  # suitabilities are excluded automatically.
  tied <- (
    max_suitability - suit <= tie_tol
  )
  
  tied[
    is.na(tied)
  ] <- FALSE
  
  # Rule 6: retain the original zone when it is among the tied winners.
  original_column <- match(
    original_zone,
    zoneID
  )
  
  retain_original <- rep(
    FALSE,
    n
  )
  
  has_original_column <- (
    valid &
      !is.na(original_column)
  )
  
  if (any(has_original_column)) {
    
    rows <- which(
      has_original_column
    )
    
    retain_original[rows] <- tied[
      cbind(
        rows,
        original_column[rows]
      )
    ]
  }
  
  retain_rows <- which(
    valid &
      retain_original
  )
  
  if (length(retain_rows) > 0) {
    assigned[retain_rows] <- original_zone[
      retain_rows
    ]
  }
  
  # Rule 7: for the remaining valid cells, draw one tied zone at random.
  # Faithful replay of the original cell-wise code: one sample.int() per cell,
  # in ascending cell order, including the single-tie case (sample.int(1, 1)
  # still consumes a draw).
  random_rows <- which(
    valid &
      !retain_original
  )
  
  if (length(random_rows) > 0) {
    for (rr in random_rows) {
      tz <- zoneID[tied[rr, ]]
      assigned[rr] <- tz[
        sample.int(
          length(tz),
          1L
        )
      ]
    }
    n_draw <- length(random_rows)
  }
  
  attr(assigned, "n_draw") <- n_draw
  assigned
}


# Check whether an existing assigned-zone map is complete and valid.
existing_output_is_valid <- function(
    output_file,
    reference_map,
    is_normal,
    zoneID,
    novel_value) {
  
  if (!file.exists(output_file)) {
    return(FALSE)
  }
  
  tryCatch({
    existing <- rast(output_file)
    
    if (!compareGeom(existing, reference_map, stopOnError = FALSE)) {
      cat("[RERUN] wrong geometry:", output_file, "\n")
      return(FALSE)
    }
    
    existing_freq <- freq(existing)
    
    if (is.null(existing_freq) || nrow(existing_freq) == 0) {
      cat("[RERUN] empty output:", output_file, "\n")
      return(FALSE)
    }
    
    existing_values <- existing_freq$value[!is.na(existing_freq$value)]
    
    allowed_value <- if (is_normal) {
      zoneID
    } else {
      c(zoneID, novel_value)
    }
    
    bad_value <- setdiff(existing_values, allowed_value)
    
    if (length(bad_value) > 0) {
      cat(
        "[RERUN] unexpected values:",
        paste(bad_value, collapse = ", "),
        "\n"
      )
      return(FALSE)
    }
    
    if (is_normal && novel_value %in% existing_values) {
      cat("[RERUN] normal output contains Zone 99\n")
      return(FALSE)
    }
    
    TRUE
    
  }, error = function(e) {
    cat(
      "[RERUN] unreadable existing output:",
      output_file,
      "|",
      conditionMessage(e),
      "\n"
    )
    FALSE
  })
}


# 2. Build jobs only from complete scenario folders =============================

jobs <- list()

for (method in method_order) {
  for (scenario in scenario_order) {
    
    input_dir <- file.path(
      dual_dir,
      method,
      scenario
    )
    
    if (!dir.exists(input_dir)) {
      next
    }
    
    files <- file.path(
      input_dir,
      paste0(
        "dual_suitability_zone",
        zoneID,
        ".tif"
      )
    )
    
    missing_zone <- zoneID[
      !file.exists(files)
    ]
    
    if (length(missing_zone) > 0) {
      
      cat(
        "[SKIP INCOMPLETE]",
        method,
        "|",
        scenario,
        "| missing zones:",
        paste(
          missing_zone,
          collapse = ", "
        ),
        "\n"
      )
      
      next
    }
    
    jobs[[
      length(jobs) + 1L
    ]] <- data.frame(
      method = method,
      scenario = scenario
    )
  }
}

if (length(jobs) == 0) {
  stop(
    "No complete rf_var or mf_var method-scenario folders were found."
  )
}

jobs <- do.call(
  rbind,
  jobs
)

print(jobs)

reference_map <- rast(
  reference_file
)


# 3. Assign zones ================================================================

for (j in seq_len(nrow(jobs))) {
  
  method <- jobs$method[j]
  scenario <- jobs$scenario[j]
  
  cat(
    "\n==============================\n",
    "ASSIGN: ",
    method,
    " | ",
    scenario,
    "\n",
    "==============================\n",
    sep = ""
  )
  
  start_time <- Sys.time()
  
  # Reset the RNG per job so each map depends only on its own inputs and not on
  # the number of draws consumed by earlier jobs.
  set.seed(base_seed)
  
  input_dir <- file.path(
    dual_dir,
    method,
    scenario
  )
  
  output_dir <- file.path(
    output_root,
    method
  )
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  files <- file.path(
    input_dir,
    paste0(
      "dual_suitability_zone",
      zoneID,
      ".tif"
    )
  )
  
  names(files) <- zoneID
  
  dual_stack <- rast(
    files
  )
  
  names(dual_stack) <- zoneID
  
  # Script 3.21 normally writes these rasters on the reference geometry.
  # Retain this guard in case a geometry mismatch is encountered.
  if (!compareGeom(
    dual_stack,
    reference_map,
    stopOnError = FALSE
  )) {
    
    cat(
      "[RESAMPLE] dual suitability rasters to original map geometry\n"
    )
    
    dual_stack <- resample(
      dual_stack,
      reference_map,
      method = "bilinear"
    )
    
    names(dual_stack) <- zoneID
  }
  
  stack_with_reference <- c(
    dual_stack,
    reference_map
  )
  
  names(stack_with_reference)[
    nlyr(stack_with_reference)
  ] <- "original_zone"
  
  is_normal <- (
    scenario == "normal"
  )
  
  output_file <- file.path(
    output_dir,
    paste0(
      "assigned_zone_",
      scenario,
      "_threshold",
      novel_threshold,
      "_tol",
      tie_tol,
      "_novel99_maskNA8_noNovelNormal.tif"
    )
  )
  
  if (
    reuse_existing_outputs &&
    existing_output_is_valid(
      output_file = output_file,
      reference_map = reference_map,
      is_normal = is_normal,
      zoneID = zoneID,
      novel_value = novel_value
    )
  ) {
    cat("[SKIP EXISTING]", output_file, "\n")
    
    rm(
      dual_stack,
      stack_with_reference
    )
    
    gc()
    next
  }
  
  # Remove a partial or invalid file left by an interrupted run.
  remove_raster_files(
    output_file
  )
  
  # Explicit block loop. The raster is traversed once: read a block, assign,
  # write the block. Cells are visited in row-major order, matching the RNG
  # consumption order of the original cell-wise code.
  out_r <- rast(reference_map)
  names(out_r) <- "assigned_zone"
  
  blocks <- writeStart(
    out_r,
    output_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "INT2S",
      gdal = "COMPRESS=LZW"
    )
  )
  
  readStart(stack_with_reference)
  
  n_draw_total <- 0L
  
  for (i in seq_len(blocks$n)) {
    
    v <- readValues(
      stack_with_reference,
      row = blocks$row[i],
      nrows = blocks$nrows[i],
      mat = TRUE
    )
    
    assigned_block <- assign_zone_block(
      v,
      is_normal = is_normal,
      zoneID = zoneID,
      novel_threshold = novel_threshold,
      novel_value = novel_value,
      tie_tol = tie_tol
    )
    
    n_draw_total <- n_draw_total +
      attr(assigned_block, "n_draw")
    
    writeValues(
      out_r,
      as.numeric(assigned_block),
      blocks$row[i],
      blocks$nrows[i]
    )
    
    cat(
      "  block",
      i,
      "/",
      blocks$n,
      "\n"
    )
  }
  
  readStop(stack_with_reference)
  writeStop(out_r)
  
  
  # 4. Validate the saved file ===================================================
  
  allowed_value <- if (is_normal) {
    zoneID
  } else {
    c(
      zoneID,
      novel_value
    )
  }
  
  saved_freq <- freq(
    rast(output_file)
  )
  
  saved_values <- saved_freq$value[
    !is.na(saved_freq$value)
  ]
  
  bad_value <- setdiff(
    saved_values,
    allowed_value
  )
  
  if (length(bad_value) > 0) {
    stop(
      "Saved raster contains unexpected values: ",
      paste(
        bad_value,
        collapse = ", "
      )
    )
  }
  
  # Explicitly verify that normal-period maps contain no novel Zone 99.
  if (
    is_normal &&
    novel_value %in% saved_values
  ) {
    stop(
      "Normal-period assigned map unexpectedly contains novel Zone ",
      novel_value,
      "."
    )
  }
  
  cat(
    "[SAVED]",
    output_file,
    "\n"
  )
  
  cat(
    "[TIEBREAK] cells drawn:",
    n_draw_total,
    "\n"
  )
  
  cat(
    "[TIME]",
    round(
      difftime(
        Sys.time(),
        start_time,
        units = "mins"
      ),
      2
    ),
    "min\n"
  )
  
  rm(
    dual_stack,
    stack_with_reference,
    out_r
  )
  
  gc()
}


cat(
  "\nCOMPLETE\n",
  "Methods: ",
  paste(
    method_order,
    collapse = ", "
  ),
  "\n",
  "Future novel threshold: ",
  novel_threshold,
  "\n",
  "Tie-breaking: original zone first; otherwise random among tied zones ",
  "(seed reset per job)\n",
  "Outputs: ",
  output_root,
  "\n",
  sep = ""
)
