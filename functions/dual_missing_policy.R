# Shared input-availability policy. No fitting or raster writes on source().
eco_missing_policy <- "complete_climate_soil_candidates_v1"

eco_input_signature <- function(files) {
  if (!length(files) || !all(file.exists(files))) return(NA_character_)
  z <- file.info(files)
  paste(paste(normalizePath(files, winslash = "/", mustWork = TRUE), z$size,
    sprintf("%.6f", as.numeric(z$mtime)), sep = "|"), collapse = "\n")
}

eco_policy_current <- function(file, inputs, kind) {
  marker <- paste0(file, ".missing_policy.rds")
  if (!file.exists(file) || !file.exists(marker) || !all(file.exists(inputs))) return(FALSE)
  tryCatch({
    z <- readRDS(marker)
    identical(z$policy, eco_missing_policy) && identical(z$kind, kind) &&
      identical(z$inputs, eco_input_signature(inputs)) &&
      identical(z$output_size, unname(file.info(file)$size)) &&
      identical(z$output_mtime, as.numeric(file.info(file)$mtime))
  }, error = function(e) FALSE)
}

eco_record_policy <- function(file, inputs, kind) {
  if (is.na(eco_input_signature(inputs))) stop("Cannot record policy with missing inputs: ", file)
  saveRDS(list(policy = eco_missing_policy, kind = kind,
    inputs = eco_input_signature(inputs), output_size = unname(file.info(file)$size),
    output_mtime = as.numeric(file.info(file)$mtime)), paste0(file, ".missing_policy.rds"))
}

eco_strict_dual <- function(climate, soil, gate = 0.20) {
  if (length(gate) != 1L || !is.finite(gate) || gate < 0 || gate > 1) stop("Invalid soil gate.")
  if (inherits(climate, "SpatRaster")) {
    if (gate == 0) return(terra::ifel(is.finite(climate), climate, NA_real_))
    return(terra::ifel(!is.finite(climate) | !is.finite(soil), NA_real_,
                        terra::ifel(soil > gate, climate, 0)))
  }
  answer <- climate
  if (gate == 0) {
    answer[!is.finite(climate)] <- NA_real_
  } else {
    if (length(climate) != length(soil)) stop("Climate/soil input sizes differ.")
    available <- is.finite(climate) & is.finite(soil)
    answer[!available] <- NA_real_
    answer[available & soil <= gate] <- 0
  }
  answer
}

eco_root <- function() {
  p <- Sys.getenv("ECOCHINA2_DIR", unset = getwd())
  p <- normalizePath(p, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, "script", "1. data.R"))) return(p)
    q <- dirname(p)
    if (q == p) stop("Run inside ecoChina2 or set ECOCHINA2_DIR.")
    p <- q
  }
}

eco_commit_raster <- function(temp_file, target, backup_root, project_root) {
  # Only completed temporary files are committed. Original files and sidecars
  # are moved within the same project volume, retaining a reversible backup.
  rel <- substring(normalizePath(target, winslash = "/", mustWork = FALSE),
                   nchar(normalizePath(project_root, winslash = "/")) + 2L)
  backup <- file.path(backup_root, rel)
  dir.create(dirname(backup), recursive = TRUE, showWarnings = FALSE)
  sidecars <- function(f) unique(c(f, paste0(f, c(".aux.xml", ".ovr", ".msk", ".missing_policy.rds")),
    paste0(tools::file_path_sans_ext(f), c(".aux.xml", ".ovr", ".tfw"))))
  old <- sidecars(target); old <- old[file.exists(old)]
  moved <- character()
  for (f in old) {
    suffix <- substring(f, nchar(target) + 1L)
    b <- if (startsWith(f, target)) paste0(backup, suffix) else file.path(dirname(backup), basename(f))
    if (file.exists(b)) b <- paste0(b, ".previous_", format(Sys.time(), "%Y%m%d_%H%M%OS6"))
    if (!file.rename(f, b)) {
      for (j in seq_along(moved)) file.rename(moved[j], names(moved)[j])
      stop("Cannot back up a locked output: ", f)
    }
    moved <- c(moved, stats::setNames(b, f))
  }
  if (!file.rename(temp_file, target)) {
    for (j in seq_along(moved)) file.rename(moved[j], names(moved)[j])
    stop("Cannot commit new output; prior files restored: ", target)
  }
  for (suffix in c(".aux.xml", ".ovr", ".msk")) {
    if (file.exists(paste0(temp_file, suffix))) {
      if (!file.rename(paste0(temp_file, suffix), paste0(target, suffix))) {
        stop("Output was committed but its new sidecar could not be moved: ", target, suffix)
      }
    }
  }
  invisible(backup)
}

eco_surface_summary <- function(surface, area, cutoff, inputs, cache_file) {
  # One numeric pass per unique surface, reused across species sharing a zone.
  signature <- paste(eco_input_signature(inputs), format(cutoff, digits = 17), eco_missing_policy)
  if (file.exists(cache_file)) {
    old <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(old) && identical(old$signature, signature)) return(old$value)
  }
  x <- c(surface, area)
  readStart <- terra::readStart; readValues <- terra::readValues; readStop <- terra::readStop
  readStart(x)
  on.exit(readStop(x), add = TRUE)
  n <- total_area <- suitable_area <- weighted <- total_score <- 0
  maximum <- -Inf
  rows <- max(1L, min(256L, floor(6e6 / (terra::ncol(x) * 2L))))
  for (r in seq.int(1L, terra::nrow(x), by = rows)) {
    v <- readValues(x, row = r, nrows = min(rows, terra::nrow(x) - r + 1L), mat = TRUE)
    good <- is.finite(v[, 1]) & is.finite(v[, 2])
    score <- v[good, 1]; a <- v[good, 2]
    n <- n + length(score); total_area <- total_area + sum(a)
    suitable_area <- suitable_area + sum(a[score >= cutoff])
    weighted <- weighted + sum(score * a); total_score <- total_score + sum(score)
    if (length(score)) maximum <- max(maximum, max(score))
  }
  value <- list(valid_pixels = n, valid_area_km2 = total_area, suitable_area_km2 = suitable_area,
    suitability_weighted_area_km2 = weighted, mean_dual_suitability = if (n) total_score / n else NA_real_,
    max_dual_suitability = if (is.finite(maximum)) maximum else NA_real_)
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(signature = signature, value = value), cache_file)
  value
}

eco_reuse_positive_population_ranks <- function(root, method, scenario, zones,
    population_stack, rank_files, rank_index_file, inputs, reference) {
  # Only zeros were removed: positive rank IDs and scores are unchanged.
  # Verify source backups and the original index before retaining them.
  report <- file.path(root, "assessment", "missing_climate_repair_v1")
  proof_file <- file.path(report, "dual_repair_summary.csv")
  if (!file.exists(proof_file) || !file.exists(rank_index_file) || !all(file.exists(rank_files))) return(FALSE)
  proof <- data.table::fread(proof_file)
  wanted_method <- method; wanted_scenario <- scenario; wanted_zones <- zones
  p <- proof[method == wanted_method & scenario == wanted_scenario & zone %in% wanted_zones]
  if (nrow(p) != length(zones) || any(p$added_finite_cells != 0 | p$altered_finite_values != 0 |
                                    p$removed_nonzero_cells != 0)) return(FALSE)
  index <- data.table::fread(rank_index_file)
  if (!all(c("source_zone", "rank_min_suitability") %in% names(index)) ||
      !identical(as.integer(index$source_zone), as.integer(zones)) || any(index$rank_min_suitability != 0)) return(FALSE)
  old_inputs <- file.path(report, "backups", "dual suit", method, scenario,
                         paste0("dual_suitability_zone", zones, ".tif"))
  if (!all(file.exists(old_inputs)) || min(file.info(rank_files)$mtime) < max(file.info(old_inputs)$mtime)) return(FALSE)
  if (!"old_dual_mtime" %in% names(p) ||
      as.numeric(min(file.info(rank_files)$mtime)) < max(p$old_dual_mtime)) return(FALSE)
  # Source/current correspondence must still match the repair's completion record.
  completed <- file.path(report, "repair_completed.rds")
  if (!file.exists(completed)) return(FALSE)
  state <- readRDS(completed)
  if (!identical(state$policy, eco_missing_policy)) return(FALSE)
  checked_inputs <- inputs[inputs %in% state$input_files]
  if (length(checked_inputs) < length(zones) + 1L) return(FALSE)
  current_lines <- strsplit(eco_input_signature(checked_inputs), "\n", fixed = TRUE)[[1]]
  recorded_lines <- strsplit(state$signature, "\n", fixed = TRUE)[[1]]
  if (!all(current_lines %in% recorded_lines)) return(FALSE)
  any_available <- terra::app(population_stack, fun = function(x) {
    if (is.matrix(x)) rowSums(is.finite(x)) > 0L else any(is.finite(x))
  })
  mask <- terra::ifel(any_available, 1, NA_real_)
  x <- terra::rast(rank_files[3])
  temp <- tempfile("rank_summary_", tmpdir = dirname(rank_files[3]), fileext = ".tif")
  terra::writeRaster(terra::mask(x, mask), temp, overwrite = FALSE,
                     datatype = "FLT4S", gdal = "COMPRESS=LZW")
  rm(x, mask, any_available); gc()
  eco_commit_raster(temp, rank_files[3], file.path(report, "backups"), root)
  for (f in rank_files[1:2]) eco_record_policy(f, inputs, "positive_population_rank_unchanged")
  eco_record_policy(rank_files[3], inputs, "positive_population_rank_summary")
  TRUE
}

eco_recompute_sensitivity <- function(climate_files, soil_files, reference, area,
    zones, scenario, cache_file, gates = c(0,.10,.15,.20,.25,.30,.40,.50), cutoff = .40) {
  inputs <- c(climate_files, soil_files, terra::sources(reference)[1])
  kind <- paste0("sensitivity_complete_", paste(gates, collapse = "_"), "_cutoff_", cutoff)
  if (eco_policy_current(cache_file, inputs, kind)) return(data.table::fread(cache_file))
  nz <- length(zones)
  input <- c(terra::rast(climate_files), terra::rast(soil_files), reference, area)
  result <- data.table::data.table(soil_gate = rep(gates, 2L),
    domain = rep(c("native", "common_all_gates"), each = length(gates)),
    valid_pixels = 0, valid_area_km2 = 0, below_0.4_pixels = 0,
    below_0.4_area_km2 = 0, sum_max_dual = 0, sum_area_max_dual = 0)
  rows <- max(1L, min(32L, floor(14e6 / (terra::ncol(input) * terra::nlyr(input)))))
  starts <- seq.int(1L, terra::nrow(input), by = rows)
  terra::readStart(input)
  on.exit(terra::readStop(input), add = TRUE)
  for (bi in seq_along(starts)) {
    start <- starts[bi]
    v <- terra::readValues(input, row = start, nrows = min(rows, terra::nrow(input) - start + 1L), mat = TRUE)
    # Restrict BEFORE the gate arithmetic, unlike the original full-grid audit.
    eligible <- v[, 2L*nz+1L] %in% zones & is.finite(v[, 2L*nz+2L])
    if (any(eligible)) {
      cv <- v[eligible, seq_len(nz), drop = FALSE]
      sv <- v[eligible, nz + seq_len(nz), drop = FALSE]
      a <- v[eligible, 2L*nz+2L]
      climate_complete <- rowSums(is.finite(cv)) == nz
      both_complete <- climate_complete & rowSums(is.finite(sv)) == nz
      for (g in seq_along(gates)) {
        candidate <- eco_strict_dual(cv, sv, gates[g])
        candidate[!is.finite(candidate)] <- -Inf
        mc <- max.col(candidate, ties.method = "first")
        score <- candidate[cbind(seq_len(nrow(candidate)), mc)]
        native <- if (gates[g] == 0) climate_complete else both_complete
        for (di in 1:2) {
          selected <- if (di == 1) native else both_complete
          below <- selected & score < cutoff
          ii <- g + (di - 1L)*length(gates)
          result[ii, `:=`(valid_pixels = valid_pixels + sum(selected),
            valid_area_km2 = valid_area_km2 + sum(a[selected]),
            below_0.4_pixels = below_0.4_pixels + sum(below),
            below_0.4_area_km2 = below_0.4_area_km2 + sum(a[below]),
            sum_max_dual = sum_max_dual + sum(score[selected]),
            sum_area_max_dual = sum_area_max_dual + sum(a[selected]*score[selected]))]
        }
      }
    }
    if (bi %% 10L == 0L || bi == length(starts)) {
      cat(sprintf("  [SENSITIVITY] %s %.0f%%\n", scenario, 100*bi/length(starts))); flush.console()
    }
  }
  result[, `:=`(scenario = scenario, novel_rule_applied = scenario != "normal",
    below_0.4_pixel_share = below_0.4_pixels/valid_pixels,
    below_0.4_area_share = below_0.4_area_km2/valid_area_km2,
    mean_max_dual = sum_max_dual/valid_pixels,
    area_weighted_mean_max_dual = sum_area_max_dual/valid_area_km2)]
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(result, cache_file)
  eco_record_policy(cache_file, inputs, kind)
  result
}
