# Audit the soil-gate sensitivity domains without refitting or modifying models.
# Run from the ecoChina2 root: source("script/audit_soil_domain.R")
# Requires terra and data.table, plus existing climate/soil suitability rasters.
# Replicates script/3.2's CSV algorithm exactly; strict finite-climate semantics
# are diagnosed separately. Never source script/3.2 from this audit.

local({
  if (!requireNamespace("terra", quietly = TRUE) ||
      !requireNamespace("data.table", quietly = TRUE)) {
    stop("Install packages terra and data.table first.")
  }
  library(terra)
  library(data.table)

  find_root <- function(path = getwd()) {
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    repeat {
      if (file.exists(file.path(path, "script", "1. data.R"))) return(path)
      parent <- dirname(path)
      if (parent == path) stop("Run inside ecoChina2 or set ECOCHINA2_DIR.")
      path <- parent
    }
  }
  task_root <- Sys.getenv("ECOCHINA2_DIR", unset = "")
  task_root <- if (nzchar(task_root)) normalizePath(task_root, winslash = "/",
                                                 mustWork = TRUE) else find_root()
  path_setting <- function(key, fallback) Sys.getenv(key, unset = fallback)
  reference_file <- path_setting("ECOCHINA2_REFERENCE_RASTER",
                                file.path(task_root, "raster", "ecosys_ori.tif"))
  climate_root <- path_setting("ECOCHINA2_AUDIT_CLIMATE_SUITABILITY_DIR",
                              file.path(task_root, "clim suitability", "mf_var"))
  soil_root <- path_setting("ECOCHINA2_AUDIT_SOIL_SUITABILITY_DIR",
                           file.path(task_root, "soil suitability", "plain_mf", "normal"))
  prior_csv <- file.path(task_root, "assessment", "soil_gate_sensitivity.csv")
  region_file <- path_setting("ECOCHINA2_REGION_RASTER", "")
  region_lookup_file <- path_setting("ECOCHINA2_REGION_LOOKUP", "")
  output_dir <- path_setting("ECOCHINA2_DOMAIN_AUDIT_OUT", file.path(
    task_root, "assessment", paste0("soil_domain_audit_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  ))
  if (dir.exists(output_dir)) {
    stop("Output directory already exists; choose a new ECOCHINA2_DOMAIN_AUDIT_OUT: ", output_dir)
  }
  zones <- c(1:7, 9:50, 52:55)
  gates <- c(0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
  primary_index <- match(0.20, gates)
  cutoff <- 0.40
  scenarios <- c("normal", "2011-2040SSP245", "2041-2070SSP245", "2071-2100SSP245",
                 "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585")
  soil_files <- file.path(soil_root, paste0("soil_suit_zone", zones, ".tif"))
  climate_files <- lapply(scenarios, function(s) file.path(
    climate_root, s, paste0("clim_suit_zone", zones, ".tif")))
  required <- c(reference_file, soil_files, unlist(climate_files),
                if (nzchar(region_file)) region_file)
  missing <- required[!file.exists(required)]
  if (length(missing)) stop("Missing input(s); this audit does not rebuild them:\n",
                            paste(missing, collapse = "\n"))
  reference <- rast(reference_file)
  if (nlyr(reference) != 1L || !nzchar(crs(reference))) {
    stop("Reference must be a single-layer raster with a defined CRS.")
  }
  checked_stack <- function(files) {
    x <- rast(files)
    if (nlyr(x) != length(files) ||
        !compareGeom(x, reference, stopOnError = FALSE)) {
      stop("Input geometry differs from the reference. This audit never resamples: ", files[1])
    }
    x
  }
  soil <- checked_stack(soil_files)
  invisible(lapply(climate_files, checked_stack))
  region <- if (nzchar(region_file)) checked_stack(region_file) else NULL
  if (!is.null(region) && nlyr(region) != 1L) stop("Region raster must have one layer.")
  region_lookup <- if (nzchar(region_lookup_file)) fread(region_lookup_file) else NULL
  if (!is.null(region_lookup) &&
      !all(c("region_id", "region_name") %in% names(region_lookup))) {
    stop("ECOCHINA2_REGION_LOOKUP must contain region_id and region_name.")
  }
  if (!is.null(region_lookup) && anyDuplicated(region_lookup$region_id)) {
    stop("Region lookup contains duplicate region_id values.")
  }
  dir.create(output_dir, recursive = TRUE)
  dir.create(file.path(output_dir, "temporary"))
  terraOptions(tempdir = file.path(output_dir, "temporary"), memfrac = 0.15)
  area_raster <- cellSize(reference, unit = "km")
  # Identical modeled-observed eligibility to script/3.2, not a prediction mask.
  group_lookup <- data.table(
    domain_code = 1:4,
    domain_group = c("common_valid", "climate_only_valid_only", "gate_020_valid_only",
                     "invalid_both")
  )
  fwrite(group_lookup, file.path(output_dir, "domain_codes.csv"))
  manifest <- data.table(path = normalizePath(required, winslash = "/"),
                         size_bytes = file.info(required)$size,
                         modified = as.character(file.info(required)$mtime))
  fwrite(manifest, file.path(output_dir, "input_manifest.csv"))
  captured_session <- capture.output(sessionInfo())
  writeLines(c("Source: audit_soil_domain.R", paste("Run UTC:", format(Sys.time(), tz = "UTC")),
               paste("Project:", task_root), "Original sensitivity semantics retained.",
               "Strict diagnostic requires finite climate AND finite soil per candidate.",
               "1-degree bins are geographic bins, not provinces or ecological regions.",
               "Below-cutoff reference cells are not assigned future novel-ecosystem labels.",
               captured_session), file.path(output_dir, "run_info.txt"))

  max_rows <- function(x) {
    x[!is.finite(x)] <- -Inf
    idx <- max.col(x, ties.method = "first")
    x[cbind(seq_len(nrow(x)), idx)]
  }
  summarize_metric <- function(selected, score, n_suitable, area, scenario,
                               gate, domain, semantics = "existing_csv") {
    n <- sum(selected)
    a <- sum(area[selected])
    below <- selected & score < cutoff
    data.table(scenario = scenario, soil_gate = gate, domain = domain,
      semantics = semantics, novel_rule_applied = scenario != "normal",
      valid_pixels = n, valid_area_km2 = a,
      below_0.4_pixels = sum(below), below_0.4_area_km2 = sum(area[below]),
      sum_max_dual = sum(score[selected]),
      sum_area_max_dual = sum(area[selected] * score[selected]),
      sum_suitable_ecotypes = sum(n_suitable[selected]),
      sum_area_suitable_ecotypes = sum(area[selected] * n_suitable[selected]))
  }
  all_metrics <- list()
  all_zone_metrics <- list()
  all_zone <- list()
  all_grid <- list()
  all_region <- list()
  all_issues <- list()
  mask_changes <- list()
  changed_file <- file.path(output_dir, "exclusive_domain_cells.csv")
  counter <- 0L
  reference_domain <- NULL

  for (s_index in seq_along(scenarios)) {
    scenario <- scenarios[s_index]
    cat("[DOMAIN AUDIT]", scenario, "\n")
    climate <- checked_stack(climate_files[[s_index]])
    input <- c(climate, soil, reference, area_raster)
    if (!is.null(region)) input <- c(input, region)
    n_zone <- length(zones)
    rows_per_block <- max(1L, floor(5e6 / (ncol(input) * nlyr(input))))
    starts <- seq.int(1L, nrow(input), by = rows_per_block)
    mask_file <- file.path(output_dir, paste0("domain_", scenario, ".tif"))
    domain_raster <- rast(reference)
    names(domain_raster) <- "domain_code"
    writeStart(domain_raster, filename = mask_file, overwrite = FALSE,
               datatype = "INT1U", gdal = "COMPRESS=LZW", NAflag = 255)
    if (!is.null(reference_domain)) readStart(reference_domain)
    readStart(input)
    for (start in starts) {
      rows <- min(rows_per_block, nrow(input) - start + 1L)
      v <- readValues(input, row = start, nrows = rows, mat = TRUE)
      cv <- v[, seq_len(n_zone), drop = FALSE]
      sv <- v[, n_zone + seq_len(n_zone), drop = FALSE]
      observed <- v[, 2L * n_zone + 1L]
      area <- v[, 2L * n_zone + 2L]
      eligible <- observed %in% zones & is.finite(area)
      scores <- matrix(NA_real_, nrow(v), length(gates))
      suitable <- matrix(0L, nrow(v), length(gates))
      strict_scores <- matrix(NA_real_, nrow(v), length(gates))
      strict_suitable <- matrix(0L, nrow(v), length(gates))
      for (g in seq_along(gates)) {
        candidate <- cv
        if (gates[g] > 0) {
          # Exact existing summary implementation: soil failure writes zero
          # even if the corresponding climate prediction was NA.
          candidate[!is.na(sv) & sv <= gates[g]] <- 0
          candidate[is.na(sv)] <- -Inf
        }
        candidate[!is.finite(candidate)] <- -Inf
        suitable[, g] <- rowSums(candidate >= cutoff)
        scores[, g] <- max_rows(candidate)
        if (gates[g] > 0) {
          candidate[!is.finite(cv) | !is.finite(sv)] <- -Inf
        }
        strict_scores[, g] <- max_rows(candidate)
        strict_suitable[, g] <- rowSums(candidate >= cutoff)
      }
      valid <- is.finite(scores) & eligible
      strict_valid <- is.finite(strict_scores) & eligible
      common <- valid[, 1L] & valid[, primary_index]
      common_all_gates <- rowSums(valid) == length(gates)
      code <- rep(NA_integer_, nrow(v))
      code[eligible] <- 4L
      code[valid[, 1L] & !valid[, primary_index]] <- 2L
      code[!valid[, 1L] & valid[, primary_index]] <- 3L
      code[common] <- 1L
      writeValues(domain_raster, code, start = start, nrows = rows)

      counter <- counter + 1L
      metrics <- list()
      zone_metrics <- list()
      for (g in seq_along(gates)) {
        for (which_domain in c("native", "common_gate0_gate020", "common_all_gates")) {
          if (which_domain == "common_gate0_gate020" &&
              !g %in% c(1L, primary_index)) next
          selected <- switch(which_domain, native = valid[, g],
                             common_gate0_gate020 = common & valid[, g],
                             common_all_gates = common_all_gates)
          metrics[[length(metrics) + 1L]] <- summarize_metric(
            selected, scores[, g], suitable[, g], area, scenario, gates[g], which_domain)
          if (any(selected)) {
            z <- data.table(reference_zone = observed[selected], area = area[selected],
                            score = scores[selected, g], ns = suitable[selected, g])
            z <- z[, .(valid_pixels = .N, valid_area_km2 = sum(area),
              below_0.4_pixels = sum(score < cutoff),
              below_0.4_area_km2 = sum(area[score < cutoff]),
              sum_max_dual = sum(score), sum_area_max_dual = sum(area * score),
              sum_suitable_ecotypes = sum(ns), sum_area_suitable_ecotypes = sum(area * ns)),
              by = reference_zone]
            z[, `:=`(scenario = scenario, soil_gate = gates[g], domain = which_domain)]
            zone_metrics[[length(zone_metrics) + 1L]] <- z
          }
        }
        metrics[[length(metrics) + 1L]] <- summarize_metric(
          strict_valid[, g], strict_scores[, g], strict_suitable[, g], area, scenario,
          gates[g], "native", "strict_finite_climate_and_soil")
        rescue <- valid[, g] & !strict_valid[, g]
        any_climate <- rowSums(is.finite(cv)) > 0L
        all_issues[[length(all_issues) + 1L]] <- data.table(
          scenario = scenario, soil_gate = gates[g],
          zero_fill_only_valid_pixels = sum(rescue),
          zero_fill_only_valid_area_km2 = sum(area[rescue]),
          zero_fill_valid_no_climate_pixels = sum(rescue & !any_climate),
          zero_fill_valid_no_climate_area_km2 = sum(area[rescue & !any_climate]),
          max_changed_under_strict_pixels = sum(eligible &
            (is.finite(scores[, g]) != is.finite(strict_scores[, g]) |
             (is.finite(scores[, g]) & is.finite(strict_scores[, g]) &
                abs(scores[, g] - strict_scores[, g]) > 1e-8))))
      }
      all_metrics[[counter]] <- rbindlist(metrics)
      all_zone_metrics[[counter]] <- rbindlist(zone_metrics)
      cell <- (start - 1) * ncol(input) + seq_len(nrow(v))
      xy <- xyFromCell(reference, cell)
      if (!is.lonlat(reference)) {
        xy <- terra::project(xy, from = crs(reference), to = "EPSG:4326")
      }
      tab <- data.table(domain_code = code, reference_zone = observed, area_km2 = area,
                        longitude = xy[, 1L], latitude = xy[, 2L])
      tab <- tab[!is.na(domain_code)]
      tab[, `:=`(longitude_bin = floor(longitude), latitude_bin = floor(latitude))]
      all_zone[[counter]] <- tab[, .(pixels = .N, area_km2 = sum(area_km2)),
                                by = .(domain_code, reference_zone)][, scenario := scenario]
      all_grid[[counter]] <- tab[, .(pixels = .N, area_km2 = sum(area_km2)),
        by = .(domain_code, longitude_bin, latitude_bin)][, scenario := scenario]
      if (!is.null(region)) {
        tab[, region_id := v[eligible, 2L * n_zone + 3L]]
        all_region[[counter]] <- tab[, .(pixels = .N, area_km2 = sum(area_km2)),
          by = .(domain_code, region_id)][, scenario := scenario]
      }
      exclusive <- which(code %in% c(2L, 3L))
      if (length(exclusive)) {
        changed <- data.table(scenario = scenario, cell = cell[exclusive],
          longitude = xy[exclusive, 1L], latitude = xy[exclusive, 2L],
          reference_zone = observed[exclusive], domain_code = code[exclusive],
          area_km2 = area[exclusive],
          finite_climate_ecotypes = rowSums(is.finite(cv[exclusive, , drop = FALSE])),
          finite_soil_ecotypes = rowSums(is.finite(sv[exclusive, , drop = FALSE])),
          strict_gate020_valid = strict_valid[exclusive, primary_index])
        changed[, coverage_reason := ifelse(domain_code == 3L,
          "no_finite_climate_zero_filled_by_soil_failure",
          ifelse(finite_soil_ecotypes == 0L, "all_soil_predictions_missing",
                 "no_finite_positive_gate_candidate"))]
        fwrite(changed, changed_file, append = file.exists(changed_file))
      }
      if (!is.null(reference_domain)) {
        ref_code <- readValues(reference_domain, row = start, nrows = rows)
        different <- eligible & !is.na(ref_code) & code != ref_code
        mask_changes[[length(mask_changes) + 1L]] <- data.table(
          scenario = scenario, different_from_normal_pixels = sum(different),
          different_from_normal_area_km2 = sum(area[different]))
      }
    }
    readStop(input)
    writeStop(domain_raster)
    if (!is.null(reference_domain)) readStop(reference_domain)
    if (scenario == "normal") reference_domain <- rast(mask_file)
    rm(input, climate)
    gc()
  }

  metric_keys <- c("scenario", "soil_gate", "domain", "semantics", "novel_rule_applied")
  result <- rbindlist(all_metrics)[, lapply(.SD, sum), by = metric_keys]
  result[, `:=`(below_0.4_pixel_share = below_0.4_pixels / valid_pixels,
    below_0.4_area_share = below_0.4_area_km2 / valid_area_km2,
    mean_max_dual = sum_max_dual / valid_pixels,
    area_weighted_mean_max_dual = sum_area_max_dual / valid_area_km2,
    mean_suitable_ecotypes = sum_suitable_ecotypes / valid_pixels,
    area_weighted_mean_suitable_ecotypes = sum_area_suitable_ecotypes / valid_area_km2)]
  setorder(result, scenario, semantics, domain, soil_gate)
  fwrite(result, file.path(output_dir, "soil_gate_native_and_common_metrics.csv"))
  zone_result <- rbindlist(all_zone_metrics)[, lapply(.SD, sum),
    by = .(scenario, soil_gate, domain, reference_zone)]
  zone_result[, `:=`(below_0.4_pixel_share = below_0.4_pixels / valid_pixels,
    below_0.4_area_share = below_0.4_area_km2 / valid_area_km2,
    mean_max_dual = sum_max_dual / valid_pixels,
    area_weighted_mean_max_dual = sum_area_max_dual / valid_area_km2,
    mean_suitable_ecotypes = sum_suitable_ecotypes / valid_pixels,
    area_weighted_mean_suitable_ecotypes = sum_area_suitable_ecotypes / valid_area_km2)]
  fwrite(zone_result, file.path(output_dir, "soil_gate_metrics_by_reference_zone.csv"))
  issues <- rbindlist(all_issues)[, lapply(.SD, sum), by = .(scenario, soil_gate)]
  fwrite(issues, file.path(output_dir, "missing_climate_zero_fill_diagnostic.csv"))
  aggregate_spatial <- function(items, keys, filename) {
    out <- rbindlist(items)[, .(pixels = sum(pixels), area_km2 = sum(area_km2)), by = keys]
    out <- merge(out, group_lookup, by = "domain_code", all.x = TRUE)
    fwrite(out, file.path(output_dir, filename))
    out
  }
  zonal <- aggregate_spatial(all_zone, c("scenario", "domain_code", "reference_zone"),
                             "domain_by_reference_zone.csv")
  aggregate_spatial(all_grid, c("scenario", "domain_code", "longitude_bin", "latitude_bin"),
                    "domain_by_1degree_geographic_bin.csv")
  national <- zonal[, .(pixels = sum(pixels), area_km2 = sum(area_km2)),
                    by = .(scenario, domain_code, domain_group)]
  fwrite(national, file.path(output_dir, "domain_national.csv"))
  zone_eligible <- zonal[, .(eligible_pixels = sum(pixels), eligible_area_km2 = sum(area_km2)),
                         by = .(scenario, reference_zone)]
  zonal <- merge(zonal, zone_eligible, by = c("scenario", "reference_zone"))
  zonal[, `:=`(share_of_reference_zone_pixels = pixels / eligible_pixels,
               share_of_reference_zone_area = area_km2 / eligible_area_km2)]
  fwrite(zonal, file.path(output_dir, "domain_by_reference_zone.csv"))
  if (length(all_region)) {
    regional <- aggregate_spatial(all_region, c("scenario", "domain_code", "region_id"),
                                  "domain_by_region.csv")
    if (!is.null(region_lookup)) {
      regional <- merge(regional, region_lookup, by = "region_id", all.x = TRUE)
      fwrite(regional, file.path(output_dir, "domain_by_region.csv"))
    }
  }
  if (length(mask_changes)) {
    changes <- rbindlist(mask_changes)[, lapply(.SD, sum), by = scenario]
    fwrite(changes, file.path(output_dir, "domain_pattern_change_from_normal.csv"))
  }

  # Confirm that input rasters actually reproduce the CSV cited in the paper.
  # A mismatch is a provenance failure; do not silently use the new numbers.
  control_failed <- FALSE
  if (file.exists(prior_csv)) {
    old <- fread(prior_csv)
    check <- result[domain == "native" & semantics == "existing_csv",
      .(scenario, soil_gate, audit_pixels = valid_pixels, audit_area_km2 = valid_area_km2,
        audit_mean_max_dual = mean_max_dual, audit_below_area_km2 = below_0.4_area_km2)]
    check <- merge(check, old[, .(scenario, soil_gate, valid_pixels, valid_area_km2,
                                  mean_max_dual, below_0.4_area_km2)],
                   by = c("scenario", "soil_gate"), all = TRUE)
    check[, matches_original_csv := !is.na(audit_pixels) & !is.na(valid_pixels) &
      audit_pixels == valid_pixels & abs(audit_area_km2 - valid_area_km2) < 0.05 &
      abs(audit_mean_max_dual - mean_max_dual) < 1e-6 &
      abs(audit_below_area_km2 - below_0.4_area_km2) < 0.05]
    fwrite(check, file.path(output_dir, "control_against_existing_soil_gate_csv.csv"))
    control_failed <- any(is.na(check$matches_original_csv) | !check$matches_original_csv)
  } else {
    warning("No assessment/soil_gate_sensitivity.csv available for provenance controls.")
  }
  unlink(file.path(output_dir, "temporary"), recursive = TRUE)
  cat("[SAVED AUDIT]", output_dir, "\n")
  if (control_failed) stop("Inputs do not reproduce existing soil_gate_sensitivity.csv; ",
                           "inspect the control CSV before interpreting any audit output.")
  cat("Domain codes 2 and 3 are separate sets; their difference is the net coverage change.\n")
  cat("No model, suitability raster, existing assessment CSV or manuscript was modified.\n")
})
