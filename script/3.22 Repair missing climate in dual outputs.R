# Repair from CACHED suitability onward. No model is loaded or fitted.
# Run before 3.23. Jobs checkpoint individually; originals are backed up.
local({
  run_repair <- function() {
    for (p in c("terra", "data.table")) if (!requireNamespace(p, quietly = TRUE)) stop("Install ", p)
    library(terra); library(data.table)
    root <- normalizePath(Sys.getenv("ECOCHINA2_DIR", getwd()), winslash = "/", mustWork = TRUE)
    while (!file.exists(file.path(root, "script", "1. data.R"))) {
      q <- dirname(root); if (q == root) stop("Run inside ecoChina2 or set ECOCHINA2_DIR."); root <- q
    }
    source(file.path(root, "functions", "dual_missing_policy.R"), local = environment())
    parse_env <- function(key, default) strsplit(Sys.getenv(key, default), ",", fixed = TRUE)[[1]]
    zones <- as.integer(parse_env("ECOCHINA2_REPAIR_ZONES", paste(c(1:7, 9:50, 52:55), collapse = ",")))
    methods <- parse_env("ECOCHINA2_REPAIR_METHODS", "rf_var,mf_var")
    scenarios <- parse_env("ECOCHINA2_REPAIR_SCENARIOS", paste(c("normal", "2011-2040SSP245",
      "2041-2070SSP245", "2071-2100SSP245", "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585"), collapse = ","))
    if (anyNA(zones) || anyDuplicated(zones) || !length(zones) ||
        any(!methods %in% c("rf_var", "mf_var"))) stop("Invalid repair settings.")
    report <- file.path(root, "assessment", "missing_climate_repair_v1")
    backup <- file.path(report, "backups")
    dir.create(file.path(report, "jobs"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(report, "temporary"), recursive = TRUE, showWarnings = FALSE)
    old_opt <- terraOptions(print = FALSE)
    terraOptions(tempdir = file.path(report, "temporary"), memfrac = 0.15)
    on.exit({
      temp <- old_opt$tempdir
      if (is.null(temp) || !dir.exists(temp)) temp <- tempdir()
      terraOptions(tempdir = temp, memfrac = old_opt$memfrac)
    }, add = TRUE)
    ref_file <- file.path(root, "raster", "ecosys_ori.tif")
    reference <- rast(ref_file)
    cell_area <- cellSize(reference, unit = "km")
    nzone <- length(zones)
    # At most ~110 MB for the widest input block on the established grid.
    rows_wide <- max(1L, min(32L, floor(14e6 / (ncol(reference) * (2 * nzone + 3)))))
    rows_dual <- max(1L, min(256L, floor(8e6 / (ncol(reference) * 4L))))
    starts_for <- function(rows) seq.int(1L, nrow(reference), by = rows)
    model_paths <- function(method, scenario) {
      soil_method <- if (method == "mf_var") "plain_mf" else "plain_rf"
      list(climate = file.path(root, "clim suitability", method, scenario,
             paste0("clim_suit_zone", zones, ".tif")),
        soil = file.path(root, "soil suitability", soil_method, "normal",
             paste0("soil_suit_zone", zones, ".tif")),
        dual = file.path(root, "dual suit", method, scenario,
             paste0("dual_suitability_zone", zones, ".tif")),
        assigned = file.path(root, "result maps", method, paste0("assigned_zone_", scenario,
             "_threshold0.4_tol1e-04_novel99_maskNA8_noNovelNormal.tif")))
    }
    check_grid <- function(files) {
      for (file in files) {
        if (!file.exists(file)) stop("Missing existing input: ", file)
        if (!compareGeom(rast(file), reference, stopOnError = FALSE)) {
          stop("Geometry mismatch; not resampled automatically: ", file)
        }
      }
    }
    jobs <- expand.grid(method = methods, scenario = scenarios, stringsAsFactors = FALSE)
    # Complete preflight before mutating any established output.
    for (i in seq_len(nrow(jobs))) {
      p <- model_paths(jobs$method[i], jobs$scenario[i])
      check_grid(c(p$climate, p$soil, p$dual, p$assigned))
    }
    # A failed/interrupted attempt must not leave a prior completion flag usable.
    unlink(file.path(report, "repair_completed.rds"))
    fwrite(data.table(code = 0:4, meaning = c("outside_reference_footprint_or_zone8",
      "all_required_climate_and_soil_candidates_available", "no_climate_candidate_available",
      "some_climate_candidates_missing_soil_may_also_be_missing", "climate_complete_but_soil_candidates_missing")),
      file.path(report, "availability_codes.csv"))

    repair_dual <- function(method, scenario, index, paths) {
      target <- paths$dual[index]
      inputs <- c(paths$climate[index], paths$soil[index], ref_file)
      if (eco_policy_current(target, inputs, "dual_gate_0.2")) return(invisible(FALSE))
      old_info <- file.info(target)
      temp <- tempfile("dual_complete_", tmpdir = dirname(target), fileext = ".tif")
      input <- c(rast(inputs[1]), rast(inputs[2]), rast(target), reference)
      output <- rast(reference); names(output) <- names(rast(target))
      writeStart(output, temp, overwrite = FALSE,
        datatype = "FLT4S", gdal = "COMPRESS=LZW")
      readStart(input)
      removed <- removed_area <- changed <- added <- nonzero_removed <- 0
      for (start in starts_for(rows_dual)) {
        nr <- min(rows_dual, nrow(reference) - start + 1L)
        v <- readValues(input, row = start, nrows = nr, mat = TRUE)
        result <- eco_strict_dual(v[, 1], v[, 2], 0.20)
        result[!is.finite(v[, 4])] <- NA_real_
        # Prove unchanged positive values throughout the ranking/map footprint,
        # including observed classes outside the modeled reference subset.
        eligible <- is.finite(v[, 4]) & v[, 4] != 8
        old_finite <- is.finite(v[, 3]); new_finite <- is.finite(result)
        removed <- removed + sum(eligible & old_finite & !new_finite)
        nonzero_removed <- nonzero_removed + sum(eligible & old_finite & !new_finite & abs(v[, 3]) > 1e-7)
        added <- added + sum(eligible & !old_finite & new_finite)
        changed <- changed + sum(eligible & old_finite & new_finite & abs(v[, 3] - result) > 1e-7,
                                  na.rm = TRUE)
        writeValues(output, result, start = start, nrows = nr)
      }
      readStop(input); writeStop(output)
      rm(input, output); gc()
      row <- data.table(method = method, scenario = scenario, zone = zones[index],
        removed_finite_cells = removed, added_finite_cells = added,
        altered_finite_values = changed, removed_nonzero_cells = nonzero_removed,
        old_dual_mtime = as.numeric(old_info$mtime), old_dual_size = old_info$size)
      fwrite(row, file.path(report, "jobs", paste0("dual_", method, "_", scenario, "_", zones[index], ".csv")))
      if (added > 0 || changed > 0 || nonzero_removed > 0) {
        unlink(temp)
        stop("Existing dual differs from the expected cached gate beyond zero-to-NA. See job CSV: ", target)
      }
      eco_commit_raster(temp, target, backup, root)
      eco_record_policy(target, inputs, "dual_gate_0.2")
      invisible(TRUE)
    }

    repair_maps <- function(method, scenario, paths) {
      status_file <- file.path(report, paste0("availability_", method, "_", scenario, ".tif"))
      inputs <- c(paths$dual, paths$climate, ref_file)
      if (eco_policy_current(status_file, inputs, "availability") &&
          eco_policy_current(paths$assigned, c(paths$dual, ref_file), "assigned_complete_candidates")) {
        return(status_file)
      }
      input <- c(rast(paths$dual), rast(paths$climate), rast(paths$assigned), reference, cell_area)
      map_temp <- tempfile("assigned_complete_", tmpdir = dirname(paths$assigned), fileext = ".tif")
      status_temp <- tempfile("availability_", tmpdir = report, fileext = ".tif")
      out_map <- rast(reference); names(out_map) <- "assigned_zone"
      out_status <- rast(reference, nlyrs = 3)
      names(out_status) <- c("availability_status", "n_dual_available", "n_climate_available")
      writeStart(out_map, map_temp, overwrite = FALSE, datatype = "INT2S", gdal = "COMPRESS=LZW")
      writeStart(out_status, status_temp, overwrite = FALSE, datatype = "INT1U", NAflag = 255, gdal = "COMPRESS=LZW")
      readStart(input)
      totals <- c(eligible_cells = 0, eligible_area_km2 = 0, climate_any_cells = 0, climate_any_area_km2 = 0,
        climate_complete_cells = 0, climate_complete_area_km2 = 0, dual_any_cells = 0, dual_any_area_km2 = 0,
        complete_cells = 0, complete_area_km2 = 0, old_assigned_cells = 0, new_assigned_cells = 0,
        removed_assigned_cells = 0, complete_but_missing_old_assignment = 0)
      grouped_pixels <- grouped_area <- numeric((max(zones) + 1L) * 10L)
      starts <- starts_for(rows_wide)
      for (b in seq_along(starts)) {
        start <- starts[b]; nr <- min(rows_wide, nrow(reference) - start + 1L)
        v <- readValues(input, row = start, nrows = nr, mat = TRUE)
        nd <- rowSums(is.finite(v[, seq_len(nzone), drop = FALSE]))
        nc <- rowSums(is.finite(v[, nzone + seq_len(nzone), drop = FALSE]))
        old <- v[, 2L * nzone + 1L]; original <- v[, 2L * nzone + 2L]; area <- v[, 2L * nzone + 3L]
        footprint <- is.finite(original) & original != 8
        eligible <- original %in% zones & is.finite(area)
        complete <- nd == nzone
        if (any(complete & nc != nzone)) stop("Invalid repaired dual availability.")
        status <- ifelse(complete, 1L, ifelse(nc == 0L, 2L, ifelse(nc < nzone, 3L, 4L)))
        status[!footprint] <- 0L
        current <- old; current[!footprint | !complete] <- NA_real_
        totals <- totals + c(sum(eligible), sum(area[eligible]),
          sum(eligible & nc > 0L), sum(area[eligible & nc > 0L]),
          sum(eligible & nc == nzone), sum(area[eligible & nc == nzone]),
          sum(eligible & nd > 0L), sum(area[eligible & nd > 0L]),
          sum(eligible & complete), sum(area[eligible & complete]),
          sum(eligible & is.finite(old)), sum(eligible & is.finite(current)),
          sum(eligible & is.finite(old) & !is.finite(current)),
          sum(eligible & complete & !is.finite(old)))
        keys <- as.integer(original[eligible]) * 10L + status[eligible] + 1L
        grouped_pixels <- grouped_pixels + tabulate(keys, nbins = length(grouped_pixels))
        if (length(keys)) {
          ss <- rowsum(area[eligible], keys, reorder = FALSE)
          ii <- as.integer(rownames(ss)); grouped_area[ii] <- grouped_area[ii] + ss[, 1]
        }
        writeValues(out_map, current, start = start, nrows = nr)
        nd[!footprint] <- NA; nc[!footprint] <- NA
        writeValues(out_status, cbind(status, nd, nc), start = start, nrows = nr)
        if (b %% 10L == 0L || b == length(starts)) {
          cat(sprintf("  [AVAILABILITY] %s %s %.0f%%\n", method, scenario, 100 * b / length(starts)))
          flush.console()
        }
      }
      readStop(input); writeStop(out_map); writeStop(out_status)
      rm(input, out_map, out_status); gc()
      summary <- as.data.table(as.list(totals)); summary[, `:=`(method = method, scenario = scenario)]
      fwrite(summary, file.path(report, "jobs", paste0("domain_", method, "_", scenario, ".csv")))
      kk <- which(grouped_pixels > 0)
      fwrite(data.table(method = method, scenario = scenario, reference_zone = (kk - 1L) %/% 10L,
        availability_status = (kk - 1L) %% 10L, pixels = grouped_pixels[kk], area_km2 = grouped_area[kk]),
        file.path(report, "jobs", paste0("zones_", method, "_", scenario, ".csv")))
      if (totals[["complete_but_missing_old_assignment"]] > 0) {
        unlink(c(map_temp, status_temp))
        stop("Some complete-input cells already lacked assignments; inspect the domain job CSV before continuing.")
      }
      eco_commit_raster(status_temp, status_file, backup, root)
      eco_record_policy(status_file, inputs, "availability")
      eco_commit_raster(map_temp, paths$assigned, backup, root)
      eco_record_policy(paths$assigned, c(paths$dual, ref_file), "assigned_complete_candidates")
      status_file
    }

    repair_rank <- function(file, status_file, paths) {
      if (!file.exists(file)) return(invisible(FALSE))
      inputs <- c(paths$dual, status_file, ref_file)
      if (eco_policy_current(file, inputs, "rank_complete_candidates")) return(invisible(FALSE))
      # Complete cells have exactly the same dual values as before. Therefore
      # their rank order, scores and summaries are preserved, not re-sorted.
      x <- rast(file); check_grid(file)
      mask <- ifel(rast(status_file)[[1]] == 1L, 1, NA_real_)
      temp <- tempfile("rank_complete_", tmpdir = dirname(file), fileext = ".tif")
      datatype <- if (basename(file) == "ranked_zone.tif") "INT2S" else "FLT4S"
      writeRaster(terra::mask(x, mask), temp, overwrite = FALSE, datatype = datatype, gdal = "COMPRESS=LZW")
      rm(x, mask); gc()
      eco_commit_raster(temp, file, backup, root)
      eco_record_policy(file, inputs, "rank_complete_candidates")
      invisible(TRUE)
    }
    for (i in seq_len(nrow(jobs))) {
      method <- jobs$method[i]; scenario <- jobs$scenario[i]
      cat(sprintf("\n[JOB %d/%d] %s %s\n", i, nrow(jobs), method, scenario)); flush.console()
      p <- model_paths(method, scenario)
      for (z in seq_along(zones)) {
        changed <- repair_dual(method, scenario, z, p)
        cat(sprintf("  [DUAL %d/%d] zone %d %s\n", z, nzone, zones[z], if (changed) "repaired" else "current"))
        flush.console()
      }
      status_file <- repair_maps(method, scenario, p)
      for (name in c("ranked_zone.tif", "ranked_suitability.tif", "ranked_summary.tif")) {
        f <- file.path(root, "dual suit ranking var", method, scenario, name)
        repair_rank(f, status_file, p)
      }
    }
    collect <- function(prefix, dest) {
      ff <- list.files(file.path(report, "jobs"), pattern = paste0("^", prefix, ".*\\.csv$"), full.names = TRUE)
      if (length(ff)) fwrite(rbindlist(lapply(ff, fread), fill = TRUE), file.path(report, dest))
    }
    collect("dual_", "dual_repair_summary.csv")
    collect("domain_", "analysis_domains.csv")
    collect("zones_", "availability_by_reference_zone.csv")
    # Reuse the verified strict-finite audit instead of repeating the eight-hour
    # sensitivity pass. This is allowed only for unchanged source rasters and
    # identical complete-candidate domains; otherwise no old table is relabeled.
    audit <- file.path(root, "assessment", "soil_domain_audit_20261101_135508")
    audit_manifest <- file.path(audit, "input_manifest.csv")
    audit_metrics <- file.path(audit, "soil_gate_native_and_common_metrics.csv")
    reused <- FALSE
    sensitivity_ready <- FALSE
    publish_sensitivity <- function(a, common) {
      commit_table <- function(table, target) {
        temp <- tempfile("sensitivity_", tmpdir = dirname(target), fileext = ".csv")
        fwrite(table, temp)
        if (file.exists(target) && identical(unname(tools::md5sum(temp)), unname(tools::md5sum(target)))) {
          unlink(temp)
        } else eco_commit_raster(temp, target, backup, root)
      }
      sensitivity <- a[, .(scenario, soil_gate, climate_only = soil_gate == 0,
        novel_rule_applied = scenario != "normal", valid_pixels,
        analogue_pixels = valid_pixels - below_0.4_pixels, below_0.4_pixels,
        valid_area_km2, analogue_area_km2 = valid_area_km2 - below_0.4_area_km2,
        below_0.4_area_km2, analogue_share = 1 - below_0.4_pixel_share,
        below_0.4_share = below_0.4_pixel_share, mean_max_dual,
        below_0.4_area_share, missing_policy = eco_missing_policy)]
      target <- file.path(root, "assessment", "soil_gate_sensitivity.csv")
      commit_table(sensitivity, target)
      common[, `:=`(missing_policy = eco_missing_policy, semantics = "complete_candidate_data")]
      commit_table(common, file.path(root, "assessment", "soil_gate_sensitivity_common_domain.csv"))
    }
    if (all(c("rf_var", "mf_var") %in% methods) && length(zones) == 53L &&
        file.exists(audit_manifest) && file.exists(audit_metrics)) {
      # The original audit wrote local wall-clock strings without a timezone.
      # Prevent fread from interpreting them as UTC before comparing mtimes.
      manifest <- fread(audit_manifest, colClasses = c(modified = "character"))
      info <- file.info(manifest$path)
      unchanged <- all(file.exists(manifest$path)) && all(info$size == manifest$size_bytes) &&
        all(abs(as.numeric(info$mtime) - as.numeric(as.POSIXct(manifest$modified,
          format = "%Y-%m-%d %H:%M:%OS", tz = ""))) < 1.1)
      metrics <- fread(audit_metrics)
      domains <- fread(file.path(report, "analysis_domains.csv"))[method == "mf_var"]
      a <- metrics[domain == "native" & semantics == "strict_finite_climate_and_soil"]
      check <- merge(domains, a[soil_gate == 0.20, .(scenario, audit_positive_cells = valid_pixels,
        audit_positive_area_km2 = valid_area_km2)], by = "scenario")
      check <- merge(check, a[soil_gate == 0, .(scenario, audit_climate_cells = valid_pixels,
        audit_climate_area_km2 = valid_area_km2)], by = "scenario")
      identical_domains <- nrow(check) == length(scenarios) &&
        all(check$complete_cells == check$audit_positive_cells) &&
        all(check$climate_complete_cells == check$audit_climate_cells) &&
        all(abs(check$complete_area_km2 - check$audit_positive_area_km2) < 0.05) &&
        all(abs(check$climate_complete_area_km2 - check$audit_climate_area_km2) < 0.05)
      common_check <- merge(metrics[domain == "common_all_gates" & semantics == "existing_csv"],
        domains[, .(scenario, complete_cells, complete_area_km2)], by = "scenario")
      identical_domains <- isTRUE(identical_domains) && nrow(common_check) == length(scenarios) * 8L &&
        all(common_check$valid_pixels == common_check$complete_cells) &&
        all(abs(common_check$valid_area_km2 - common_check$complete_area_km2) < 0.05)
      fwrite(check, file.path(report, "audit_reuse_domain_check.csv"))
      if (isTRUE(unchanged) && isTRUE(identical_domains) && nrow(a) == length(scenarios) * 8L) {
        common <- metrics[domain == "common_all_gates" & semantics == "existing_csv"]
        publish_sensitivity(a, common)
        reused <- TRUE
        sensitivity_ready <- TRUE
      }
    }
    if (length(zones) == 53L && "mf_var" %in% methods && !reused) {
      cat("[SENSITIVITY] Source/domain checks prevent audit reuse; compute national summaries only.\n")
      tables <- lapply(scenarios, function(s) {
        p <- model_paths("mf_var", s)
        eco_recompute_sensitivity(p$climate, p$soil, reference, cell_area, zones, s,
          file.path(report, "jobs", paste0("sensitivity_", s, ".csv")))
      })
      tables <- rbindlist(tables)
      publish_sensitivity(tables[domain == "native"], tables[domain == "common_all_gates"])
      sensitivity_ready <- TRUE
    }
    writeLines(c(paste("Policy:", eco_missing_policy), paste("Run:", Sys.time()),
      paste("Sensitivity reused with source/domain checks:", reused),
      paste("Sensitivity ready:", sensitivity_ready),
      "Complete-input assignments were retained exactly, including prior tie outcomes.",
      "No fitted models, climate or soil suitability rasters were changed.", capture.output(sessionInfo())),
      file.path(report, "run_info.txt"))
    completion_inputs <- c(ref_file, unlist(lapply(seq_len(nrow(jobs)), function(i) {
      p <- model_paths(jobs$method[i], jobs$scenario[i])
      c(p$climate, p$soil, p$dual, p$assigned)
    }), use.names = FALSE))
    saveRDS(list(policy = eco_missing_policy, input_files = completion_inputs,
      signature = eco_input_signature(completion_inputs), zones = zones, methods = methods,
      scenarios = scenarios, sensitivity_reused = reused, sensitivity_ready = sensitivity_ready), file.path(report, "repair_completed.rds"))
    cat("\n[REPAIR COMPLETE]", report, "\n")
    cat("Next: source('script/3.23 Refresh results after missing-climate repair.R')\n")
  }
  run_repair()
})
