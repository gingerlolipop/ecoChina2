# Pixel-wise Top-k dual-suitability analysis (k = 1, 2, 3, 4, 5)
# ==============================================================================
# This script ranks the 53 current ecotypes, calculates reference-period
# agreement, future rank retention, matched reference-to-future comparisons,
# analogue availability, and population/species Top-k areas. It does not refit
# models or redefine novel ecosystem. Zone 99 always means that every current
# ecotype has dual suitability < 0.4.
#
# Top-1 uses the assigned map from script 4 and therefore preserves its 1e-4
# tie rule. Raw rank rasters use deterministic suitability order (zone ID breaks
# exact ties). Top-k diagnostics are reference-conditioned sensitivity analyses,
# not independent future vegetation forecasts.

library(terra)
library(data.table)

rm(list = ls())
gc()


# 0. Paths and settings =========================================================

find_project_root <- function() {
  configured <- Sys.getenv("ECOCHINA2_DIR", unset = "")
  if (nzchar(configured)) return(normalizePath(configured, mustWork = TRUE))

  current <- normalizePath(getwd(), mustWork = TRUE)
  if (file.exists(file.path(current, "data", "ecosys_ori.tif")) ||
      file.exists(file.path(current, "raster", "ecosys_ori.tif"))) return(current)
  if (basename(current) == "script") return(dirname(current))
  stop("Run from the repository root or set ECOCHINA2_DIR.")
}

base_dir <- find_project_root()

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "true" else "false")
  tolower(trimws(value)) %in% c("1", "true", "yes", "y")
}

modeled_zones <- c(1:7, 9:50, 52:55)
k_values <- 1:5
dual_threshold <- 0.4
rank_min_suitability <- 0
novel_value <- 99L
force_ranking <- env_flag("ECOCHINA2_FORCE_RANKING", FALSE)

scenarios <- c(
  "normal",
  "2011-2040SSP245", "2041-2070SSP245", "2071-2100SSP245",
  "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585"
)
future_scenarios <- setdiff(scenarios, "normal")

reference_candidates <- c(
  file.path(base_dir, "raster", "ecosys_ori.tif"),
  file.path(base_dir, "data", "ecosys_ori.tif")
)
reference_file <- reference_candidates[file.exists(reference_candidates)][1]
if (is.na(reference_file)) stop("Missing data/ecosys_ori.tif.")
dual_root <- file.path(base_dir, "dual suit", "mf_var")
map_root <- file.path(base_dir, "result maps", "mf_var")
rank_root <- file.path(base_dir, "dual suit ranking var", "mf_var")
output_root <- file.path(base_dir, "assessment", "topk")
table_dir <- file.path(output_root, "tables")
cache_root <- file.path(output_root, "cache")
reference_map_cache <- file.path(cache_root, "reference_topk_maps")
future_map_cache <- file.path(cache_root, "future_topk_maps")
summary_cache <- file.path(cache_root, "change_summaries")
legacy_reference_cache <- file.path(
  base_dir, "assessment_var", "future_topk_matched", "cache", "reference_topk_maps"
)
legacy_summary_cache <- file.path(
  base_dir, "assessment_var", "future_topk_matched", "cache", "change_summaries"
)
tmp_dir <- file.path(base_dir, "tmp_dual_suit_ranking_var")

for (directory in c(
  table_dir, cache_root, reference_map_cache,
  future_map_cache, summary_cache, tmp_dir
)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

terraOptions(tempdir = tmp_dir, memfrac = 0.20)
detected_cores <- parallel::detectCores()
if (is.na(detected_cores)) detected_cores <- 1L
rank_cores <- min(4L, max(1L, detected_cores - 1L))


# 1. Helpers ===================================================================

first_existing <- function(paths, label, required = TRUE) {
  found <- paths[file.exists(paths)]
  if (length(found)) return(found[[1]])
  if (required) stop("Missing ", label, ":\n", paste(paths, collapse = "\n"))
  NA_character_
}

dual_file <- function(scenario, zone) {
  first_existing(
    c(
      file.path(dual_root, scenario, paste0("dual_suitability_zone", zone, ".tif")),
      file.path(dual_root, scenario, paste0("dual_zone", zone, ".tif"))
    ),
    paste("dual suitability", scenario, "zone", zone)
  )
}

assigned_map_file <- function(scenario) {
  first_existing(
    c(
      file.path(
        map_root,
        paste0(
          "assigned_zone_", scenario,
          "_threshold0.4_tol1e-04_novel99_maskNA8_noNovelNormal.tif"
        )
      ),
      file.path(map_root, paste0("assigned_zone_", scenario, ".tif"))
    ),
    paste("assigned map", scenario)
  )
}

scenario_fields <- function(scenario) {
  if (scenario == "normal") return(list(period = "1961-1990", ssp = "Reference"))
  list(
    period = sub("SSP.*$", "", scenario),
    ssp = sub("^.*(SSP[0-9]+)$", "\\1", scenario)
  )
}

check_geometry <- function(x, template, label) {
  if (!compareGeom(x, template, stopOnError = FALSE)) {
    stop(label, " does not match data/ecosys_ori.tif geometry.")
  }
  invisible(TRUE)
}

rank_paths <- function(scenario) {
  output_dir <- file.path(rank_root, scenario)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  list(
    output_dir = output_dir,
    zone = file.path(output_dir, "ranked_zone.tif"),
    suitability = file.path(output_dir, "ranked_suitability.tif"),
    summary = file.path(output_dir, "ranked_summary.tif")
  )
}

valid_rank_outputs <- function(paths, template) {
  if (!all(file.exists(unlist(paths[c("zone", "suitability", "summary")])))) return(FALSE)
  tryCatch({
    zone <- rast(paths$zone)
    suitability <- rast(paths$suitability)
    summary <- rast(paths$summary)
    expected_zone_names <- paste0("rank", 1:5, "_zone")
    expected_suitability_names <- paste0("rank", 1:5, "_suit")
    zone_names <- names(zone)
    suitability_names <- names(suitability)
    summary_names <- names(summary)
    compareGeom(zone, template, stopOnError = FALSE) &&
      compareGeom(suitability, template, stopOnError = FALSE) &&
      compareGeom(summary, template, stopOnError = FALSE) &&
      nlyr(zone) >= 5L && nlyr(suitability) >= 5L && nlyr(summary) >= 4L &&
      identical(zone_names[1:5], expected_zone_names) &&
      identical(suitability_names[1:5], expected_suitability_names) &&
      any(c("n_zone_above_threshold", "n_suitable") %in% summary_names)
  }, error = function(e) FALSE)
}

raster_valid <- function(file, template, minimum_layers = 1L) {
  if (!file.exists(file)) return(FALSE)
  tryCatch({
    x <- rast(file)
    nlyr(x) >= minimum_layers &&
      compareGeom(x, template, stopOnError = FALSE)
  }, error = function(e) FALSE)
}

files_current <- function(outputs, inputs) {
  outputs <- outputs[!is.na(outputs) & nzchar(outputs)]
  if (!length(outputs) || !all(file.exists(outputs))) return(FALSE)
  inputs <- inputs[!is.na(inputs) & nzchar(inputs) & file.exists(inputs)]
  if (!length(inputs)) return(TRUE)
  isTRUE(min(file.info(outputs)$mtime) >= max(file.info(inputs)$mtime))
}

reference_cache_file <- function(k, legacy = FALSE) {
  root <- if (legacy) legacy_reference_cache else reference_map_cache
  file.path(
    root,
    paste0("mf_var_normal_top", k, "_common_domain_v2.tif")
  )
}

matched_summary_file <- function(scenario, k, legacy = FALSE) {
  root <- if (legacy) legacy_summary_cache else summary_cache
  suffix <- if (k == 1L) "" else "_common_domain_v2"
  file.path(
    root,
    paste0("mf_var_", scenario, "_top", k, suffix, "_matched_change.csv")
  )
}

table_has_all_k <- function(file, expected_k = k_values) {
  if (!file.exists(file)) return(FALSE)
  tryCatch({
    table <- fread(file, select = "k")
    setequal(sort(unique(table$k[!is.na(table$k)])), expected_k)
  }, error = function(e) FALSE)
}

normalize_matched_summary <- function(table, scenario, fields, k) {
  if ("pixel_count" %in% names(table) && !"pixels" %in% names(table)) {
    setnames(table, "pixel_count", "pixels")
  }
  if ("change_code" %in% names(table) && !"code" %in% names(table)) {
    setnames(table, "change_code", "code")
  }
  if (!"code" %in% names(table)) {
    table[, code := fcase(
      grepl("novel", change_class, ignore.case = TRUE), 3L,
      grepl("different|changed", change_class, ignore.case = TRUE), 2L,
      default = 1L
    )]
  }
  table[, `:=`(
    code = as.integer(code),
    pixels = as.numeric(pixels),
    area_km2 = as.numeric(area_km2)
  )]
  table <- table[code %in% 1:3, .(
    pixels = sum(pixels, na.rm = TRUE),
    area_km2 = sum(area_km2, na.rm = TRUE)
  ), by = code]
  table <- merge(data.table(code = 1:3), table, by = "code", all.x = TRUE)
  table[is.na(pixels), pixels := 0]
  table[is.na(area_km2), area_km2 := 0]
  table[, change_class := c("stable", "changed", "novel")[as.integer(code)]]
  table[, `:=`(
    scenario = scenario,
    period = fields$period,
    ssp = fields$ssp,
    k = as.integer(k),
    pixels = as.numeric(pixels),
    area_km2 = as.numeric(area_km2)
  )]
  table[, .(code, pixels, area_km2, change_class, scenario, period, ssp, k)]
}

matched_cache_valid <- function(file, inputs = character()) {
  if (!file.exists(file)) return(FALSE)
  tryCatch({
    table <- fread(file)
    columns <- names(table)
    has_code <- any(c("code", "change_code", "change_class") %in% columns)
    has_pixels <- any(c("pixels", "pixel_count") %in% columns)
    if (!(has_code && has_pixels && "area_km2" %in% columns)) return(FALSE)
    codes <- if ("code" %in% columns) {
      as.integer(table$code)
    } else if ("change_code" %in% columns) {
      as.integer(table$change_code)
    } else {
      fcase(
        grepl("novel", table$change_class, ignore.case = TRUE), 3L,
        grepl("different|changed", table$change_class, ignore.case = TRUE), 2L,
        default = 1L
      )
    }
    pixel_column <- if ("pixels" %in% columns) "pixels" else "pixel_count"
    nrow(table) == 3L && !anyDuplicated(codes) && setequal(codes, 1:3) &&
      all(is.finite(as.numeric(table[[pixel_column]]))) &&
      all(as.numeric(table[[pixel_column]]) >= 0) &&
      all(is.finite(as.numeric(table$area_km2))) &&
      all(as.numeric(table$area_km2) >= 0) && files_current(file, inputs)
  }, error = function(e) FALSE)
}

mask_stats <- function(mask, cell_area) {
  pixels <- global(mask, "sum", na.rm = TRUE)[1, 1]
  area <- global(ifel(mask, cell_area, 0), "sum", na.rm = TRUE)[1, 1]
  data.table(
    pixels = ifelse(is.na(pixels), 0, as.numeric(pixels)),
    area_km2 = ifelse(is.na(area), 0, as.numeric(area))
  )
}

coded_area <- function(code_raster, cell_area) {
  pixel <- as.data.table(freq(code_raster, value = TRUE))
  if (!nrow(pixel)) return(data.table())
  setnames(pixel, tail(names(pixel), 2), c("code", "pixels"))
  pixel <- pixel[, .(code = as.integer(code), pixels = as.numeric(pixels))]

  area <- as.data.table(zonal(cell_area, code_raster, fun = "sum", na.rm = TRUE))
  if (!nrow(area)) return(data.table())
  setnames(area, names(area)[1:2], c("code", "area_km2"))
  area[, code := as.integer(code)]
  merge(pixel, area, by = "code", all = TRUE)
}

reference_hit <- function(reference, ranked_zone, k) {
  hit <- !is.na(ranked_zone[[1]]) & ranked_zone[[1]] == reference
  if (k > 1L) {
    for (rank_index in 2:k) {
      hit <- hit |
        (!is.na(ranked_zone[[rank_index]]) & ranked_zone[[rank_index]] == reference)
    }
  }
  hit
}

build_reference_topk_map <- function(reference, assigned, ranked_zone, k) {
  if (k == 1L) return(assigned)
  hit <- reference_hit(reference, ranked_zone, k)
  ifel(hit, reference, assigned)
}

first_suitable_rank <- function(ranked_zone, ranked_suitability, zone_values) {
  first_rank <- ranked_zone[[1]] * NA_real_
  for (rank_index in k_values) {
    hit <- ranked_zone[[rank_index]] %in% zone_values &
      ranked_suitability[[rank_index]] >= dual_threshold
    first_rank <- ifel(is.na(first_rank) & hit, rank_index, first_rank)
  }
  names(first_rank) <- "first_suitable_rank"
  first_rank
}

rank_area_rows <- function(first_rank, cell_area, identifiers = list()) {
  exact <- coded_area(first_rank, cell_area)
  exact_area <- setNames(rep(0, 5), k_values)
  exact_pixels <- setNames(rep(0, 5), k_values)
  if (nrow(exact)) {
    keep <- exact$code %in% k_values
    exact_area[as.character(exact$code[keep])] <- exact$area_km2[keep]
    exact_pixels[as.character(exact$code[keep])] <- exact$pixels[keep]
  }

  result <- rbindlist(lapply(k_values, function(k) {
    row <- data.table(
      k = k,
      pixels = sum(exact_pixels[seq_len(k)]),
      area_km2 = sum(exact_area[seq_len(k)])
    )
    for (name in names(identifiers)) row[, (name) := identifiers[[name]]]
    row
  }), fill = TRUE)
  result[]
}

make_rank_cell <- function(zone_ids, rank_min, threshold, top_k) {
  force(zone_ids)
  force(rank_min)
  force(threshold)
  force(top_k)

  function(x) {
    finite <- which(is.finite(x))
    if (!length(finite)) return(rep(NA_real_, 2L * top_k + 6L))

    rankable <- finite[x[finite] > rank_min]
    ordered <- if (length(rankable)) {
      rankable[order(-x[rankable], zone_ids[rankable])]
    } else integer()
    top <- head(ordered, top_k)

    zone_out <- rep(NA_real_, top_k)
    suitability_out <- rep(NA_real_, top_k)
    if (length(top)) {
      zone_out[seq_along(top)] <- zone_ids[top]
      suitability_out[seq_along(top)] <- x[top]
    }

    highest <- if (length(ordered)) x[ordered[[1]]] else NA_real_
    second <- if (length(ordered) >= 2L) x[ordered[[2]]] else NA_real_
    margin <- if (is.na(second)) NA_real_ else highest - second

    c(
      zone_out,
      suitability_out,
      n_zone_ranked = length(rankable),
      n_zone_above_threshold = sum(x[finite] >= threshold),
      top1_minus_top2 = margin,
      top1_suit = highest,
      top2_suit = second,
      novel_by_threshold = as.numeric(!any(x[finite] >= threshold))
    )
  }
}


# 2. Rank all 53 dual-suitability surfaces =====================================

reference <- rast(reference_file)
cell_area_cache <- NULL
get_cell_area <- function() {
  if (is.null(cell_area_cache)) {
    area <- cellSize(reference, unit = "km")
    names(area) <- "cell_area_km2"
    cell_area_cache <<- area
  }
  cell_area_cache
}
ranking_mask_cache <- NULL
get_ranking_mask <- function() {
  if (is.null(ranking_mask_cache)) {
    # Preserve the original ranking domain: exclude reference NA and Zone 8.
    ranking_mask_cache <<- ifel(!is.na(reference) & reference != 8, 1, NA)
  }
  ranking_mask_cache
}

rank_index <- list()

for (scenario in scenarios) {
  cat("\n[RANK]", scenario, "\n")
  paths <- rank_paths(scenario)

  if (valid_rank_outputs(paths, reference) && !force_ranking) {
    cat("[REUSE]", paths$output_dir, "\n")
  } else {
    existing_rank_files <- unlist(paths[c("zone", "suitability", "summary")])
    if (any(file.exists(existing_rank_files)) && !force_ranking) {
      stop(
        "Existing legacy rank output is incomplete or invalid; it was not overwritten:\n",
        paste(existing_rank_files[file.exists(existing_rank_files)], collapse = "\n"),
        "\nSet ECOCHINA2_FORCE_RANKING=true only when these ranks must be rebuilt."
      )
    }

    files <- vapply(modeled_zones, function(z) dual_file(scenario, z), character(1))
    dual_stack <- rast(files)
    names(dual_stack) <- paste0("zone", modeled_zones)
    if (!compareGeom(dual_stack, reference, stopOnError = FALSE)) {
      dual_stack <- resample(dual_stack, reference, method = "bilinear")
      names(dual_stack) <- paste0("zone", modeled_zones)
    }
    dual_stack <- mask(dual_stack, get_ranking_mask())

    rank_cell <- make_rank_cell(
      modeled_zones, rank_min_suitability, dual_threshold, 10L
    )

    ranked_all <- app(
      dual_stack,
      rank_cell,
      cores = rank_cores,
      filename = file.path(paths$output_dir, "rank_all.tmp.tif"),
      overwrite = TRUE,
      wopt = list(gdal = c("COMPRESS=LZW"))
    )
    names(ranked_all) <- c(
      paste0("rank", 1:10, "_zone"),
      paste0("rank", 1:10, "_suit"),
      "n_zone_ranked", "n_zone_above_threshold", "top1_minus_top2",
      "top1_suit", "top2_suit", "novel_by_threshold"
    )

    writeRaster(
      ranked_all[[1:10]], paths$zone, overwrite = force_ranking,
      datatype = "INT2S", NAflag = -32768,
      wopt = list(gdal = c("COMPRESS=LZW"))
    )
    writeRaster(
      ranked_all[[11:20]], paths$suitability, overwrite = force_ranking,
      datatype = "FLT4S", wopt = list(gdal = c("COMPRESS=LZW"))
    )
    writeRaster(
      ranked_all[[21:26]], paths$summary, overwrite = force_ranking,
      datatype = "FLT4S", wopt = list(gdal = c("COMPRESS=LZW"))
    )

    unlink(file.path(paths$output_dir, "rank_all.tmp.tif"))
    rm(dual_stack, ranked_all)
    gc()
  }

  rank_index[[length(rank_index) + 1L]] <- data.table(
    scenario = scenario,
    ranked_zone = paths$zone,
    ranked_suitability = paths$suitability,
    rank_summary = paths$summary
  )
}

rank_index <- rbindlist(rank_index)
fwrite(rank_index, file.path(table_dir, "rank_index.csv"))


# 3. Reference-period Top-k agreement and confusion ============================

normal_assigned_file <- assigned_map_file("normal")
normal_rank_file <- rank_paths("normal")$zone
normal_assigned <- rast(normal_assigned_file)
normal_ranked <- rast(normal_rank_file)[[1:5]]
check_geometry(normal_assigned, reference, "Normal assigned map")

valid_reference_cache <- NULL
get_valid_reference <- function() {
  if (is.null(valid_reference_cache)) {
    valid_reference_cache <<-
      reference %in% modeled_zones & !is.na(normal_assigned)
  }
  valid_reference_cache
}

reference_topk_files <- list()

for (k in k_values) {
  if (k == 1L) {
    reference_topk_files[[as.character(k)]] <- normal_assigned_file
    next
  }

  legacy_file <- reference_cache_file(k, legacy = TRUE)
  output_file <- reference_cache_file(k)
  reference_cache_inputs <- c(reference_file, normal_assigned_file, normal_rank_file)
  if (raster_valid(legacy_file, reference) &&
      files_current(legacy_file, reference_cache_inputs)) {
    reference_topk_files[[as.character(k)]] <- legacy_file
  } else if (raster_valid(output_file, reference) &&
             files_current(output_file, reference_cache_inputs)) {
    reference_topk_files[[as.character(k)]] <- output_file
  } else {
    diagnostic_map <- build_reference_topk_map(
      reference, normal_assigned, normal_ranked, k
    )
    diagnostic_map <- ifel(get_valid_reference(), diagnostic_map, NA)
    names(diagnostic_map) <- paste0("reference_top", k)
    writeRaster(
      diagnostic_map, output_file, overwrite = TRUE,
      datatype = "INT2S", NAflag = -32768,
      wopt = list(gdal = "COMPRESS=LZW")
    )
    reference_topk_files[[as.character(k)]] <- output_file
    rm(diagnostic_map)
  }
}

reference_agreement_file <- file.path(table_dir, "reference_agreement.csv")
reference_confusion_file <- file.path(table_dir, "reference_confusion.csv")

if (table_has_all_k(reference_agreement_file) &&
    table_has_all_k(reference_confusion_file) &&
    files_current(
      c(reference_agreement_file, reference_confusion_file),
      c(reference_file, unlist(reference_topk_files))
    )) {
  reference_agreement <- fread(reference_agreement_file)
  reference_confusion <- fread(reference_confusion_file)
  cat("[REUSE TABLE] reference Top-k agreement and confusion\n")
} else {
  valid_reference <- get_valid_reference()
  reference_total <- mask_stats(valid_reference, get_cell_area())
  reference_agreement_results <- list()
  reference_confusion_results <- list()

  for (k in k_values) {
    diagnostic_map <- rast(reference_topk_files[[as.character(k)]])
    diagnostic_map <- ifel(valid_reference, diagnostic_map, NA)
    matched <- valid_reference & diagnostic_map == reference
    matched_stats <- mask_stats(matched, get_cell_area())
    reference_agreement_results[[length(reference_agreement_results) + 1L]] <- data.table(
      k = k,
      total_pixels = reference_total$pixels,
      matched_pixels = matched_stats$pixels,
      pixel_agreement = matched_stats$pixels / reference_total$pixels,
      total_area_km2 = reference_total$area_km2,
      matched_area_km2 = matched_stats$area_km2,
      area_agreement = matched_stats$area_km2 / reference_total$area_km2
    )
    confusion_code <- ifel(valid_reference, reference * 1000L + diagnostic_map, NA)
    flow <- coded_area(confusion_code, get_cell_area())
    if (nrow(flow)) {
      flow[, `:=`(
        k = k,
        observed_zone = as.integer(code %/% 1000L),
        assigned_zone = as.integer(code %% 1000L)
      )]
      reference_confusion_results[[length(reference_confusion_results) + 1L]] <- flow
    }
    rm(diagnostic_map, matched, confusion_code, flow)
  }

  reference_agreement <- rbindlist(reference_agreement_results)
  reference_confusion <- rbindlist(reference_confusion_results, fill = TRUE)
  fwrite(reference_agreement, reference_agreement_file)
  fwrite(reference_confusion, reference_confusion_file)
}


# 4. Future rank retention and matched Top-k comparison ========================

future_retention <- list()
future_retention_by_zone <- list()
matched_change <- list()

for (scenario in future_scenarios) {
  cat("\n[FUTURE TOP-K]", scenario, "\n")
  fields <- scenario_fields(scenario)
  future_assigned_file <- assigned_map_file(scenario)
  future_rank_file <- rank_paths(scenario)$zone
  future_assigned <- rast(future_assigned_file)
  future_ranked <- rast(future_rank_file)[[1:5]]
  check_geometry(future_assigned, reference, paste("Future assigned", scenario))

  valid <- novel <- stable <- NULL

  retention_cache_file <- file.path(
    cache_root, paste0("mf_var_", scenario, "_future_retention.csv")
  )
  retention_zone_cache_file <- file.path(
    cache_root, paste0("mf_var_", scenario, "_future_retention_by_zone.csv")
  )

  if (table_has_all_k(retention_cache_file) &&
      table_has_all_k(retention_zone_cache_file) &&
      files_current(
        c(retention_cache_file, retention_zone_cache_file),
        c(reference_file, normal_assigned_file, future_assigned_file, future_rank_file)
      )) {
    future_retention[[length(future_retention) + 1L]] <- fread(retention_cache_file)
    future_retention_by_zone[[length(future_retention_by_zone) + 1L]] <-
      fread(retention_zone_cache_file)
    cat("[REUSE TABLE] future retention |", scenario, "\n")
  } else {
    valid <- !is.na(normal_assigned) & !is.na(future_assigned)
    novel <- valid & future_assigned == novel_value
    stable <- valid & !novel & future_assigned == normal_assigned
    scenario_retention <- list()
    scenario_retention_by_zone <- list()
    first_raw_rank <- future_ranked[[1]] * NA_real_
    for (rank_index_value in k_values) {
      hit <- !is.na(future_ranked[[rank_index_value]]) &
        future_ranked[[rank_index_value]] == normal_assigned
      first_raw_rank <- ifel(is.na(first_raw_rank) & hit, rank_index_value, first_raw_rank)
    }

    retention_class <- ifel(
      novel, 7L,
      ifel(stable, 1L,
        ifel(first_raw_rank == 1L, 2L,
          ifel(first_raw_rank == 2L, 3L,
            ifel(first_raw_rank == 3L, 4L,
              ifel(first_raw_rank == 4L, 5L,
                ifel(first_raw_rank == 5L, 6L, 8L)
              )
            )
          )
        )
      )
    )
    class_labels <- data.table(
      class_code = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
      retention_class = c(
        "stable_top1", "former_raw_rank1_after_tie", "former_rank2",
        "former_rank3", "former_rank4", "former_rank5", "novel", "below_top5"
      )
    )
    class_table <- coded_area(retention_class, get_cell_area())
    if (nrow(class_table)) {
      setnames(class_table, "code", "class_code")
      class_table <- merge(class_table, class_labels, by = "class_code", all.x = TRUE)
      class_table[, `:=`(scenario = scenario, period = fields$period, ssp = fields$ssp)]
      scenario_retention[[length(scenario_retention) + 1L]] <- class_table
    }

    valid_stats <- mask_stats(valid, get_cell_area())
    novel_stats <- mask_stats(novel, get_cell_area())
    for (k in k_values) {
      retained <- if (k == 1L) stable else
        stable | reference_hit(normal_assigned, future_ranked, k)
      retained <- valid & !novel & retained
      retained_stats <- mask_stats(retained, get_cell_area())
      scenario_retention[[length(scenario_retention) + 1L]] <- data.table(
        class_code = NA_integer_, retention_class = "cumulative_topk",
        scenario = scenario, period = fields$period, ssp = fields$ssp, k = k,
        pixels = retained_stats$pixels, area_km2 = retained_stats$area_km2,
        total_pixels = valid_stats$pixels, total_area_km2 = valid_stats$area_km2,
        retained_pixel_share = retained_stats$pixels / valid_stats$pixels,
        retained_area_share = retained_stats$area_km2 / valid_stats$area_km2,
        novel_pixels = novel_stats$pixels, novel_area_km2 = novel_stats$area_km2
      )
      status <- ifel(novel, 3L, ifel(retained, 1L, 2L))
      zone_status <- ifel(valid, normal_assigned * 10L + status, NA)
      by_zone <- coded_area(zone_status, get_cell_area())
      if (nrow(by_zone)) {
        by_zone[, `:=`(
          normal_zone = as.integer(code %/% 10L),
          status_code = as.integer(code %% 10L),
          status = c("retained", "not_retained", "novel")[status_code],
          scenario = scenario, period = fields$period, ssp = fields$ssp, k = k
        )]
        scenario_retention_by_zone[[length(scenario_retention_by_zone) + 1L]] <- by_zone
      }
      rm(retained, status, zone_status, by_zone)
    }

    scenario_retention <- rbindlist(scenario_retention, fill = TRUE)
    scenario_retention_by_zone <- rbindlist(scenario_retention_by_zone, fill = TRUE)
    fwrite(scenario_retention, retention_cache_file)
    fwrite(scenario_retention_by_zone, retention_zone_cache_file)
    future_retention[[length(future_retention) + 1L]] <- scenario_retention
    future_retention_by_zone[[length(future_retention_by_zone) + 1L]] <-
      scenario_retention_by_zone
    rm(first_raw_rank, retention_class)
  }

  for (k in k_values) {
    legacy_file <- matched_summary_file(scenario, k, legacy = TRUE)
    output_file <- matched_summary_file(scenario, k)
    future_map_inputs <- c(
      reference_file,
      normal_assigned_file,
      future_assigned_file,
      if (k > 1L) future_rank_file else character()
    )
    matched_inputs <- c(
      future_map_inputs,
      reference_topk_files[[as.character(k)]]
    )
    cache_file <- if (matched_cache_valid(legacy_file, matched_inputs)) {
      legacy_file
    } else if (matched_cache_valid(output_file, matched_inputs)) {
      output_file
    } else {
      NA_character_
    }
    future_topk_dir <- file.path(future_map_cache, scenario)
    dir.create(future_topk_dir, recursive = TRUE, showWarnings = FALSE)
    future_topk_file <- file.path(future_topk_dir, paste0("future_top", k, ".tif"))
    future_topk_current <- k > 1L &&
      raster_valid(future_topk_file, reference) &&
      files_current(future_topk_file, future_map_inputs)

    # A valid matched summary needs no raster algebra. For k > 1, still create
    # a missing display map once so script 11 can draw the complete 1:5 panel.
    if (!is.na(cache_file) && (k == 1L || future_topk_current)) {
      change_table <- normalize_matched_summary(
        fread(cache_file), scenario, fields, k
      )
      matched_change[[length(matched_change) + 1L]] <- change_table
      cat("[REUSE MATCHED]", scenario, "| k =", k, "\n")
      next
    }

    if (k == 1L) {
      # Top-1 is already the assigned map from script 4. Keep it in memory for
      # the matched summary instead of writing a duplicate national raster.
      if (is.null(valid)) {
        valid <- !is.na(normal_assigned) & !is.na(future_assigned)
      }
      future_topk <- ifel(valid, future_assigned, NA)
    } else if (future_topk_current) {
      future_topk <- rast(future_topk_file)
    } else {
      if (is.null(valid)) {
        valid <- !is.na(normal_assigned) & !is.na(future_assigned)
      }
      observed_hit <- reference_hit(reference, future_ranked, k)
      future_topk <- ifel(
        future_assigned == novel_value,
        novel_value,
        ifel(observed_hit, reference, future_assigned)
      )
      future_topk <- ifel(valid, future_topk, NA)
      writeRaster(
        future_topk,
        future_topk_file,
        overwrite = TRUE,
        datatype = "INT2S",
        NAflag = -32768,
        wopt = list(gdal = "COMPRESS=LZW")
      )
    }

    if (!is.na(cache_file)) {
      change_table <- normalize_matched_summary(
        fread(cache_file), scenario, fields, k
      )
      cat("[REUSE MATCHED]", scenario, "| k =", k, "\n")
    } else {
      reference_topk <- rast(reference_topk_files[[as.character(k)]])
      change_class_raster <- ifel(
        future_topk == novel_value,
        3L,
        ifel(future_topk == reference_topk, 1L, 2L)
      )
      change_table <- coded_area(change_class_raster, get_cell_area())
      change_table <- normalize_matched_summary(change_table, scenario, fields, k)
      fwrite(change_table, output_file)
      rm(reference_topk, change_class_raster)
    }
    matched_change[[length(matched_change) + 1L]] <- change_table
    rm(future_topk)
  }

  rm(
    list = intersect(
      c("future_assigned", "future_ranked", "valid", "novel", "stable",
        "hit", "observed_hit"),
      ls(envir = environment())
    ),
    envir = environment()
  )
  gc()
}

future_retention <- rbindlist(future_retention, fill = TRUE)
future_retention_by_zone <- rbindlist(future_retention_by_zone, fill = TRUE)
matched_change <- rbindlist(matched_change, fill = TRUE)

expected_matched_rows <- length(future_scenarios) * length(k_values) * 3L
if (nrow(matched_change) != expected_matched_rows ||
    anyDuplicated(matched_change[, .(scenario, k, code)])) {
  stop("Matched Top-k summary is incomplete or contains duplicate scenario/k/code rows.")
}

fwrite(future_retention, file.path(table_dir, "future_retention.csv"))
fwrite(future_retention_by_zone, file.path(table_dir, "future_retention_by_zone.csv"))
fwrite(matched_change, file.path(table_dir, "matched_change.csv"))


# 5. Suitable-analogue availability and source-zone Top-k area =================

analogue_availability <- list()
zone_topk_area <- list()
species_topk_area <- list()

population_file <- first_existing(
  c(
    file.path(
      base_dir, "future tree niche dual suitability var", "tables",
      "dual_population_projection_lookup_var.csv"
    ),
    file.path(base_dir, "future tree niche var", "tables", "population_projection_lookup_var.csv"),
    file.path(base_dir, "species_zone_population_long.csv"),
    file.path(base_dir, "data", "species_zone_population_long.csv"),
    file.path(base_dir, "data", "processed", "species_zone_population_long.csv")
  ),
  "population lookup",
  required = FALSE
)

population <- NULL
if (!is.na(population_file)) {
  population <- fread(population_file)
  if (!"source_zone" %in% names(population)) {
    if (!"Zone" %in% names(population)) {
      stop("Population lookup needs source_zone or Zone.")
    }
    population[, source_zone := as.integer(gsub("[^0-9]", "", as.character(Zone)))]
  }
  if (!"Species" %in% names(population)) stop("Population lookup needs Species.")
  if (!"projection_eligible" %in% names(population)) {
    population[, projection_eligible := if ("projected" %in% names(population)) {
      as.logical(projected)
    } else {
      source_zone %in% modeled_zones
    }]
  }
  population <- population[projection_eligible == TRUE]
  if (anyDuplicated(population[, .(Species, source_zone)])) {
    stop("Population lookup contains duplicated Species x source_zone rows.")
  }
}

analogue_file <- file.path(table_dir, "analogue_availability.csv")
zone_area_file <- file.path(table_dir, "source_zone_area.csv")
species_area_file <- file.path(table_dir, "species_area.csv")
population_area_file <- file.path(table_dir, "population_area.csv")
area_cache_dir <- file.path(cache_root, "area")
dir.create(area_cache_dir, recursive = TRUE, showWarnings = FALSE)

area_cache_file <- function(kind, scenario) {
  file.path(area_cache_dir, paste0(kind, "_", scenario, ".csv"))
}

complete_area_rows <- function(rows, k_column, expected_rows,
                               key_columns, numeric_columns,
                               expected_species = NULL,
                               expected_zones = NULL) {
  if (is.null(rows) || !nrow(rows) ||
      !all(c("scenario", k_column, key_columns, numeric_columns) %in% names(rows))) {
    return(FALSE)
  }
  observed_k <- suppressWarnings(as.numeric(rows[[k_column]]))
  if (any(!is.finite(observed_k)) || any(observed_k != floor(observed_k)) ||
      !setequal(sort(unique(observed_k)), k_values)) return(FALSE)
  if (!is.null(expected_rows) && nrow(rows) != expected_rows) return(FALSE)
  if (anyNA(rows[, ..key_columns]) || anyDuplicated(rows[, ..key_columns])) return(FALSE)
  if ("source_zone" %in% key_columns) {
    observed_zone <- suppressWarnings(as.numeric(rows$source_zone))
    if (any(!is.finite(observed_zone)) ||
        any(observed_zone != floor(observed_zone))) return(FALSE)
  }
  numeric_ok <- vapply(numeric_columns, function(column) {
    values <- suppressWarnings(as.numeric(rows[[column]]))
    integer_ok <- !grepl("pixel", column, ignore.case = TRUE) ||
      all(values == floor(values))
    all(is.finite(values)) && all(values >= 0) && integer_ok
  }, logical(1))
  if (!all(numeric_ok)) return(FALSE)
  if (!is.null(expected_species)) {
    if (!"Species" %in% names(rows) ||
        !setequal(sort(unique(as.character(rows$Species))), expected_species)) {
      return(FALSE)
    }
  }
  if (!is.null(expected_zones)) {
    if (!"source_zone" %in% names(rows) ||
        !setequal(sort(unique(as.integer(rows$source_zone))), expected_zones)) {
      return(FALSE)
    }
  }
  TRUE
}

read_area_cache <- function(cache_file, global_file, scenario_value, inputs,
                            k_column, expected_rows, key_columns,
                            numeric_columns, expected_species = NULL,
                            expected_zones = NULL) {
  for (candidate in unique(c(cache_file, global_file))) {
    if (!files_current(candidate, inputs)) next
    table <- tryCatch(fread(candidate), error = function(e) NULL)
    if (is.null(table) || !"scenario" %in% names(table)) next
    rows <- table[table[["scenario"]] == scenario_value]
    if (!complete_area_rows(
      rows, k_column, expected_rows, key_columns, numeric_columns,
      expected_species, expected_zones
    )) next
    if (!identical(candidate, cache_file)) fwrite(rows, cache_file)
    return(rows)
  }
  NULL
}

has_population <- !is.null(population) && nrow(population) > 0L
expected_species <- if (has_population) {
  sort(unique(as.character(population$Species)))
} else {
  character()
}

analogue_availability <- list()
zone_topk_area <- list()
species_topk_area <- list()
population_topk_area <- list()

for (scenario in scenarios) {
  fields <- scenario_fields(scenario)
  paths <- rank_paths(scenario)
  analogue_inputs <- c(reference_file, paths$summary)
  zone_inputs <- c(reference_file, paths$zone, paths$suitability)
  population_inputs <- c(zone_inputs, population_file)

  scenario_analogue_file <- area_cache_file("analogue", scenario)
  scenario_zone_file <- area_cache_file("source_zone", scenario)
  scenario_species_file <- area_cache_file("species", scenario)
  scenario_population_file <- area_cache_file("population", scenario)

  scenario_analogue <- read_area_cache(
    scenario_analogue_file, analogue_file, scenario, analogue_inputs,
    "minimum_analogues", length(k_values), "minimum_analogues",
    c("pixels_below", "area_below_km2")
  )
  scenario_zone <- read_area_cache(
    scenario_zone_file, zone_area_file, scenario, zone_inputs,
    "k", length(modeled_zones) * length(k_values), c("source_zone", "k"),
    c("pixels", "area_km2"), expected_zones = modeled_zones
  )
  scenario_species <- if (has_population) {
    read_area_cache(
      scenario_species_file, species_area_file, scenario, population_inputs,
      "k", length(expected_species) * length(k_values), c("Species", "k"),
      c("pixels", "area_km2"), expected_species
    )
  } else {
    data.table()
  }
  scenario_population <- if (has_population) {
    read_area_cache(
      scenario_population_file, population_area_file, scenario, population_inputs,
      "k", nrow(population) * length(k_values),
      c("Species", "source_zone", "k"), c("pixels", "area_km2"),
      expected_species, sort(unique(as.integer(population$source_zone)))
    )
  } else {
    data.table()
  }

  need_analogue <- is.null(scenario_analogue)
  need_zone <- is.null(scenario_zone)
  need_species <- has_population && is.null(scenario_species)
  need_population <- has_population && is.null(scenario_population)

  if (!any(c(need_analogue, need_zone, need_species, need_population))) {
    cat("[REUSE AREA]", scenario, "\n")
  } else {
    cat(
      "[UPDATE AREA]", scenario, "|",
      paste(c("analogue", "zone", "species", "population")[
        c(need_analogue, need_zone, need_species, need_population)
      ], collapse = ", "), "\n"
    )

    if (need_analogue) {
      rank_summary <- rast(paths$summary)
      summary_names <- names(rank_summary)
      suitable_name <- if ("n_zone_above_threshold" %in% summary_names) {
        "n_zone_above_threshold"
      } else if ("n_suitable" %in% summary_names) {
        "n_suitable"
      } else {
        stop("Rank summary lacks the number-above-threshold layer: ", paths$summary)
      }
      n_suitable <- rank_summary[[suitable_name]]
      scenario_analogue <- rbindlist(lapply(k_values, function(k) {
        scarce <- !is.na(n_suitable) & n_suitable < k
        stats <- mask_stats(scarce, get_cell_area())
        data.table(
          scenario = scenario, period = fields$period, ssp = fields$ssp,
          minimum_analogues = k, pixels_below = stats$pixels,
          area_below_km2 = stats$area_km2,
          definition = paste0(
            "number of ecotypes with dual suitability >= ", dual_threshold, " is < k"
          )
        )
      }))
      fwrite(scenario_analogue, scenario_analogue_file)
      rm(rank_summary, n_suitable)
    }

    if (need_zone || need_species) {
      ranked_zone <- rast(paths$zone)[[paste0("rank", 1:5, "_zone")]]
      ranked_suitability <- rast(paths$suitability)[[paste0("rank", 1:5, "_suit")]]
    }

    if (need_zone) {
      scenario_zone <- rbindlist(lapply(modeled_zones, function(zone) {
        first_rank <- first_suitable_rank(ranked_zone, ranked_suitability, zone)
        rank_area_rows(
          first_rank, get_cell_area(),
          list(
            scenario = scenario, period = fields$period, ssp = fields$ssp,
            source_zone = zone
          )
        )
      }))
      fwrite(scenario_zone, scenario_zone_file)
    }

    if (need_species) {
      scenario_species <- rbindlist(lapply(expected_species, function(species_name) {
        source_zones <- unique(population[Species == species_name]$source_zone)
        first_rank <- first_suitable_rank(
          ranked_zone, ranked_suitability, source_zones
        )
        rank_area_rows(
          first_rank, get_cell_area(),
          list(
            scenario = scenario, period = fields$period, ssp = fields$ssp,
            Species = species_name, n_populations = length(source_zones)
          )
        )
      }))
      fwrite(scenario_species, scenario_species_file)
    }

    if (need_population) {
      scenario_population <- merge(
        population, scenario_zone,
        by = "source_zone", allow.cartesian = TRUE
      )
      fwrite(scenario_population, scenario_population_file)
    }

    if (exists("ranked_zone", inherits = FALSE)) {
      rm(ranked_zone, ranked_suitability)
    }
    gc()
  }

  analogue_availability[[length(analogue_availability) + 1L]] <- scenario_analogue
  zone_topk_area[[length(zone_topk_area) + 1L]] <- scenario_zone
  if (has_population) {
    species_topk_area[[length(species_topk_area) + 1L]] <- scenario_species
    population_topk_area[[length(population_topk_area) + 1L]] <- scenario_population
  }
}

analogue_availability <- rbindlist(analogue_availability)
zone_topk_area <- rbindlist(zone_topk_area)
fwrite(analogue_availability, analogue_file)
fwrite(zone_topk_area, zone_area_file)

if (has_population) {
  species_topk_area <- rbindlist(species_topk_area)
  population_topk_area <- rbindlist(population_topk_area)
  fwrite(species_topk_area, species_area_file)
  fwrite(population_topk_area, population_area_file)
}


# 6. Reference zone metrics and optional multiclass comparison =================

zone_metric_results <- list()

for (k_value in k_values) {
  flow <- reference_confusion[reference_confusion[["k"]] == k_value]
  if (!nrow(flow)) next

  for (zone in modeled_zones) {
    tp <- sum(flow[observed_zone == zone & assigned_zone == zone]$pixels)
    fn <- sum(flow[observed_zone == zone & assigned_zone != zone]$pixels)
    fp <- sum(flow[observed_zone != zone & assigned_zone == zone]$pixels)
    tn <- sum(flow[observed_zone != zone & assigned_zone != zone]$pixels)

    precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_

    zone_metric_results[[length(zone_metric_results) + 1L]] <- data.table(
      k = k_value,
      zoneID = zone,
      precision = precision,
      sensitivity = sensitivity,
      specificity = specificity,
      f1 = if (is.finite(precision + sensitivity) && precision + sensitivity > 0) {
        2 * precision * sensitivity / (precision + sensitivity)
      } else NA_real_,
      tss = sensitivity + specificity - 1
    )
  }
}

reference_zone_metrics <- rbindlist(zone_metric_results)
fwrite(reference_zone_metrics, file.path(table_dir, "reference_zone_metrics.csv"))

multiclass_metric_file <- first_existing(
  c(
    file.path(base_dir, "assessment", "multiclass_rf", "reference_zone_metrics.csv"),
    file.path(base_dir, "multiclass", "tables", "reference_zone_metrics.csv"),
    file.path(base_dir, "multiclass", "tables", "zone_metrics.csv")
  ),
  "multiclass zone metrics",
  required = FALSE
)

if (!is.na(multiclass_metric_file)) {
  multiclass_metrics <- fread(multiclass_metric_file)
  zone_column <- intersect(c("zoneID", "zone"), names(multiclass_metrics))[[1]]
  setnames(multiclass_metrics, zone_column, "zoneID")
  comparison <- merge(
    reference_zone_metrics,
    multiclass_metrics[, .(zoneID, multiclass_f1 = f1)],
    by = "zoneID",
    all.x = TRUE
  )
  fwrite(comparison, file.path(table_dir, "f1_vs_multiclass.csv"))
}


# 7. Reproducibility checks and settings =======================================

novel_check <- matched_change[change_class == "novel", .(
  pixel_range = max(pixels) - min(pixels),
  area_range_km2 = max(area_km2) - min(area_km2)
), by = .(scenario)]

if (nrow(novel_check) && any(
  novel_check$pixel_range != 0 |
    novel_check$area_range_km2 > 1e-6
)) {
  stop("Novel ecosystem area changed with k; check the matched Top-k logic.")
}
fwrite(novel_check, file.path(table_dir, "novel_invariance_check.csv"))

settings <- data.table(
  setting = c(
    "workflow", "k_values", "dual_threshold", "novel_value",
    "ranking_rule", "top1_rule", "matched_rule"
  ),
  value = c(
    "selected-variable Multi-Forest",
    paste(k_values, collapse = ","),
    dual_threshold,
    novel_value,
    paste0(
      "finite dual suitability > ", rank_min_suitability,
      "; descending, zone ID breaks exact ties"
    ),
    "assigned map from script 4, including the 1e-4 tie rule",
    "same k on reference and future sides; observed reference ecotype restored when within first k ranks"
  )
)
fwrite(settings, file.path(table_dir, "settings.csv"))

cat(
  "\nCOMPLETE\n",
  "Reference Top-k rows:", nrow(reference_agreement), "\n",
  "Future matched rows:", nrow(matched_change), "\n",
  "Output:", output_root, "\n"
)
