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

safe_ratio <- function(numerator, denominator) {
  if (length(numerator) != 1L || length(denominator) != 1L ||
      !is.finite(numerator) || !is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  as.numeric(numerator / denominator)
}

summary_metadata_matches <- function(table, scenario, fields) {
  required <- c("scenario", "period", "ssp")
  if (!all(required %in% names(table)) || !nrow(table)) return(FALSE)
  isTRUE(all(
    as.character(table[["scenario"]]) == as.character(scenario)[[1]]
  )) && isTRUE(all(
    as.character(table[["period"]]) == as.character(fields$period)[[1]]
  )) && isTRUE(all(
    as.character(table[["ssp"]]) == as.character(fields$ssp)[[1]]
  ))
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

matched_bundle_file <- function(scenario) {
  file.path(
    summary_cache,
    paste0("mf_var_", scenario, "_top1-5_matched_change.csv")
  )
}

table_has_all_k <- function(file, expected_k = k_values) {
  if (!file.exists(file)) return(FALSE)
  tryCatch({
    table <- fread(file, select = "k")
    setequal(sort(unique(table$k[!is.na(table$k)])), expected_k)
  }, error = function(e) FALSE)
}

reference_confusion_cache_valid <- function(file, agreement_file = NA_character_) {
  if (!file.exists(file)) return(FALSE)
  tryCatch({
    table <- fread(file)
    required <- c("k", "observed_zone", "assigned_zone", "pixels", "area_km2")
    if (!all(required %in% names(table))) return(FALSE)
    table <- table[
      k %in% k_values & observed_zone %in% modeled_zones &
        assigned_zone %in% modeled_zones
    ]
    valid <- nrow(table) > 0L && setequal(sort(unique(table$k)), k_values) &&
      all(is.finite(as.numeric(table$pixels))) &&
      all(as.numeric(table$pixels) >= 0) &&
      all(is.finite(as.numeric(table$area_km2))) &&
      all(as.numeric(table$area_km2) >= 0)
    if (!valid) return(FALSE)
    
    if (!is.na(agreement_file) && file.exists(agreement_file)) {
      agreement <- fread(agreement_file)
      required_agreement <- c("k", "total_pixels", "total_area_km2")
      if (!all(required_agreement %in% names(agreement))) return(FALSE)
      totals <- table[, .(
        confusion_pixels = sum(as.numeric(pixels)),
        confusion_area_km2 = sum(as.numeric(area_km2))
      ), by = k]
      totals <- merge(
        agreement[, .(
          k = as.integer(k),
          total_pixels = as.numeric(total_pixels),
          total_area_km2 = as.numeric(total_area_km2)
        )],
        totals, by = "k", all = TRUE
      )
      tolerance <- pmax(1e-6, abs(totals$total_area_km2) * 1e-8)
      if (nrow(totals) != length(k_values) || any(
        !is.finite(totals$total_pixels) |
        totals$total_pixels != totals$confusion_pixels |
        !is.finite(totals$total_area_km2) |
        !is.finite(totals$confusion_area_km2) |
        abs(totals$total_area_km2 - totals$confusion_area_km2) > tolerance
      )) return(FALSE)
    }
    TRUE
  }, error = function(e) FALSE)
}

normalize_matched_summary <- function(table, scenario, fields, k) {
  # Capture canonical metadata outside data.table's column scope. Legacy cache
  # files can contain an older scenario label, but their file path and the
  # current loop identify the scenario unambiguously.
  scenario_value <- as.character(scenario)[[1]]
  period_value <- as.character(fields$period)[[1]]
  ssp_value <- as.character(fields$ssp)[[1]]
  k_value <- as.integer(k)[[1]]
  table <- copy(as.data.table(table))
  if ("k" %in% names(table)) {
    input_k <- suppressWarnings(as.integer(table[["k"]]))
    if (anyNA(input_k) || any(!is.finite(input_k)) ||
        !all(input_k == k_value)) {
      stop(
        "Matched Top-k input mixes ranks or does not match requested k = ",
        k_value, "."
      )
    }
  }
  if ("scenario" %in% names(table)) {
    input_scenario <- as.character(table[["scenario"]])
    if (anyNA(input_scenario) || length(unique(input_scenario)) != 1L) {
      stop("Matched Top-k input has missing or mixed scenario labels.")
    }
  }
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
  if (any(
    !is.finite(table$code) |
    !is.finite(table$pixels) | table$pixels < 0 |
    !is.finite(table$area_km2) | table$area_km2 < 0
  )) {
    stop("Matched Top-k summary contains invalid code, pixel, or area values.")
  }
  table <- table[code %in% 1:3, .(
    pixels = sum(pixels),
    area_km2 = sum(area_km2)
  ), by = code]
  table <- merge(data.table(code = 1:3), table, by = "code", all.x = TRUE)
  table[is.na(pixels), pixels := 0]
  table[is.na(area_km2), area_km2 := 0]
  if (any((table$pixels == 0) != (table$area_km2 == 0))) {
    stop("Matched Top-k pixel and area counts disagree about an empty class.")
  }
  table[, change_class := c("stable", "changed", "novel")[as.integer(code)]]
  table[, `:=`(
    scenario = scenario_value,
    period = period_value,
    ssp = ssp_value,
    k = k_value,
    pixels = as.numeric(pixels),
    area_km2 = as.numeric(area_km2)
  )]
  table[, .(code, pixels, area_km2, change_class, scenario, period, ssp, k)]
}

matched_cache_valid <- function(file, inputs = character(),
                                expected_pixels = NA_real_,
                                expected_area_km2 = NA_real_,
                                expected_novel_pixels = NA_real_,
                                expected_novel_area_km2 = NA_real_,
                                expected_k = NA_integer_) {
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
    if (is.finite(expected_k) && "k" %in% columns) {
      file_k <- suppressWarnings(as.integer(table$k))
      if (anyNA(file_k) || any(!is.finite(file_k)) ||
          !all(file_k == expected_k)) return(FALSE)
    }
    valid <- nrow(table) == 3L && !anyDuplicated(codes) && setequal(codes, 1:3) &&
      all(is.finite(as.numeric(table[[pixel_column]]))) &&
      all(as.numeric(table[[pixel_column]]) >= 0) &&
      all(is.finite(as.numeric(table$area_km2))) &&
      all(as.numeric(table$area_km2) >= 0) &&
      all(
        (as.numeric(table[[pixel_column]]) == 0) ==
          (as.numeric(table$area_km2) == 0)
      ) && files_current(file, inputs)
    if (!valid) return(FALSE)
    
    if (is.finite(expected_pixels) &&
        sum(as.numeric(table[[pixel_column]])) != expected_pixels) return(FALSE)
    if (is.finite(expected_area_km2)) {
      tolerance <- max(1e-6, abs(expected_area_km2) * 1e-8)
      if (abs(sum(as.numeric(table$area_km2)) - expected_area_km2) > tolerance) {
        return(FALSE)
      }
    }
    novel_index <- which(codes == 3L)
    if (length(novel_index) != 1L) return(FALSE)
    if (is.finite(expected_novel_pixels) &&
        as.numeric(table[[pixel_column]][novel_index]) !=
        expected_novel_pixels) return(FALSE)
    if (is.finite(expected_novel_area_km2)) {
      novel_tolerance <- max(
        1e-6, abs(expected_novel_area_km2) * 1e-8
      )
      if (abs(as.numeric(table$area_km2[novel_index]) -
              expected_novel_area_km2) > novel_tolerance) return(FALSE)
    }
    TRUE
  }, error = function(e) FALSE)
}

assert_matched_summary <- function(table, scenario, k,
                                   expected_pixels, expected_area_km2,
                                   expected_novel_pixels,
                                   expected_novel_area_km2) {
  scenario_value <- as.character(scenario)[[1]]
  k_value <- as.integer(k)[[1]]
  required <- c(
    "code", "pixels", "area_km2", "change_class", "scenario", "k"
  )
  if (!all(required %in% names(table)) || nrow(table) != 3L ||
      anyDuplicated(table$code) || !setequal(as.integer(table$code), 1:3) ||
      any(as.character(table$scenario) != scenario_value) ||
      any(as.integer(table$k) != k_value) ||
      any(!is.finite(table$pixels) | table$pixels < 0) ||
      any(!is.finite(table$area_km2) | table$area_km2 < 0)) {
    stop(
      "Matched Top-k summary has invalid rows before writing: ",
      scenario_value, " | k = ", k_value
    )
  }
  
  pixel_total <- sum(as.numeric(table$pixels))
  area_total <- sum(as.numeric(table$area_km2))
  area_tolerance <- max(1e-6, abs(expected_area_km2) * 1e-8)
  novel_row <- table[as.integer(table$code) == 3L]
  novel_area_tolerance <- max(
    1e-6, abs(expected_novel_area_km2) * 1e-8
  )
  if (!is.finite(expected_pixels) || pixel_total != expected_pixels ||
      !is.finite(expected_area_km2) ||
      abs(area_total - expected_area_km2) > area_tolerance ||
      nrow(novel_row) != 1L ||
      !is.finite(expected_novel_pixels) ||
      as.numeric(novel_row$pixels) != expected_novel_pixels ||
      !is.finite(expected_novel_area_km2) ||
      abs(as.numeric(novel_row$area_km2) - expected_novel_area_km2) >
      novel_area_tolerance) {
    stop(
      "Matched Top-k summary failed its common-domain or novel check before ",
      "writing: ", scenario_value, " | k = ", k_value
    )
  }
  invisible(TRUE)
}

normalize_matched_bundle <- function(table, scenario, fields,
                                     expected_pixels, expected_area_km2,
                                     expected_novel_pixels,
                                     expected_novel_area_km2) {
  if (!"k" %in% names(table)) {
    stop("Matched Top-k bundle has no k column.")
  }
  result <- lapply(k_values, function(k_value) {
    # `which()` evaluates the comparison outside data.table column scope.
    row_index <- which(as.integer(table[["k"]]) == k_value)
    if (length(row_index) != 3L) {
      stop(
        "Matched Top-k bundle must contain three rows for k = ", k_value, "."
      )
    }
    one_k_input <- as.data.table(
      as.data.frame(table)[row_index, , drop = FALSE]
    )
    one_k <- normalize_matched_summary(
      one_k_input, scenario, fields, k_value
    )
    assert_matched_summary(
      one_k, scenario, k_value,
      expected_pixels, expected_area_km2,
      expected_novel_pixels, expected_novel_area_km2
    )
    one_k
  })
  rbindlist(result, use.names = TRUE, fill = TRUE)
}

# Recover the complete per-k matched table without reading rasters when a
# previous run already wrote assessment/topk/tables/matched_change.csv. The
# first clean version used terra::freq(value = TRUE), which retained the exact
# code-1 pixel count and all zonal areas but wrote zero pixels for codes 2/3.
# That one recognizable fingerprint can be repaired from the validated domain
# and novel totals; every other inconsistency fails closed and triggers the
# one-pass raster repair.
recover_matched_bundle <- function(table, scenario, fields,
                                   expected_pixels, expected_area_km2,
                                   expected_novel_pixels,
                                   expected_novel_area_km2) {
  tryCatch({
    table <- copy(as.data.table(table))
    if ("pixel_count" %in% names(table) && !"pixels" %in% names(table)) {
      setnames(table, "pixel_count", "pixels")
    }
    if ("change_code" %in% names(table) && !"code" %in% names(table)) {
      setnames(table, "change_code", "code")
    }
    if (!"code" %in% names(table) && "change_class" %in% names(table)) {
      table[, code := fcase(
        grepl("novel", change_class, ignore.case = TRUE), 3L,
        grepl("different|changed", change_class, ignore.case = TRUE), 2L,
        default = 1L
      )]
    }
    required <- c("scenario", "k", "code", "pixels", "area_km2")
    if (!all(required %in% names(table))) return(NULL)
    
    scenario_value <- as.character(scenario)[[1]]
    scenario_rows <- which(
      !is.na(table[["scenario"]]) &
        as.character(table[["scenario"]]) == scenario_value
    )
    if (!length(scenario_rows)) return(NULL)
    table <- as.data.table(
      as.data.frame(table)[scenario_rows, , drop = FALSE]
    )
    if (nrow(table) != length(k_values) * 3L) return(NULL)
    repaired_legacy <- FALSE
    result <- lapply(k_values, function(k_value) {
      row_index <- which(as.integer(table[["k"]]) == k_value)
      if (length(row_index) != 3L) stop("invalid scenario/k cardinality")
      one_k <- as.data.table(
        as.data.frame(table)[row_index, , drop = FALSE]
      )
      one_k[, `:=`(
        code = as.integer(code),
        pixels = as.numeric(pixels),
        area_km2 = as.numeric(area_km2),
        k = as.integer(k)
      )]
      if (anyNA(one_k$code) || anyDuplicated(one_k$code) ||
          !setequal(one_k$code, 1:3) ||
          any(!is.finite(one_k$area_km2) | one_k$area_km2 < 0)) {
        stop("invalid scenario/k codes or areas")
      }
      
      normalized <- NULL
      if (all(is.finite(one_k$pixels)) && all(one_k$pixels >= 0)) {
        normalized <- tryCatch({
          candidate <- normalize_matched_summary(
            one_k, scenario_value, fields, k_value
          )
          assert_matched_summary(
            candidate, scenario_value, k_value,
            expected_pixels, expected_area_km2,
            expected_novel_pixels, expected_novel_area_km2
          )
          candidate
        }, error = function(e) NULL)
      }
      if (!is.null(normalized)) return(normalized)
      
      setorder(one_k, code)
      missing_class_fingerprint <-
        is.finite(one_k$pixels[[1]]) && one_k$pixels[[1]] >= 0 &&
        one_k$pixels[[1]] == floor(one_k$pixels[[1]]) &&
        all(is.na(one_k$pixels[2:3]) | one_k$pixels[2:3] == 0)
      area_tolerance <- max(1e-6, abs(expected_area_km2) * 1e-8)
      novel_area_tolerance <- max(
        1e-6, abs(expected_novel_area_km2) * 1e-8
      )
      areas_match <-
        abs(sum(one_k$area_km2) - expected_area_km2) <= area_tolerance &&
        abs(one_k[code == 3L]$area_km2 - expected_novel_area_km2) <=
        novel_area_tolerance
      changed_pixels <- expected_pixels - one_k$pixels[[1]] -
        expected_novel_pixels
      repaired_counts_valid <-
        is.finite(expected_pixels) && expected_pixels > 0 &&
        expected_pixels == floor(expected_pixels) &&
        is.finite(expected_novel_pixels) && expected_novel_pixels >= 0 &&
        expected_novel_pixels == floor(expected_novel_pixels) &&
        is.finite(changed_pixels) && changed_pixels >= 0 &&
        changed_pixels == floor(changed_pixels)
      if (!missing_class_fingerprint || !areas_match ||
          !repaired_counts_valid) {
        stop("not a safely recoverable legacy freq(value = TRUE) table")
      }
      one_k[code == 2L, pixels := changed_pixels]
      one_k[code == 3L, pixels := expected_novel_pixels]
      candidate <- normalize_matched_summary(
        one_k, scenario_value, fields, k_value
      )
      assert_matched_summary(
        candidate, scenario_value, k_value,
        expected_pixels, expected_area_km2,
        expected_novel_pixels, expected_novel_area_km2
      )
      repaired_legacy <<- TRUE
      candidate
    })
    result <- rbindlist(result, use.names = TRUE, fill = TRUE)
    attr(result, "legacy_repaired") <- repaired_legacy
    result
  }, error = function(e) NULL)
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
  # `value` in terra::freq() is a numeric value to count. Using value = TRUE
  # counts only raster value 1; bylayer = FALSE returns all codes and counts.
  pixel <- as.data.table(freq(code_raster, bylayer = FALSE))
  if (!nrow(pixel)) return(data.table())
  value_column <- grep("^value$|^value_", names(pixel), ignore.case = TRUE, value = TRUE)[1]
  count_column <- if ("count" %in% names(pixel)) {
    "count"
  } else {
    tail(names(pixel), 1)
  }
  if (is.na(value_column) || is.na(count_column)) {
    stop("Could not identify value/count columns returned by terra::freq().")
  }
  pixel <- pixel[, .(
    code = as.integer(get(value_column)),
    pixels = as.numeric(get(count_column))
  )]
  # terra::freq() can include the NA category outside a masked analysis domain.
  # It must not enter a confusion table because comparisons with it propagate NA.
  pixel <- pixel[is.finite(code) & is.finite(pixels) & pixels >= 0]
  
  area <- as.data.table(zonal(cell_area, code_raster, fun = "sum", na.rm = TRUE))
  if (!nrow(area)) return(data.table())
  setnames(area, names(area)[1:2], c("code", "area_km2"))
  area[, `:=`(code = as.integer(code), area_km2 = as.numeric(area_km2))]
  area <- area[is.finite(code) & is.finite(area_km2) & area_km2 >= 0]
  
  result <- merge(pixel, area, by = "code", all = TRUE)
  if (!nrow(result) || anyDuplicated(result$code) || any(
    !is.finite(result$code) |
    !is.finite(result$pixels) | result$pixels < 0 |
    result$pixels != floor(result$pixels) |
    !is.finite(result$area_km2) | result$area_km2 < 0
  )) {
    stop(
      "terra::freq() and terra::zonal() returned inconsistent coded-raster summaries."
    )
  }
  result[]
}

clean_reference_confusion <- function(table) {
  required <- c("k", "observed_zone", "assigned_zone", "pixels", "area_km2")
  if (!all(required %in% names(table))) {
    stop("Reference confusion table lacks: ",
         paste(setdiff(required, names(table)), collapse = ", "))
  }
  
  table <- copy(table)
  table[, `:=`(
    k = as.integer(k),
    observed_zone = as.integer(observed_zone),
    assigned_zone = as.integer(assigned_zone),
    pixels = as.numeric(pixels),
    area_km2 = as.numeric(area_km2)
  )]
  
  valid_key <- table$k %in% k_values &
    table$observed_zone %in% modeled_zones &
    table$assigned_zone %in% modeled_zones
  removed <- sum(!valid_key)
  table <- table[valid_key]
  
  if (!nrow(table) || any(
    !is.finite(table$pixels) | table$pixels < 0 |
    !is.finite(table$area_km2) | table$area_km2 < 0
  )) {
    stop("Reference confusion table contains invalid in-domain counts.")
  }
  
  table <- table[, .(
    pixels = sum(pixels),
    area_km2 = sum(area_km2)
  ), by = .(k, observed_zone, assigned_zone)]
  
  list(table = table, removed = removed)
}

clean_future_retention <- function(table) {
  required <- c(
    "scenario", "k", "class_code", "retention_class", "pixels", "area_km2",
    "total_pixels", "total_area_km2", "novel_pixels", "novel_area_km2"
  )
  if (!all(required %in% names(table))) {
    stop("Future-retention table lacks: ",
         paste(setdiff(required, names(table)), collapse = ", "))
  }
  valid_labels <- c(
    "stable_top1", "former_raw_rank1_after_tie", "former_rank2",
    "former_rank3", "former_rank4", "former_rank5", "novel",
    "below_top5", "cumulative_topk"
  )
  table <- copy(table)[retention_class %in% valid_labels]
  table[, `:=`(
    k = as.integer(k),
    class_code = as.integer(class_code),
    pixels = as.numeric(pixels),
    area_km2 = as.numeric(area_km2),
    total_pixels = as.numeric(total_pixels),
    total_area_km2 = as.numeric(total_area_km2),
    novel_pixels = as.numeric(novel_pixels),
    novel_area_km2 = as.numeric(novel_area_km2)
  )]
  if (any(!is.finite(table$pixels) | table$pixels < 0 |
          !is.finite(table$area_km2) | table$area_km2 < 0)) {
    stop("Future-retention table contains invalid counts.")
  }
  
  cumulative <- table[retention_class == "cumulative_topk"]
  classes <- table[retention_class != "cumulative_topk"]
  if (uniqueN(table$scenario) != 1L ||
      any(!(classes$class_code %in% 1:8)) ||
      anyDuplicated(classes[, .(scenario, class_code)])) {
    stop("Future-retention class keys are invalid or duplicated.")
  }
  class_totals <- merge(
    data.table(class_code = 1:8),
    classes[, .(class_code, class_pixels = pixels, class_area_km2 = area_km2)],
    by = "class_code", all.x = TRUE
  )
  class_totals[is.na(class_pixels), class_pixels := 0]
  class_totals[is.na(class_area_km2), class_area_km2 := 0]
  expected_retained <- rbindlist(lapply(k_values, function(k_value) {
    retained_codes <- if (k_value == 1L) 1L else c(1L, 2L:(k_value + 1L))
    data.table(
      k = k_value,
      expected_retained_pixels = sum(
        class_totals[class_code %in% retained_codes]$class_pixels
      ),
      expected_retained_area_km2 = sum(
        class_totals[class_code %in% retained_codes]$class_area_km2
      )
    )
  }))
  cumulative_check <- merge(cumulative, expected_retained, by = "k", all = TRUE)
  retained_area_tolerance <- pmax(
    1e-6, abs(cumulative_check$expected_retained_area_km2) * 1e-8
  )
  novel_class_pixels <- class_totals[class_code == 7L]$class_pixels
  novel_class_area <- class_totals[class_code == 7L]$class_area_km2
  if (nrow(cumulative) != length(k_values) ||
      anyDuplicated(cumulative[, .(scenario, k)]) ||
      !setequal(cumulative$k, k_values) ||
      any(
        !is.finite(cumulative$total_pixels) | cumulative$total_pixels <= 0 |
        !is.finite(cumulative$total_area_km2) | cumulative$total_area_km2 <= 0 |
        !is.finite(cumulative$novel_pixels) | cumulative$novel_pixels < 0 |
        !is.finite(cumulative$novel_area_km2) | cumulative$novel_area_km2 < 0
      ) ||
      uniqueN(cumulative$total_pixels) != 1L ||
      uniqueN(cumulative$total_area_km2) != 1L ||
      sum(classes$pixels) != cumulative$total_pixels[[1]] ||
      abs(sum(classes$area_km2) - cumulative$total_area_km2[[1]]) >
      max(1e-6, cumulative$total_area_km2[[1]] * 1e-8) ||
      any(cumulative_check$pixels != cumulative_check$expected_retained_pixels) ||
      any(abs(
        cumulative_check$area_km2 -
        cumulative_check$expected_retained_area_km2
      ) > retained_area_tolerance) ||
      any(cumulative$novel_pixels != novel_class_pixels) ||
      any(abs(cumulative$novel_area_km2 - novel_class_area) >
          max(1e-6, abs(novel_class_area) * 1e-8)) ||
      any(diff(cumulative[order(k)]$pixels) < 0) ||
      any(diff(cumulative[order(k)]$area_km2) < -1e-6)) {
    stop("Future-retention classes do not form one complete common domain.")
  }
  
  table[retention_class == "cumulative_topk", `:=`(
    retained_pixel_share = mapply(safe_ratio, pixels, total_pixels),
    retained_area_share = mapply(safe_ratio, area_km2, total_area_km2)
  )]
  table[]
}

clean_future_retention_by_zone <- function(
    table, expected_pixels = NA_real_, expected_area_km2 = NA_real_) {
  required <- c(
    "scenario", "k", "normal_zone", "status_code", "status", "pixels", "area_km2"
  )
  if (!all(required %in% names(table))) {
    stop("Zone-retention table lacks: ",
         paste(setdiff(required, names(table)), collapse = ", "))
  }
  table <- copy(table)
  table[, `:=`(
    k = as.integer(k),
    normal_zone = as.integer(normal_zone),
    status_code = as.integer(status_code),
    pixels = as.numeric(pixels),
    area_km2 = as.numeric(area_km2)
  )]
  table <- table[
    k %in% k_values & normal_zone %in% modeled_zones & status_code %in% 1:3
  ]
  table[, status := c("retained", "not_retained", "novel")[status_code]]
  if (any(!is.finite(table$pixels) | table$pixels < 0 |
          !is.finite(table$area_km2) | table$area_km2 < 0) ||
      uniqueN(table$scenario) != 1L ||
      anyDuplicated(table[, .(scenario, k, normal_zone, status_code)]) ||
      !setequal(table$k, k_values)) {
    stop("Zone-retention table contains invalid or duplicated counts.")
  }
  totals <- table[, .(
    pixels = sum(pixels),
    area_km2 = sum(area_km2)
  ), by = k]
  if (is.finite(expected_pixels) && any(totals$pixels != expected_pixels)) {
    stop("Zone-retention pixel totals do not match the common domain.")
  }
  if (is.finite(expected_area_km2)) {
    tolerance <- max(1e-6, abs(expected_area_km2) * 1e-8)
    if (any(abs(totals$area_km2 - expected_area_km2) > tolerance)) {
      stop("Zone-retention area totals do not match the common domain.")
    }
  }
  table[]
}

future_retention_cache_valid <- function(file, inputs = character()) {
  if (!file.exists(file) || !files_current(file, inputs)) return(FALSE)
  tryCatch({
    clean_future_retention(fread(file))
    TRUE
  }, error = function(e) FALSE)
}

future_retention_by_zone_cache_valid <- function(
    file, inputs = character(), expected_pixels = NA_real_,
    expected_area_km2 = NA_real_) {
  if (!file.exists(file) || !files_current(file, inputs)) return(FALSE)
  tryCatch({
    clean_future_retention_by_zone(
      fread(file), expected_pixels, expected_area_km2
    )
    TRUE
  }, error = function(e) FALSE)
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

# Recover reference agreement and confusion in one sequential read. The former
# implementation traversed the national raster several times for every k, which
# made recovery from an invalid summary cache unnecessarily expensive.
reference_summary_repair_rows <- function(reference, reference_topk,
                                          cell_area) {
  if (nlyr(reference_topk) != length(k_values) ||
      !compareGeom(reference, reference_topk, stopOnError = FALSE) ||
      !compareGeom(reference, cell_area, stopOnError = FALSE)) {
    stop("Reference-summary repair needs matching five-layer Top-k geometry.")
  }
  
  repair_stack <- c(reference, reference_topk, cell_area)
  n_zones <- length(modeled_zones)
  n_flows <- n_zones * n_zones
  flow_pixels <- matrix(0, nrow = n_flows, ncol = length(k_values))
  flow_area <- flow_pixels
  matched_pixels <- numeric(length(k_values))
  matched_area <- numeric(length(k_values))
  total_pixels <- 0
  total_area <- 0
  
  rows_per_block <- max(1L, floor(250000 / ncol(repair_stack)))
  starts <- seq.int(1L, nrow(repair_stack), by = rows_per_block)
  report_every <- max(1L, ceiling(length(starts) / 10L))
  reading <- FALSE
  readStart(repair_stack)
  reading <- TRUE
  on.exit({
    if (reading) try(readStop(repair_stack), silent = TRUE)
  }, add = TRUE)
  
  for (block_index in seq_along(starts)) {
    start_row <- starts[[block_index]]
    n_rows <- min(rows_per_block, nrow(repair_stack) - start_row + 1L)
    values <- as.matrix(readValues(
      repair_stack, row = start_row, nrows = n_rows, mat = TRUE
    ))
    reference_values <- as.numeric(values[, 1L])
    topk_values <- values[, 1L + seq_along(k_values), drop = FALSE]
    area_values <- as.numeric(values[, length(k_values) + 2L])
    # Match the established common domain exactly: a modeled reference class
    # with a non-NA Top-1 assignment. Infinite or non-modeled assignments are
    # errors, not cells to silently drop from the denominator.
    valid <- reference_values %in% modeled_zones & !is.na(topk_values[, 1L])
    
    if (any(valid & rowSums(
      !is.finite(topk_values) | !(topk_values %in% modeled_zones)
    ) > 0L)) {
      stop("Reference Top-k maps contain an invalid value in the common domain.")
    }
    for (rank_index in 2:length(k_values)) {
      impossible_value <- valid &
        topk_values[, rank_index] != reference_values &
        topk_values[, rank_index] != topk_values[, 1L]
      reversed_recovery <- valid &
        topk_values[, rank_index - 1L] == reference_values &
        topk_values[, rank_index] != reference_values
      if (any(impossible_value) || any(reversed_recovery)) {
        stop(
          "Reference Top-k maps are stale or violate cumulative recovery at k = ",
          rank_index, "."
        )
      }
    }
    if (any(valid & (!is.finite(area_values) | area_values < 0))) {
      stop("Cell-area values are invalid inside the reference domain.")
    }
    
    safe_area <- area_values
    safe_area[!is.finite(safe_area)] <- 0
    observed_index <- match(reference_values, modeled_zones)
    total_pixels <- total_pixels + sum(valid)
    total_area <- total_area + sum(safe_area[valid])
    
    for (rank_index in seq_along(k_values)) {
      assigned_values <- as.integer(topk_values[, rank_index])
      assigned_index <- match(assigned_values, modeled_zones)
      flow_code <- ifelse(
        valid,
        (observed_index - 1L) * n_zones + assigned_index,
        0L
      )
      flow_pixels[, rank_index] <- flow_pixels[, rank_index] +
        tabulate(flow_code, nbins = n_flows)
      flow_area[, rank_index] <- flow_area[, rank_index] +
        weighted_tabulate(flow_code, safe_area, n_flows)
      matched <- valid & assigned_values == reference_values
      matched_pixels[[rank_index]] <-
        matched_pixels[[rank_index]] + sum(matched)
      matched_area[[rank_index]] <-
        matched_area[[rank_index]] + sum(safe_area[matched])
    }
    
    if (block_index %% report_every == 0L || block_index == length(starts)) {
      cat(
        "[REFERENCE REPAIR]",
        sprintf("%d%%", round(100 * block_index / length(starts))), "\n"
      )
    }
  }
  
  readStop(repair_stack)
  reading <- FALSE
  flow_pixel_totals <- colSums(flow_pixels)
  flow_area_totals <- colSums(flow_area)
  flow_area_tolerance <- max(1e-6, abs(total_area) * 1e-8)
  if (!is.finite(total_pixels) || total_pixels <= 0 ||
      !is.finite(total_area) || total_area <= 0 ||
      any(flow_pixel_totals != total_pixels) ||
      any(abs(flow_area_totals - total_area) > flow_area_tolerance) ||
      any(diff(matched_pixels) < 0) ||
      any(diff(matched_area) < -flow_area_tolerance)) {
    stop("One-pass reference summary failed its domain or cumulative checks.")
  }
  agreement <- data.table(
    k = k_values,
    total_pixels = total_pixels,
    matched_pixels = matched_pixels,
    pixel_agreement = vapply(matched_pixels, function(value) {
      safe_ratio(value, total_pixels)
    }, numeric(1)),
    total_area_km2 = total_area,
    matched_area_km2 = matched_area,
    area_agreement = vapply(matched_area, function(value) {
      safe_ratio(value, total_area)
    }, numeric(1))
  )
  confusion <- rbindlist(lapply(seq_along(k_values), function(rank_index) {
    present <- which(
      flow_pixels[, rank_index] > 0 | flow_area[, rank_index] > 0
    )
    data.table(
      code = as.integer(
        modeled_zones[((present - 1L) %/% n_zones) + 1L] * 1000L +
          modeled_zones[((present - 1L) %% n_zones) + 1L]
      ),
      pixels = flow_pixels[present, rank_index],
      area_km2 = flow_area[present, rank_index],
      k = k_values[[rank_index]],
      observed_zone = modeled_zones[((present - 1L) %/% n_zones) + 1L],
      assigned_zone = modeled_zones[((present - 1L) %% n_zones) + 1L]
    )
  }), use.names = TRUE, fill = TRUE)
  
  list(agreement = agreement, confusion = confusion)
}

# Repair counts and areas from existing rank rasters in one sequential read, so
# recovery from the old terra::freq(value = TRUE) bug never falls back to one
# national zonal calculation for every zone and species.
rank_area_repair_rows <- function(ranked_zone, ranked_suitability, cell_area,
                                  species_zones = list(),
                                  species_area_names = names(species_zones),
                                  label = "") {
  if (nlyr(ranked_zone) != length(k_values) ||
      nlyr(ranked_suitability) != length(k_values) ||
      !compareGeom(ranked_zone, ranked_suitability, stopOnError = FALSE) ||
      !compareGeom(ranked_zone, cell_area, stopOnError = FALSE)) {
    stop("Area repair needs matching five-layer ranks and cell-area geometry.")
  }
  
  rank_stack <- c(ranked_zone, ranked_suitability, cell_area)
  zone_exact <- matrix(
    0, nrow = length(modeled_zones), ncol = length(k_values),
    dimnames = list(as.character(modeled_zones), as.character(k_values))
  )
  zone_area_exact <- zone_exact
  species_names <- names(species_zones)
  normalized_species_zones <- lapply(species_zones, function(zones) {
    sort(unique(as.integer(zones[zones %in% modeled_zones])))
  })
  species_signature <- if (length(species_names)) {
    vapply(normalized_species_zones, paste, collapse = ",", character(1))
  } else {
    character()
  }
  group_signatures <- unique(species_signature)
  species_group <- match(species_signature, group_signatures)
  species_area_names <- intersect(
    as.character(species_area_names), species_names
  )
  area_group_indices <- unique(
    species_group[species_names %in% species_area_names]
  )
  group_zones <- lapply(group_signatures, function(signature) {
    normalized_species_zones[[match(signature, species_signature)]]
  })
  group_exact <- matrix(
    0, nrow = length(group_signatures), ncol = length(k_values),
    dimnames = list(group_signatures, as.character(k_values))
  )
  group_area_exact <- group_exact
  membership <- matrix(
    FALSE, nrow = max(modeled_zones) + 1L, ncol = length(group_signatures)
  )
  if (length(group_signatures)) {
    for (group_index in seq_along(group_signatures)) {
      membership[group_zones[[group_index]] + 1L, group_index] <- TRUE
    }
  }
  group_chunks <- if (length(group_signatures)) {
    split(
      seq_along(group_signatures),
      ceiling(seq_along(group_signatures) / 16L)
    )
  } else {
    list()
  }
  
  # About 250,000 cells x 11 layers keeps the temporary matrices modest. Species
  # that share the same source-zone set are tallied together in small chunks.
  rows_per_block <- max(1L, floor(250000 / ncol(rank_stack)))
  starts <- seq.int(1L, nrow(rank_stack), by = rows_per_block)
  report_every <- max(1L, ceiling(length(starts) / 10L))
  
  reading <- FALSE
  readStart(rank_stack)
  reading <- TRUE
  on.exit({
    if (reading) try(readStop(rank_stack), silent = TRUE)
  }, add = TRUE)
  
  for (block_index in seq_along(starts)) {
    start_row <- starts[[block_index]]
    n_rows <- min(rows_per_block, nrow(rank_stack) - start_row + 1L)
    values <- readValues(
      rank_stack, row = start_row, nrows = n_rows, mat = TRUE
    )
    values <- as.matrix(values)
    zone_values <- values[, seq_along(k_values), drop = FALSE]
    suitability_values <- values[
      , length(k_values) + seq_along(k_values), drop = FALSE
    ]
    area_values <- as.numeric(values[, 2L * length(k_values) + 1L])
    valid_rank <- is.finite(zone_values) & is.finite(suitability_values) &
      suitability_values >= dual_threshold
    if (any(rowSums(valid_rank) > 0L &
            (!is.finite(area_values) | area_values < 0))) {
      stop("Cell-area values are invalid inside the ranked analysis domain.")
    }
    safe_area_values <- area_values
    safe_area_values[!is.finite(safe_area_values)] <- 0
    
    for (rank_index in seq_along(k_values)) {
      rank_cells <- which(valid_rank[, rank_index])
      zone_index <- match(
        as.integer(zone_values[rank_cells, rank_index]),
        modeled_zones
      )
      keep_zone <- !is.na(zone_index)
      zone_index <- zone_index[keep_zone]
      rank_cells <- rank_cells[keep_zone]
      zone_exact[, rank_index] <- zone_exact[, rank_index] +
        tabulate(zone_index, nbins = length(modeled_zones))
      if (length(zone_index)) {
        area_by_zone <- rowsum(
          area_values[rank_cells], zone_index, reorder = FALSE
        )
        area_zone_index <- as.integer(rownames(area_by_zone))
        zone_area_exact[area_zone_index, rank_index] <-
          zone_area_exact[area_zone_index, rank_index] + area_by_zone[, 1]
      }
    }
    
    if (length(group_signatures)) {
      zone_rows <- matrix(as.integer(zone_values) + 1L, nrow = nrow(zone_values))
      zone_rows[
        is.na(zone_rows) | zone_rows < 1L | zone_rows > nrow(membership)
      ] <- 1L
      for (group_chunk in group_chunks) {
        area_positions <- which(group_chunk %in% area_group_indices)
        first_rank <- matrix(
          0L, nrow = nrow(zone_values), ncol = length(group_chunk)
        )
        for (rank_index in seq_along(k_values)) {
          membership_hit <- membership[
            zone_rows[, rank_index], group_chunk, drop = FALSE
          ]
          membership_hit <- membership_hit & valid_rank[, rank_index]
          hit <- first_rank == 0L & membership_hit
          first_rank[hit] <- rank_index
        }
        for (rank_index in seq_along(k_values)) {
          group_exact[group_chunk, rank_index] <-
            group_exact[group_chunk, rank_index] +
            colSums(first_rank == rank_index)
          if (length(area_positions)) {
            area_groups <- group_chunk[area_positions]
            group_area_exact[area_groups, rank_index] <-
              group_area_exact[area_groups, rank_index] +
              colSums(
                (first_rank[, area_positions, drop = FALSE] == rank_index) *
                  safe_area_values
              )
          }
        }
      }
    }
    
    if (block_index %% report_every == 0L || block_index == length(starts)) {
      cat(
        "[ONE-PASS REPAIR]", label, "|",
        sprintf("%d%%", round(100 * block_index / length(starts))), "\n"
      )
    }
  }
  
  readStop(rank_stack)
  reading <- FALSE
  if (length(group_signatures)) {
    skipped_area_groups <- setdiff(
      seq_along(group_signatures), area_group_indices
    )
    if (length(skipped_area_groups)) {
      group_area_exact[skipped_area_groups, ] <- NA_real_
    }
  }
  
  zone_rows <- rbindlist(lapply(seq_along(modeled_zones), function(index) {
    data.table(
      source_zone = modeled_zones[[index]],
      k = k_values,
      pixels = as.numeric(cumsum(zone_exact[index, ])),
      area_km2 = as.numeric(cumsum(zone_area_exact[index, ]))
    )
  }))
  species_rows <- if (length(species_names)) {
    rbindlist(lapply(seq_along(species_names), function(index) {
      data.table(
        Species = species_names[[index]],
        k = k_values,
        pixels = as.numeric(cumsum(group_exact[species_group[[index]], ])),
        area_km2 = as.numeric(cumsum(
          group_area_exact[species_group[[index]], ]
        ))
      )
    }))
  } else {
    data.table()
  }
  
  list(zone = zone_rows, species = species_rows)
}

weighted_tabulate <- function(codes, weights, nbins) {
  result <- numeric(nbins)
  keep <- is.finite(codes) & codes >= 1L & codes <= nbins &
    is.finite(weights) & weights >= 0
  if (!any(keep)) return(result)
  grouped <- rowsum(
    weights[keep], as.integer(codes[keep]), reorder = FALSE
  )
  result[as.integer(rownames(grouped))] <- grouped[, 1]
  result
}

# Rebuild all corrupted future summary tables in a single sequential read. This
# replaces dozens of separate national freq/zonal/global passes while leaving
# the cached assigned, rank, and Top-k rasters untouched.
future_summary_repair_rows <- function(
    reference, normal_assigned, normal_ranked,
    future_assigned, future_ranked, cell_area, scenario, fields) {
  rasters <- list(
    reference, normal_assigned, normal_ranked,
    future_assigned, future_ranked, cell_area
  )
  if (nlyr(normal_ranked) != length(k_values) ||
      nlyr(future_ranked) != length(k_values) ||
      !all(vapply(rasters[-1], function(x) {
        compareGeom(reference, x, stopOnError = FALSE)
      }, logical(1)))) {
    stop("One-pass future-summary repair received incompatible raster geometry.")
  }
  
  repair_stack <- c(
    reference, normal_assigned, normal_ranked,
    future_assigned, future_ranked, cell_area
  )
  class_pixels <- numeric(8L)
  class_area <- numeric(8L)
  retained_pixels <- numeric(length(k_values))
  retained_area <- numeric(length(k_values))
  status_nbins <- max(modeled_zones) * 10L + 3L
  status_pixels <- matrix(
    0, nrow = status_nbins, ncol = length(k_values)
  )
  status_area <- status_pixels
  matched_pixels <- matrix(0, nrow = 3L, ncol = length(k_values))
  matched_area <- matched_pixels
  total_pixels <- 0
  total_area <- 0
  novel_pixels <- 0
  novel_area <- 0
  
  rows_per_block <- max(1L, floor(250000 / ncol(repair_stack)))
  starts <- seq.int(1L, nrow(repair_stack), by = rows_per_block)
  report_every <- max(1L, ceiling(length(starts) / 10L))
  reading <- FALSE
  readStart(repair_stack)
  reading <- TRUE
  on.exit({
    if (reading) try(readStop(repair_stack), silent = TRUE)
  }, add = TRUE)
  
  for (block_index in seq_along(starts)) {
    start_row <- starts[[block_index]]
    n_rows <- min(rows_per_block, nrow(repair_stack) - start_row + 1L)
    values <- as.matrix(readValues(
      repair_stack, row = start_row, nrows = n_rows, mat = TRUE
    ))
    
    reference_values <- as.integer(values[, 1L])
    normal_values <- as.integer(values[, 2L])
    normal_rank_values <- values[, 3L:7L, drop = FALSE]
    future_values <- as.integer(values[, 8L])
    future_rank_values <- values[, 9L:13L, drop = FALSE]
    area_values <- as.numeric(values[, 14L])
    
    valid <- reference_values %in% modeled_zones &
      is.finite(normal_values) & is.finite(future_values)
    if (any(valid & !(normal_values %in% modeled_zones))) {
      stop("Normal assigned map has a non-modeled value in the reference domain.")
    }
    if (any(valid & !(future_values %in% c(modeled_zones, novel_value)))) {
      stop("Future assigned map has an unexpected value in the reference domain.")
    }
    if (any(valid & (!is.finite(area_values) | area_values < 0))) {
      stop("Cell-area values are invalid inside the future comparison domain.")
    }
    safe_area <- area_values
    safe_area[!is.finite(safe_area)] <- 0
    
    novel <- valid & future_values == novel_value
    stable <- valid & !novel & future_values == normal_values
    first_normal_rank <- integer(length(reference_values))
    normal_reference_hit <- rep(FALSE, length(reference_values))
    future_reference_hit <- rep(FALSE, length(reference_values))
    
    for (rank_index in seq_along(k_values)) {
      hit_normal <- valid & first_normal_rank == 0L &
        is.finite(future_rank_values[, rank_index]) &
        future_rank_values[, rank_index] == normal_values
      first_normal_rank[hit_normal] <- rank_index
    }
    
    retention_code <- integer(length(reference_values))
    retention_code[valid] <- 8L
    has_normal_rank <- valid & first_normal_rank > 0L
    retention_code[has_normal_rank] <- first_normal_rank[has_normal_rank] + 1L
    retention_code[stable] <- 1L
    retention_code[novel] <- 7L
    class_pixels <- class_pixels + tabulate(retention_code, nbins = 8L)
    class_area <- class_area + weighted_tabulate(
      retention_code, safe_area, 8L
    )
    
    total_pixels <- total_pixels + sum(valid)
    total_area <- total_area + sum(safe_area[valid])
    novel_pixels <- novel_pixels + sum(novel)
    novel_area <- novel_area + sum(safe_area[novel])
    
    for (rank_index in seq_along(k_values)) {
      retention_hit <- if (rank_index == 1L) {
        stable
      } else {
        stable | (first_normal_rank > 0L & first_normal_rank <= rank_index)
      }
      retained <- valid & !novel & retention_hit
      retained_pixels[[rank_index]] <-
        retained_pixels[[rank_index]] + sum(retained)
      retained_area[[rank_index]] <-
        retained_area[[rank_index]] + sum(safe_area[retained])
      
      status <- integer(length(reference_values))
      status[valid] <- 2L
      status[retained] <- 1L
      status[novel] <- 3L
      status_code <- ifelse(valid, normal_values * 10L + status, 0L)
      status_pixels[, rank_index] <- status_pixels[, rank_index] +
        tabulate(status_code, nbins = status_nbins)
      status_area[, rank_index] <- status_area[, rank_index] +
        weighted_tabulate(status_code, safe_area, status_nbins)
      
      normal_reference_hit <- normal_reference_hit |
        (
          valid & is.finite(normal_rank_values[, rank_index]) &
            normal_rank_values[, rank_index] == reference_values
        )
      future_reference_hit <- future_reference_hit |
        (
          valid & is.finite(future_rank_values[, rank_index]) &
            future_rank_values[, rank_index] == reference_values
        )
      if (rank_index == 1L) {
        reference_topk_values <- normal_values
        future_topk_values <- future_values
      } else {
        reference_topk_values <- ifelse(
          normal_reference_hit, reference_values, normal_values
        )
        future_topk_values <- ifelse(
          novel,
          novel_value,
          ifelse(future_reference_hit, reference_values, future_values)
        )
      }
      change_code <- integer(length(reference_values))
      change_code[valid] <- 2L
      change_code[valid & future_topk_values == reference_topk_values] <- 1L
      change_code[novel] <- 3L
      matched_pixels[, rank_index] <- matched_pixels[, rank_index] +
        tabulate(change_code, nbins = 3L)
      matched_area[, rank_index] <- matched_area[, rank_index] +
        weighted_tabulate(change_code, safe_area, 3L)
    }
    
    if (block_index %% report_every == 0L || block_index == length(starts)) {
      cat(
        "[SUMMARY REPAIR]", scenario, "|",
        sprintf("%d%%", round(100 * block_index / length(starts))), "\n"
      )
    }
  }
  
  readStop(repair_stack)
  reading <- FALSE
  
  matched_pixel_totals <- colSums(matched_pixels)
  matched_area_totals <- colSums(matched_area)
  matched_area_tolerance <- max(1e-6, abs(total_area) * 1e-8)
  if (!is.finite(total_pixels) || total_pixels <= 0 ||
      !is.finite(total_area) || total_area <= 0 ||
      any(matched_pixel_totals != total_pixels) ||
      any(abs(matched_area_totals - total_area) > matched_area_tolerance) ||
      any(matched_pixels[3L, ] != novel_pixels) ||
      any(abs(matched_area[3L, ] - novel_area) > matched_area_tolerance)) {
    stop(
      "One-pass future summary failed to conserve its common-domain totals: ",
      scenario
    )
  }
  
  class_labels <- data.table(
    class_code = 1:8,
    retention_class = c(
      "stable_top1", "former_raw_rank1_after_tie", "former_rank2",
      "former_rank3", "former_rank4", "former_rank5", "novel", "below_top5"
    )
  )
  class_table <- copy(class_labels)
  class_table[, `:=`(
    pixels = class_pixels,
    area_km2 = class_area,
    scenario = scenario,
    period = fields$period,
    ssp = fields$ssp
  )]
  cumulative_table <- data.table(
    class_code = NA_integer_,
    retention_class = "cumulative_topk",
    scenario = scenario,
    period = fields$period,
    ssp = fields$ssp,
    k = k_values,
    pixels = retained_pixels,
    area_km2 = retained_area,
    total_pixels = total_pixels,
    total_area_km2 = total_area,
    retained_pixel_share = vapply(retained_pixels, function(value) {
      safe_ratio(value, total_pixels)
    }, numeric(1)),
    retained_area_share = vapply(retained_area, function(value) {
      safe_ratio(value, total_area)
    }, numeric(1)),
    novel_pixels = novel_pixels,
    novel_area_km2 = novel_area
  )
  retention_table <- rbindlist(
    list(class_table, cumulative_table), use.names = TRUE, fill = TRUE
  )
  
  zone_tables <- lapply(seq_along(k_values), function(rank_index) {
    present <- which(
      status_pixels[, rank_index] > 0 | status_area[, rank_index] > 0
    )
    data.table(
      code = present,
      pixels = status_pixels[present, rank_index],
      area_km2 = status_area[present, rank_index],
      normal_zone = as.integer(present %/% 10L),
      status_code = as.integer(present %% 10L),
      status = c("retained", "not_retained", "novel")[
        as.integer(present %% 10L)
      ],
      scenario = scenario,
      period = fields$period,
      ssp = fields$ssp,
      k = k_values[[rank_index]]
    )
  })
  retention_by_zone_table <- rbindlist(zone_tables, use.names = TRUE, fill = TRUE)
  
  matched_table <- rbindlist(lapply(seq_along(k_values), function(rank_index) {
    data.table(
      code = 1:3,
      pixels = matched_pixels[, rank_index],
      area_km2 = matched_area[, rank_index],
      change_class = c("stable", "changed", "novel"),
      scenario = scenario,
      period = fields$period,
      ssp = fields$ssp,
      k = k_values[[rank_index]]
    )
  }))
  
  list(
    retention = retention_table,
    retention_by_zone = retention_by_zone_table,
    matched = matched_table
  )
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
    reference_confusion_cache_valid(
      reference_confusion_file, reference_agreement_file
    ) &&
    files_current(
      c(reference_agreement_file, reference_confusion_file),
      c(reference_file, unlist(reference_topk_files))
    )) {
  reference_agreement <- fread(reference_agreement_file)
  reference_confusion <- fread(reference_confusion_file)
  cat("[REUSE TABLE] reference Top-k agreement and confusion\n")
} else {
  reference_topk <- rast(unlist(
    reference_topk_files[as.character(k_values)], use.names = FALSE
  ))
  repaired_reference <- reference_summary_repair_rows(
    reference, reference_topk, get_cell_area()
  )
  reference_agreement <- repaired_reference$agreement
  reference_confusion <- repaired_reference$confusion
  rm(reference_topk, repaired_reference)
}

cleaned_reference_confusion <- clean_reference_confusion(reference_confusion)
reference_confusion <- cleaned_reference_confusion$table
if (cleaned_reference_confusion$removed > 0L) {
  cat(
    "[CLEAN TABLE] removed", cleaned_reference_confusion$removed,
    "out-of-domain reference-confusion rows\n"
  )
}

reference_totals <- reference_confusion[, .(
  confusion_pixels = sum(pixels),
  confusion_area_km2 = sum(area_km2)
), by = k]
reference_expected <- reference_agreement[, .(
  k = as.integer(k),
  total_pixels = as.numeric(total_pixels),
  total_area_km2 = as.numeric(total_area_km2)
)]
reference_totals <- merge(reference_expected, reference_totals, by = "k", all = TRUE)
area_tolerance <- pmax(1e-6, abs(reference_totals$total_area_km2) * 1e-8)
if (nrow(reference_totals) != length(k_values) ||
    !setequal(reference_totals$k, k_values) ||
    any(
      !is.finite(reference_totals$total_pixels) |
      !is.finite(reference_totals$confusion_pixels) |
      reference_totals$total_pixels != reference_totals$confusion_pixels |
      !is.finite(reference_totals$total_area_km2) |
      !is.finite(reference_totals$confusion_area_km2) |
      abs(
        reference_totals$total_area_km2 -
        reference_totals$confusion_area_km2
      ) > area_tolerance
    )) {
  stop("Reference confusion totals do not match reference agreement totals.")
}
fwrite(reference_agreement, reference_agreement_file)
fwrite(reference_confusion, reference_confusion_file)


# 4. Future rank retention and matched Top-k comparison ========================

future_retention <- list()
future_retention_by_zone <- list()
matched_change <- list()
scenario_domain_totals <- list()
matched_master_cache_file <- file.path(table_dir, "matched_change.csv")
matched_master_cache <- if (file.exists(matched_master_cache_file)) {
  tryCatch(fread(matched_master_cache_file), error = function(e) NULL)
} else {
  NULL
}

for (scenario in future_scenarios) {
  cat("\n[FUTURE TOP-K]", scenario, "\n")
  scenario_value <- as.character(scenario)[[1]]
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
  retention_inputs <- c(
    reference_file, normal_assigned_file, future_assigned_file, future_rank_file
  )
  
  retention_ok <- future_retention_cache_valid(
    retention_cache_file, retention_inputs
  )
  cached_retention <- if (retention_ok) {
    clean_future_retention(fread(retention_cache_file))
  } else {
    NULL
  }
  cached_totals <- if (!is.null(cached_retention)) {
    cached_retention[
      retention_class == "cumulative_topk",
      .(
        expected_pixels = total_pixels[[1]],
        expected_area_km2 = total_area_km2[[1]]
      )
    ][1]
  } else {
    data.table(expected_pixels = NA_real_, expected_area_km2 = NA_real_)
  }
  retention_zone_ok <- retention_ok && future_retention_by_zone_cache_valid(
    retention_zone_cache_file, retention_inputs,
    cached_totals$expected_pixels, cached_totals$expected_area_km2
  )
  scenario_summary_repair <- NULL
  
  if (retention_ok && retention_zone_ok) {
    scenario_retention <- cached_retention
    scenario_retention_by_zone <- clean_future_retention_by_zone(
      fread(retention_zone_cache_file),
      cached_totals$expected_pixels, cached_totals$expected_area_km2
    )
    cat("[REUSE TABLE] future retention |", scenario, "\n")
  } else {
    scenario_summary_repair <- future_summary_repair_rows(
      reference, normal_assigned, normal_ranked,
      future_assigned, future_ranked, get_cell_area(), scenario, fields
    )
    scenario_retention <- scenario_summary_repair$retention
    scenario_retention_by_zone <- scenario_summary_repair$retention_by_zone
    scenario_retention <- clean_future_retention(scenario_retention)
    scenario_totals <- scenario_retention[
      retention_class == "cumulative_topk",
      .(
        expected_pixels = total_pixels[[1]],
        expected_area_km2 = total_area_km2[[1]]
      )
    ][1]
    scenario_retention_by_zone <- clean_future_retention_by_zone(
      scenario_retention_by_zone,
      scenario_totals$expected_pixels, scenario_totals$expected_area_km2
    )
  }
  
  # Old result sets used more than one scenario naming convention. Keep their
  # validated numeric summaries, but canonicalize labels from the current loop
  # before joining scenarios. This is a metadata repair and reads no raster.
  write_retention_cache <- !is.null(scenario_summary_repair) ||
    !summary_metadata_matches(scenario_retention, scenario_value, fields) ||
    !summary_metadata_matches(
      scenario_retention_by_zone, scenario_value, fields
    )
  scenario_retention[, `:=`(
    scenario = scenario_value,
    period = as.character(fields$period)[[1]],
    ssp = as.character(fields$ssp)[[1]]
  )]
  scenario_retention_by_zone[, `:=`(
    scenario = scenario_value,
    period = as.character(fields$period)[[1]],
    ssp = as.character(fields$ssp)[[1]]
  )]
  if (write_retention_cache) {
    fwrite(scenario_retention, retention_cache_file)
    fwrite(scenario_retention_by_zone, retention_zone_cache_file)
  }
  
  scenario_totals <- scenario_retention[
    retention_class == "cumulative_topk",
    .(
      expected_pixels = total_pixels[[1]],
      expected_area_km2 = total_area_km2[[1]],
      expected_novel_pixels = novel_pixels[[1]],
      expected_novel_area_km2 = novel_area_km2[[1]]
    )
  ][1]
  
  # Treat retention + five matched tables as one cache bundle. If an older
  # retention table passed its internal checks but was built on a different
  # common domain, its matched files will fail these totals. Repair that one
  # scenario once here, rather than running five separate national summaries
  # and discovering the mismatch only after all scenarios finish.
  summary_bundle_inputs <- c(
    reference_file, normal_assigned_file, future_assigned_file,
    future_rank_file, unname(unlist(reference_topk_files))
  )
  scenario_bundle_file <- matched_bundle_file(scenario_value)
  scenario_matched_bundle <- NULL
  matched_bundle_ready <- vapply(k_values, function(k_value) {
    current_file <- matched_summary_file(scenario_value, k_value)
    legacy_file <- matched_summary_file(scenario_value, k_value, legacy = TRUE)
    current_inputs <- c(
      reference_file,
      normal_assigned_file,
      future_assigned_file,
      if (k_value > 1L) future_rank_file else character(),
      reference_topk_files[[as.character(k_value)]]
    )
    matched_cache_valid(
      current_file, current_inputs,
      scenario_totals$expected_pixels, scenario_totals$expected_area_km2,
      scenario_totals$expected_novel_pixels,
      scenario_totals$expected_novel_area_km2, k_value
    ) || matched_cache_valid(
      legacy_file, current_inputs,
      scenario_totals$expected_pixels, scenario_totals$expected_area_km2,
      scenario_totals$expected_novel_pixels,
      scenario_totals$expected_novel_area_km2, k_value
    )
  }, logical(1))
  
  # If a retention repair already traversed the rasters, checkpoint its
  # correct 15-row matched result immediately, before any per-k file handling.
  force_bundle_matched <- !is.null(scenario_summary_repair)
  if (!is.null(scenario_summary_repair)) {
    scenario_matched_bundle <- normalize_matched_bundle(
      scenario_summary_repair$matched, scenario_value, fields,
      scenario_totals$expected_pixels,
      scenario_totals$expected_area_km2,
      scenario_totals$expected_novel_pixels,
      scenario_totals$expected_novel_area_km2
    )
    fwrite(scenario_matched_bundle, scenario_bundle_file)
  }
  
  # The user's earlier run reached Section 5 and therefore left a global
  # matched_change.csv written before the later 5x cache corruption. Recover
  # from that table (or a newer all-k checkpoint) first; only a fully validated
  # scenario x k bundle is accepted. This path performs no raster summary.
  if (!all(matched_bundle_ready) && is.null(scenario_summary_repair)) {
    candidate_bundles <- list()
    if (file.exists(scenario_bundle_file) &&
        files_current(scenario_bundle_file, summary_bundle_inputs)) {
      candidate_bundles[["all-k checkpoint"]] <- tryCatch(
        fread(scenario_bundle_file), error = function(e) NULL
      )
    }
    if (!is.null(matched_master_cache) &&
        files_current(matched_master_cache_file, summary_bundle_inputs)) {
      candidate_bundles[["global matched table"]] <- matched_master_cache
    }
    candidate_bundles <- candidate_bundles[
      !vapply(candidate_bundles, is.null, logical(1))
    ]
    if (length(candidate_bundles)) {
      for (candidate_name in names(candidate_bundles)) {
        candidate <- recover_matched_bundle(
          candidate_bundles[[candidate_name]], scenario_value, fields,
          scenario_totals$expected_pixels,
          scenario_totals$expected_area_km2,
          scenario_totals$expected_novel_pixels,
          scenario_totals$expected_novel_area_km2
        )
        if (!is.null(candidate)) {
          scenario_matched_bundle <- candidate
          fwrite(scenario_matched_bundle, scenario_bundle_file)
          for (recovered_k in k_values) {
            recovered_index <- which(
              as.integer(scenario_matched_bundle[["k"]]) == recovered_k
            )
            fwrite(
              scenario_matched_bundle[recovered_index],
              matched_summary_file(scenario_value, recovered_k)
            )
          }
          matched_bundle_ready[] <- TRUE
          repair_note <- if (isTRUE(attr(candidate, "legacy_repaired"))) {
            " (legacy pixel counts reconstructed)"
          } else {
            ""
          }
          cat(
            "[RECOVER MATCHED WITHOUT RASTER]", scenario_value, "|",
            candidate_name, repair_note, "\n"
          )
          break
        }
      }
    }
  }
  
  if (!all(matched_bundle_ready) && is.null(scenario_summary_repair)) {
    cat("[REPAIR SUMMARY BUNDLE]", scenario_value, "\n")
    scenario_summary_repair <- future_summary_repair_rows(
      reference, normal_assigned, normal_ranked,
      future_assigned, future_ranked, get_cell_area(), scenario_value, fields
    )
    scenario_retention <- clean_future_retention(
      scenario_summary_repair$retention
    )
    scenario_retention[, `:=`(
      scenario = scenario_value,
      period = as.character(fields$period)[[1]],
      ssp = as.character(fields$ssp)[[1]]
    )]
    scenario_totals <- scenario_retention[
      retention_class == "cumulative_topk",
      .(
        expected_pixels = total_pixels[[1]],
        expected_area_km2 = total_area_km2[[1]],
        expected_novel_pixels = novel_pixels[[1]],
        expected_novel_area_km2 = novel_area_km2[[1]]
      )
    ][1]
    scenario_retention_by_zone <- clean_future_retention_by_zone(
      scenario_summary_repair$retention_by_zone,
      scenario_totals$expected_pixels, scenario_totals$expected_area_km2
    )
    scenario_retention_by_zone[, `:=`(
      scenario = scenario_value,
      period = as.character(fields$period)[[1]],
      ssp = as.character(fields$ssp)[[1]]
    )]
    fwrite(scenario_retention, retention_cache_file)
    fwrite(scenario_retention_by_zone, retention_zone_cache_file)
    scenario_matched_bundle <- normalize_matched_bundle(
      scenario_summary_repair$matched, scenario_value, fields,
      scenario_totals$expected_pixels,
      scenario_totals$expected_area_km2,
      scenario_totals$expected_novel_pixels,
      scenario_totals$expected_novel_area_km2
    )
    fwrite(scenario_matched_bundle, scenario_bundle_file)
    force_bundle_matched <- TRUE
  }
  
  scenario_domain_totals[[scenario_value]] <- data.table(
    scenario = scenario_value,
    expected_pixels = as.numeric(scenario_totals$expected_pixels),
    expected_area_km2 = as.numeric(scenario_totals$expected_area_km2)
  )
  future_retention[[length(future_retention) + 1L]] <- scenario_retention
  future_retention_by_zone[[length(future_retention_by_zone) + 1L]] <-
    scenario_retention_by_zone
  
  for (k in k_values) {
    # Do not use a bare `k` as the right-hand side of a data.table row filter:
    # it can resolve to the table's `k` column and turn the filter into k == k.
    k_value <- as.integer(k)[[1]]
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
    # Prefer the canonical cache just written by this script. Legacy results are
    # a fallback only, so a repaired table is never shadowed on the next run.
    cache_file <- if (force_bundle_matched) {
      NA_character_
    } else if (matched_cache_valid(
      output_file, matched_inputs,
      scenario_totals$expected_pixels, scenario_totals$expected_area_km2,
      scenario_totals$expected_novel_pixels,
      scenario_totals$expected_novel_area_km2, k_value
    )) {
      output_file
    } else if (matched_cache_valid(
      legacy_file, matched_inputs,
      scenario_totals$expected_pixels, scenario_totals$expected_area_km2,
      scenario_totals$expected_novel_pixels,
      scenario_totals$expected_novel_area_km2, k_value
    )) {
      legacy_file
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
      assert_matched_summary(
        change_table, scenario_value, k_value,
        scenario_totals$expected_pixels,
        scenario_totals$expected_area_km2,
        scenario_totals$expected_novel_pixels,
        scenario_totals$expected_novel_area_km2
      )
      matched_change[[length(matched_change) + 1L]] <- change_table
      cat("[REUSE MATCHED]", scenario, "| k =", k, "\n")
      next
    }
    
    if (k == 1L) {
      # Top-1 is already the assigned map from script 4. Keep it in memory for
      # the matched summary instead of writing a duplicate national raster.
      if (is.null(valid)) {
        valid <- get_valid_reference() & !is.na(future_assigned)
      }
      future_topk <- ifel(valid, future_assigned, NA)
    } else if (future_topk_current) {
      future_topk <- rast(future_topk_file)
    } else {
      if (is.null(valid)) {
        valid <- get_valid_reference() & !is.na(future_assigned)
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
    } else if (!is.null(scenario_summary_repair)) {
      repair_index <- which(
        as.integer(scenario_summary_repair$matched[["k"]]) == k_value
      )
      if (length(repair_index) != 3L) {
        stop(
          "Future matched repair must contain exactly three change classes: ",
          scenario_value, " | k = ", k_value
        )
      }
      change_table <- normalize_matched_summary(
        scenario_summary_repair$matched[repair_index],
        scenario, fields, k_value
      )
      cat("[REPAIR MATCHED]", scenario, "| k =", k, "\n")
    } else {
      if (is.null(valid)) {
        valid <- get_valid_reference() & !is.na(future_assigned)
      }
      reference_topk <- rast(reference_topk_files[[as.character(k)]])
      change_class_raster <- ifel(
        valid,
        ifel(
          future_topk == novel_value,
          3L,
          ifel(future_topk == reference_topk, 1L, 2L)
        ),
        NA
      )
      change_table <- coded_area(change_class_raster, get_cell_area())
      change_table <- normalize_matched_summary(change_table, scenario, fields, k)
      rm(reference_topk, change_class_raster)
    }
    assert_matched_summary(
      change_table, scenario_value, k_value,
      scenario_totals$expected_pixels,
      scenario_totals$expected_area_km2,
      scenario_totals$expected_novel_pixels,
      scenario_totals$expected_novel_area_km2
    )
    if (is.na(cache_file)) fwrite(change_table, output_file)
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

future_retention <- rbindlist(future_retention, use.names = TRUE, fill = TRUE)
future_retention_by_zone <- rbindlist(
  future_retention_by_zone, use.names = TRUE, fill = TRUE
)
matched_change <- rbindlist(matched_change, use.names = TRUE, fill = TRUE)

expected_matched_rows <- length(future_scenarios) * length(k_values) * 3L
if (nrow(matched_change) != expected_matched_rows ||
    !setequal(unique(as.character(matched_change$scenario)), future_scenarios) ||
    !setequal(unique(as.integer(matched_change$k)), k_values) ||
    anyDuplicated(matched_change[, .(scenario, k, code)]) ||
    any(
      !is.finite(matched_change$pixels) | matched_change$pixels < 0 |
      !is.finite(matched_change$area_km2) | matched_change$area_km2 < 0
    )) {
  stop("Matched Top-k summary is incomplete or contains duplicate scenario/k/code rows.")
}

matched_totals <- matched_change[, .(
  total_pixels = sum(pixels),
  total_area_km2 = sum(area_km2)
), by = .(scenario, k)]
expected_totals <- rbindlist(
  scenario_domain_totals, use.names = TRUE, fill = TRUE
)
if (nrow(expected_totals) != length(future_scenarios) ||
    anyDuplicated(expected_totals$scenario) ||
    !setequal(expected_totals$scenario, future_scenarios)) {
  stop("Scenario-domain totals are incomplete or duplicated.")
}
matched_totals <- merge(matched_totals, expected_totals, by = "scenario", all.x = TRUE)
matched_area_tolerance <- pmax(
  1e-6, abs(matched_totals$expected_area_km2) * 1e-8
)
matched_totals[, `:=`(
  pixel_delta = total_pixels - expected_pixels,
  area_delta_km2 = total_area_km2 - expected_area_km2
)]
bad_matched_totals <- matched_totals[
  !is.finite(matched_totals$expected_pixels) |
    matched_totals$total_pixels != matched_totals$expected_pixels |
    !is.finite(matched_totals$expected_area_km2) |
    abs(matched_totals$total_area_km2 - matched_totals$expected_area_km2) >
    matched_area_tolerance
]
if (nrow(bad_matched_totals)) {
  print(bad_matched_totals[, .(
    scenario, k, total_pixels, expected_pixels, pixel_delta,
    total_area_km2, expected_area_km2, area_delta_km2
  )])
  stop(
    "Matched Top-k totals do not match the common reference/future domain; ",
    "the offending scenario/k rows are printed above."
  )
}

matched_novel <- matched_change[code == 3L, .(
  scenario, k, matched_novel_pixels = pixels,
  matched_novel_area_km2 = area_km2
)]
expected_novel <- future_retention[
  retention_class == "cumulative_topk",
  .(
    scenario, k,
    expected_novel_pixels = novel_pixels,
    expected_novel_area_km2 = novel_area_km2
  )
]
novel_consistency <- merge(
  expected_novel, matched_novel, by = c("scenario", "k"), all = TRUE
)
novel_area_tolerance <- pmax(
  1e-6, abs(novel_consistency$expected_novel_area_km2) * 1e-8
)
if (nrow(novel_consistency) != length(future_scenarios) * length(k_values) ||
    any(
      !is.finite(novel_consistency$expected_novel_pixels) |
      !is.finite(novel_consistency$matched_novel_pixels) |
      novel_consistency$expected_novel_pixels !=
      novel_consistency$matched_novel_pixels |
      !is.finite(novel_consistency$expected_novel_area_km2) |
      !is.finite(novel_consistency$matched_novel_area_km2) |
      abs(
        novel_consistency$expected_novel_area_km2 -
        novel_consistency$matched_novel_area_km2
      ) > novel_area_tolerance
    )) {
  stop("Matched novel counts do not match the independently summarized novel mask.")
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
  population[, source_zone := as.integer(source_zone)]
  if (anyNA(population[, .(Species, source_zone)]) ||
      any(!population$source_zone %in% modeled_zones)) {
    stop("Eligible population rows need non-missing Species and modeled source_zone.")
  }
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
  repair_zone_cache <- FALSE
  if (is.null(scenario_zone)) {
    # Even when the old bug left some numeric cells as NA, its row/key structure
    # is reusable. Repair every count and area in one sequential rank scan.
    scenario_zone <- read_area_cache(
      scenario_zone_file, zone_area_file, scenario, zone_inputs,
      "k", length(modeled_zones) * length(k_values), c("source_zone", "k"),
      character(), expected_zones = modeled_zones
    )
    repair_zone_cache <- !is.null(scenario_zone)
  }
  scenario_species <- if (has_population) {
    read_area_cache(
      scenario_species_file, species_area_file, scenario, population_inputs,
      "k", length(expected_species) * length(k_values), c("Species", "k"),
      c("pixels", "area_km2"), expected_species
    )
  } else {
    data.table()
  }
  repair_species_cache <- FALSE
  if (has_population && is.null(scenario_species)) {
    scenario_species <- read_area_cache(
      scenario_species_file, species_area_file, scenario, population_inputs,
      "k", length(expected_species) * length(k_values), c("Species", "k"),
      character(), expected_species
    )
    repair_species_cache <- !is.null(scenario_species)
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
  need_population <- has_population && (
    is.null(scenario_population) || need_zone || repair_zone_cache
  )
  
  if (!any(c(
    need_analogue, need_zone, need_species, need_population,
    repair_zone_cache, repair_species_cache
  ))) {
    cat("[REUSE AREA]", scenario, "\n")
  } else {
    update_labels <- c(
      if (need_analogue) "analogue",
      if (need_zone) "zone area",
      if (repair_zone_cache) "zone counts/area repair",
      if (need_species) "species area",
      if (repair_species_cache) "species counts/area repair",
      if (need_population) "population table"
    )
    cat(
      "[UPDATE AREA]", scenario, "|",
      paste(update_labels, collapse = ", "), "\n"
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
    
    if (need_zone || need_species || repair_zone_cache || repair_species_cache) {
      ranked_zone <- rast(paths$zone)[[paste0("rank", 1:5, "_zone")]]
      ranked_suitability <- rast(paths$suitability)[[paste0("rank", 1:5, "_suit")]]
    }
    
    if (need_zone || need_species || repair_zone_cache || repair_species_cache) {
      species_zone_groups <- if (need_species || repair_species_cache) {
        split(
          as.integer(population$source_zone),
          as.character(population$Species)
        )
      } else {
        list()
      }
      species_area_names <- if (need_species) {
        expected_species
      } else if (repair_species_cache) {
        if (!"area_km2" %in% names(scenario_species)) {
          expected_species
        } else {
          unique(as.character(scenario_species[
            !is.finite(as.numeric(area_km2)) | as.numeric(area_km2) < 0,
            Species
          ]))
        }
      } else {
        character()
      }
      repaired_area <- rank_area_repair_rows(
        ranked_zone, ranked_suitability, get_cell_area(),
        species_zones = species_zone_groups,
        species_area_names = species_area_names,
        label = scenario
      )
      
      if (need_zone) {
        scenario_zone <- copy(repaired_area$zone)
        scenario_zone[, `:=`(
          scenario = scenario, period = fields$period, ssp = fields$ssp
        )]
        fwrite(scenario_zone, scenario_zone_file)
      } else if (repair_zone_cache) {
        scenario_zone[
          repaired_area$zone,
          on = .(source_zone, k),
          `:=`(pixels = i.pixels, area_km2 = i.area_km2)
        ]
        if (!complete_area_rows(
          scenario_zone, "k", length(modeled_zones) * length(k_values),
          c("source_zone", "k"), c("pixels", "area_km2"),
          expected_zones = modeled_zones
        )) {
          stop("One-pass repair failed for source-zone area: ", scenario)
        }
        fwrite(scenario_zone, scenario_zone_file)
      }
      
      if (need_species) {
        n_population_lookup <- population[, .(
          n_populations = uniqueN(source_zone)
        ), by = Species]
        scenario_species <- merge(
          repaired_area$species, n_population_lookup,
          by = "Species", all.x = TRUE
        )
        scenario_species[, `:=`(
          scenario = scenario, period = fields$period, ssp = fields$ssp
        )]
        fwrite(scenario_species, scenario_species_file)
      } else if (repair_species_cache) {
        scenario_species[
          repaired_area$species,
          on = .(Species, k),
          pixels := i.pixels
        ]
        if (length(species_area_names)) {
          scenario_species[
            repaired_area$species[Species %in% species_area_names],
            on = .(Species, k),
            area_km2 := i.area_km2
          ]
        }
        if (!complete_area_rows(
          scenario_species, "k", length(expected_species) * length(k_values),
          c("Species", "k"), c("pixels", "area_km2"), expected_species
        )) {
          stop("One-pass repair failed for species area: ", scenario)
        }
        fwrite(scenario_species, scenario_species_file)
      }
      rm(repaired_area, species_zone_groups, species_area_names)
      if (exists("n_population_lookup", inherits = FALSE)) {
        rm(n_population_lookup)
      }
    }
    
    if (need_population) {
      scenario_population <- merge(
        population, scenario_zone,
        by = "source_zone", allow.cartesian = TRUE
      )
      fwrite(scenario_population, scenario_population_file)
    }
    
    if (!complete_area_rows(
      scenario_zone, "k", length(modeled_zones) * length(k_values),
      c("source_zone", "k"), c("pixels", "area_km2"),
      expected_zones = modeled_zones
    )) {
      stop("Source-zone area checkpoint is incomplete after update: ", scenario)
    }
    if (has_population && !complete_area_rows(
      scenario_species, "k", length(expected_species) * length(k_values),
      c("Species", "k"), c("pixels", "area_km2"), expected_species
    )) {
      stop("Species area checkpoint is incomplete after update: ", scenario)
    }
    if (has_population && !complete_area_rows(
      scenario_population, "k", nrow(population) * length(k_values),
      c("Species", "source_zone", "k"), c("pixels", "area_km2"),
      expected_species, sort(unique(as.integer(population$source_zone)))
    )) {
      stop("Population area checkpoint is incomplete after update: ", scenario)
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

analogue_availability <- rbindlist(
  analogue_availability, use.names = TRUE, fill = TRUE
)
zone_topk_area <- rbindlist(zone_topk_area, use.names = TRUE, fill = TRUE)
fwrite(analogue_availability, analogue_file)
fwrite(zone_topk_area, zone_area_file)

if (has_population) {
  species_topk_area <- rbindlist(
    species_topk_area, use.names = TRUE, fill = TRUE
  )
  population_topk_area <- rbindlist(
    population_topk_area, use.names = TRUE, fill = TRUE
  )
  fwrite(species_topk_area, species_area_file)
  fwrite(population_topk_area, population_area_file)
}
cat("[CHECKPOINT] Top-k area tables saved; later failures will reuse them.\n")


# 6. Reference zone metrics and optional multiclass comparison =================

zone_metric_results <- list()

for (k_value in k_values) {
  flow <- reference_confusion[reference_confusion[["k"]] == k_value]
  if (!nrow(flow)) next
  if (any(
    !is.finite(flow$pixels) | flow$pixels < 0 |
    !(flow$observed_zone %in% modeled_zones) |
    !(flow$assigned_zone %in% modeled_zones)
  )) {
    stop("Invalid reference-confusion rows remain at k = ", k_value, ".")
  }
  
  for (zone in modeled_zones) {
    tp <- sum(flow[observed_zone == zone & assigned_zone == zone]$pixels)
    fn <- sum(flow[observed_zone == zone & assigned_zone != zone]$pixels)
    fp <- sum(flow[observed_zone != zone & assigned_zone == zone]$pixels)
    tn <- sum(flow[observed_zone != zone & assigned_zone != zone]$pixels)
    
    precision <- safe_ratio(tp, tp + fp)
    sensitivity <- safe_ratio(tp, tp + fn)
    specificity <- safe_ratio(tn, tn + fp)
    
    zone_metric_results[[length(zone_metric_results) + 1L]] <- data.table(
      k = k_value,
      zoneID = zone,
      precision = precision,
      sensitivity = sensitivity,
      specificity = specificity,
      f1 = safe_ratio(2 * tp, 2 * tp + fp + fn),
      tss = if (all(is.finite(c(sensitivity, specificity)))) {
        sensitivity + specificity - 1
      } else {
        NA_real_
      }
    )
  }
}

reference_zone_metrics <- rbindlist(zone_metric_results)
if (nrow(reference_zone_metrics) != length(k_values) * length(modeled_zones) ||
    anyDuplicated(reference_zone_metrics[, .(k, zoneID)])) {
  stop("Reference zone metrics are incomplete or duplicated.")
}
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
  tryCatch({
    multiclass_metrics <- fread(multiclass_metric_file)
    zone_candidates <- intersect(c("zoneID", "zone"), names(multiclass_metrics))
    if (!length(zone_candidates) || !"f1" %in% names(multiclass_metrics)) {
      stop("table needs zoneID/zone and f1 columns")
    }
    zone_column <- zone_candidates[[1]]
    if (zone_column != "zoneID") {
      setnames(multiclass_metrics, zone_column, "zoneID")
    }
    multiclass_metrics[, zoneID := as.integer(zoneID)]
    multiclass_metrics[, f1 := as.numeric(f1)]
    if (nrow(multiclass_metrics) != length(modeled_zones) ||
        anyNA(multiclass_metrics$zoneID) ||
        anyDuplicated(multiclass_metrics$zoneID) ||
        !setequal(multiclass_metrics$zoneID, modeled_zones) ||
        any(!is.finite(multiclass_metrics$f1))) {
      stop("table is incomplete or has invalid zone/F1 values")
    }
    comparison <- merge(
      reference_zone_metrics,
      multiclass_metrics[, .(zoneID, multiclass_f1 = f1)],
      by = "zoneID",
      all.x = TRUE
    )
    fwrite(comparison, file.path(table_dir, "f1_vs_multiclass.csv"))
  }, error = function(e) {
    warning(
      "[SKIP OPTIONAL] multiclass comparison was not written: ",
      conditionMessage(e),
      ". Script 7 may still be writing ", multiclass_metric_file,
      call. = FALSE
    )
  })
}


# 7. Reproducibility checks and settings =======================================

novel_check <- matched_change[change_class == "novel", .(
  pixel_range = max(pixels) - min(pixels),
  area_range_km2 = max(area_km2) - min(area_km2),
  area_tolerance_km2 = max(1e-6, max(abs(area_km2)) * 1e-8)
), by = .(scenario)]

if (nrow(novel_check) && (any(!is.finite(novel_check$pixel_range)) ||
                          any(!is.finite(novel_check$area_range_km2)) || any(
                            novel_check$pixel_range != 0 |
                            novel_check$area_range_km2 > 2 * novel_check$area_tolerance_km2
                          ))) {
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
