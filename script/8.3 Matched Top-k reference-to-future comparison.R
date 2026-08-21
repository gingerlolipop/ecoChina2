# 8.3 Matched Top-k reference-to-future comparison
# ==============================================================================
# Purpose
# -------
# Compare reference-period and future ecotype-niche assignments using the same
# Top-k rule on both sides of each comparison:
#   Top-1 reference map versus Top-1 future map
#   Top-3 reference map versus Top-3 future map
#   Top-5 reference map versus Top-5 future map
#
# Top-1 is the existing assigned map and therefore retains the 1e-4 tie rule
# used in script 4.1. For Top-3 and Top-5, the observed reference-map ecotype is
# restored when it occurs among the first k dual-suitability ranks for that
# period; otherwise the period-specific Top-1 assignment is retained. Future
# cells already assigned to Zone 99 remain Zone 99.
#
# This is a rank-based sensitivity analysis. It does not refit models, rerun
# projections, recalculate dual suitability, or change the absolute definition
# of novel ecosystem: Zone 99 still requires all 53 current ecotypes to have
# dual suitability < 0.4.
#
# Outputs
# -------
# assessment_var/future_topk_matched/
#   matched_topk_change_summary.csv
#   matched_top1_consistency_check.csv
#   cache/reference_topk_maps/*.tif
#   cache/change_summaries/*.csv
# ==============================================================================

library(terra)
library(data.table)

rm(list = ls())
gc()


# 0. Paths and settings =========================================================

base_dir <- "H:/Jing/ecoChina2"

reference_file <- file.path(
  base_dir,
  "raster",
  "ecosys_ori.tif"
)

ranking_root <- file.path(
  base_dir,
  "dual suit ranking var"
)

result_map_root <- file.path(
  base_dir,
  "result maps"
)

output_dir <- file.path(
  base_dir,
  "assessment_var",
  "future_topk_matched"
)

reference_cache_dir <- file.path(
  output_dir,
  "cache",
  "reference_topk_maps"
)

summary_cache_dir <- file.path(
  output_dir,
  "cache",
  "change_summaries"
)

for (directory in c(
  output_dir,
  reference_cache_dir,
  summary_cache_dir
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

model_zoneID <- c(
  1:7,
  9:50,
  52:55
)

method_order <- c(
  "rf_var",
  "mf_var"
)

method_labels <- c(
  rf_var = "Plain RF",
  mf_var = "Plain MF RF"
)

future_order <- c(
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

topk_values <- c(
  1L,
  3L,
  5L
)

dual_threshold <- 0.4
tie_tol <- 1e-4
novel_value <- 99L

change_lookup <- data.table(
  change_code = 1:3,
  change_class = c(
    "Reference ecotype retained",
    "Different current ecotype assigned",
    "Novel ecosystem"
  )
)


# 1. Helpers ====================================================================

require_file <- function(file) {
  if (!file.exists(file)) {
    stop(
      "Missing required file: ",
      file
    )
  }
  
  file
}


assigned_map_file <- function(
    method,
    scenario) {
  file.path(
    result_map_root,
    method,
    paste0(
      "assigned_zone_",
      scenario,
      "_threshold",
      dual_threshold,
      "_tol",
      tie_tol,
      "_novel",
      novel_value,
      "_maskNA8_noNovelNormal.tif"
    )
  )
}


ranked_zone_file <- function(
    method,
    scenario) {
  file.path(
    ranking_root,
    method,
    scenario,
    "ranked_zone.tif"
  )
}


scenario_fields <- function(scenario) {
  data.table(
    period = sub(
      "SSP.*$",
      "",
      scenario
    ),
    ssp = sub(
      "^.*(SSP[0-9]+)$",
      "\\1",
      scenario
    )
  )
}


cache_is_current <- function(
    cache_files,
    input_files) {
  cache_files <- unlist(cache_files)
  input_files <- unlist(input_files)
  
  if (!all(file.exists(cache_files))) {
    return(FALSE)
  }
  
  if (!all(file.exists(input_files))) {
    return(FALSE)
  }
  
  min(file.info(cache_files)$mtime) >=
    max(file.info(input_files)$mtime)
}


check_geometry <- function(
    raster,
    template,
    label) {
  same_geometry <- compareGeom(
    raster,
    template,
    stopOnError = FALSE
  )
  
  if (!isTRUE(same_geometry)) {
    stop(
      "Raster geometry does not match for ",
      label
    )
  }
  
  invisible(TRUE)
}


reference_in_topk <- function(
    reference_zone,
    ranked_zone,
    k) {
  hit <- ifel(
    is.na(ranked_zone[[1]]),
    FALSE,
    ranked_zone[[1]] == reference_zone
  )
  
  if (k > 1L) {
    for (rank_index in 2:k) {
      this_hit <- ifel(
        is.na(ranked_zone[[rank_index]]),
        FALSE,
        ranked_zone[[rank_index]] == reference_zone
      )
      
      hit <- hit | this_hit
    }
  }
  
  hit
}


build_topk_zone_map <- function(
    reference_zone,
    assigned_zone,
    ranked_zone,
    k,
    keep_novel = FALSE) {
  if (k == 1L) {
    return(assigned_zone)
  }
  
  hit <- reference_in_topk(
    reference_zone,
    ranked_zone,
    k
  )
  
  # The observed raster also contains classes that were not among the 53
  # modeled ecotypes. Such cells cannot restore their observed class under a
  # Top-k rule, but they must remain in the same comparison domain as Top-1.
  # They therefore retain the period-specific Top-1 assignment.
  retain_reference <- ifel(
    is.na(reference_zone),
    FALSE,
    assigned_zone == reference_zone |
      hit
  )
  
  result <- ifel(
    is.na(assigned_zone),
    NA,
    ifel(
      retain_reference,
      reference_zone,
      assigned_zone
    )
  )
  
  if (keep_novel) {
    result <- ifel(
      assigned_zone == novel_value,
      novel_value,
      result
    )
  }
  
  result
}


summarize_change_map <- function(
    change_map,
    cell_area) {
  area_table <- as.data.table(
    zonal(
      cell_area,
      change_map,
      fun = "sum",
      na.rm = TRUE
    )
  )
  
  setnames(
    area_table,
    names(area_table)[1:2],
    c(
      "change_code",
      "area_km2"
    )
  )
  
  pixel_table <- as.data.table(
    freq(
      change_map,
      bylayer = FALSE
    )
  )
  
  value_column <- grep(
    "^value",
    names(pixel_table),
    value = TRUE
  )[1]
  
  if (is.na(value_column)) {
    stop(
      "Could not identify the value column returned by terra::freq()."
    )
  }
  
  count_column <- if (
    "count" %in% names(pixel_table)
  ) {
    "count"
  } else {
    names(pixel_table)[ncol(pixel_table)]
  }
  
  pixel_table <- pixel_table[
    !is.na(get(value_column)),
    .(
      change_code = as.integer(get(value_column)),
      pixel_count = as.numeric(get(count_column))
    )
  ]
  
  summary_table <- merge(
    copy(change_lookup),
    area_table,
    by = "change_code",
    all.x = TRUE
  )
  
  summary_table <- merge(
    summary_table,
    pixel_table,
    by = "change_code",
    all.x = TRUE
  )
  
  summary_table[
    is.na(area_km2),
    area_km2 := 0
  ]
  
  summary_table[
    is.na(pixel_count),
    pixel_count := 0
  ]
  
  summary_table[
    ,
    `:=`(
      total_area_km2 = sum(area_km2),
      total_pixels = sum(pixel_count)
    )
  ]
  
  summary_table[
    ,
    `:=`(
      area_share = area_km2 / total_area_km2,
      pixel_share = pixel_count / total_pixels
    )
  ]
  
  summary_table[]
}


# 2. Reference raster and matched reference-period maps =========================

require_file(reference_file)

reference_zone <- rast(reference_file)

reference_modeled <- subst(
  reference_zone,
  from = model_zoneID,
  to = model_zoneID,
  others = NA
)

# Cell area depends only on the common raster geometry and is calculated once.
cell_area_km2 <- cellSize(
  reference_modeled,
  unit = "km"
)

reference_topk_files <- list()

for (method in method_order) {
  normal_assigned_file <- require_file(
    assigned_map_file(
      method,
      "normal"
    )
  )
  
  normal_rank_file <- require_file(
    ranked_zone_file(
      method,
      "normal"
    )
  )
  
  reference_topk_files[[method]] <- list(
    top1 = normal_assigned_file
  )
  
  for (k in c(3L, 5L)) {
    cache_file <- file.path(
      reference_cache_dir,
      paste0(
        method,
        "_normal_top",
        k,
        "_common_domain_v2",
        ".tif"
      )
    )
    
    if (!cache_is_current(
      cache_file,
      c(
        reference_file,
        normal_assigned_file,
        normal_rank_file
      )
    )) {
      normal_assigned <- rast(normal_assigned_file)
      normal_ranked <- rast(normal_rank_file)[[1:k]]
      
      check_geometry(
        normal_assigned,
        reference_modeled,
        paste(method, "normal assigned map")
      )
      
      check_geometry(
        normal_ranked[[1]],
        reference_modeled,
        paste(method, "normal rank raster")
      )
      
      normal_topk <- build_topk_zone_map(
        reference_zone = reference_modeled,
        assigned_zone = normal_assigned,
        ranked_zone = normal_ranked,
        k = k,
        keep_novel = FALSE
      )
      
      writeRaster(
        normal_topk,
        cache_file,
        overwrite = TRUE,
        datatype = "INT2S",
        gdal = c(
          "COMPRESS=LZW",
          "PREDICTOR=2"
        )
      )
      
      rm(
        normal_assigned,
        normal_ranked,
        normal_topk
      )
      gc()
    }
    
    reference_topk_files[[method]][[
      paste0("top", k)
    ]] <- cache_file
  }
}


# 3. Matched reference-to-future comparisons ===================================

summary_files <- character()

for (method in method_order) {
  for (scenario in future_order) {
    future_assigned_file <- require_file(
      assigned_map_file(
        method,
        scenario
      )
    )
    
    future_rank_file <- require_file(
      ranked_zone_file(
        method,
        scenario
      )
    )
    
    future_assigned <- NULL
    future_ranked <- NULL
    
    for (k in topk_values) {
      reference_topk_file <- reference_topk_files[[method]][[
        paste0("top", k)
      ]]
      
      summary_file <- file.path(
        summary_cache_dir,
        paste0(
          method,
          "_",
          scenario,
          "_top",
          k,
          if (
            k == 1L
          ) {
            ""
          } else {
            "_common_domain_v2"
          },
          "_matched_change.csv"
        )
      )
      
      summary_files <- c(
        summary_files,
        summary_file
      )
      
      relevant_inputs <- c(
        reference_file,
        reference_topk_file,
        future_assigned_file
      )
      
      if (k > 1L) {
        relevant_inputs <- c(
          relevant_inputs,
          future_rank_file
        )
      }
      
      if (cache_is_current(
        summary_file,
        relevant_inputs
      )) {
        next
      }
      
      if (is.null(future_assigned)) {
        future_assigned <- rast(
          future_assigned_file
        )
        
        check_geometry(
          future_assigned,
          reference_modeled,
          paste(method, scenario, "assigned map")
        )
      }
      
      reference_topk <- rast(
        reference_topk_file
      )
      
      if (k == 1L) {
        future_topk <- future_assigned
      } else {
        if (is.null(future_ranked)) {
          future_ranked <- rast(
            future_rank_file
          )[[1:max(topk_values)]]
          
          check_geometry(
            future_ranked[[1]],
            reference_modeled,
            paste(method, scenario, "rank raster")
          )
        }
        
        future_topk <- build_topk_zone_map(
          reference_zone = reference_modeled,
          assigned_zone = future_assigned,
          ranked_zone = future_ranked,
          k = k,
          keep_novel = TRUE
        )
      }
      
      change_map <- ifel(
        is.na(reference_topk) |
          is.na(future_topk),
        NA,
        ifel(
          future_topk == novel_value,
          3L,
          ifel(
            future_topk == reference_topk,
            1L,
            2L
          )
        )
      )
      
      scenario_info <- scenario_fields(
        scenario
      )
      
      summary_table <- summarize_change_map(
        change_map = change_map,
        cell_area = cell_area_km2
      )
      
      summary_table[
        ,
        `:=`(
          method = method,
          method_label = unname(method_labels[method]),
          scenario = scenario,
          period = scenario_info$period,
          ssp = scenario_info$ssp,
          k = as.integer(k),
          rank_label = paste0("Top-", k),
          comparison = paste0(
            "reference Top-",
            k,
            " versus future Top-",
            k
          ),
          topk_rule = if (
            k == 1L
          ) {
            "existing assigned maps with 1e-4 tie rule"
          } else {
            paste0(
              "observed ecotype restored when within first ",
              k,
              " ranks; otherwise period-specific Top-1 retained"
            )
          },
          novel_definition = paste0(
            "Zone 99: all 53 current ecotypes have dual suitability < ",
            dual_threshold
          )
        )
      ]
      
      setcolorder(
        summary_table,
        c(
          "method",
          "method_label",
          "scenario",
          "period",
          "ssp",
          "k",
          "rank_label",
          "change_code",
          "change_class",
          "area_km2",
          "total_area_km2",
          "area_share",
          "pixel_count",
          "total_pixels",
          "pixel_share",
          "comparison",
          "topk_rule",
          "novel_definition"
        )
      )
      
      fwrite(
        summary_table,
        summary_file
      )
      
      rm(
        reference_topk,
        future_topk,
        change_map,
        summary_table
      )
      gc()
    }
    
    rm(
      future_assigned,
      future_ranked
    )
    gc()
  }
}


# 4. Combine cached summaries ===================================================

missing_summaries <- summary_files[
  !file.exists(summary_files)
]

if (length(missing_summaries) > 0L) {
  stop(
    "Missing cached summary files: ",
    paste(
      missing_summaries,
      collapse = "; "
    )
  )
}

matched_summary <- rbindlist(
  lapply(
    summary_files,
    fread
  ),
  use.names = TRUE,
  fill = TRUE
)

matched_summary[
  ,
  `:=`(
    method_order_index = match(
      method,
      method_order
    ),
    scenario_order_index = match(
      scenario,
      future_order
    )
  )
]

setorder(
  matched_summary,
  method_order_index,
  scenario_order_index,
  k,
  change_code
)

matched_summary[
  ,
  c(
    "method_order_index",
    "scenario_order_index"
  ) := NULL
]

fwrite(
  matched_summary,
  file.path(
    output_dir,
    "matched_topk_change_summary.csv"
  )
)


# 5. Top-1 regression check against the existing transition table ==============

existing_transition_file <- file.path(
  base_dir,
  "visualization var threshold0.4",
  "tables",
  "figure_transition_share_data_var.csv"
)

if (file.exists(existing_transition_file)) {
  existing_transition <- fread(
    existing_transition_file
  )
  
  existing_transition[
    ,
    change_code := fcase(
      grepl("^Reference", transition_type), 1L,
      grepl("^Different", transition_type), 2L,
      default = 3L
    )
  ]
  
  top1_check <- merge(
    matched_summary[
      k == 1L,
      .(
        method,
        scenario,
        change_code,
        new_area_km2 = area_km2,
        new_area_share = area_share
      )
    ],
    existing_transition[
      ,
      .(
        method,
        scenario,
        change_code,
        existing_area_km2 = area_km2,
        existing_area_share = area_share
      )
    ],
    by = c(
      "method",
      "scenario",
      "change_code"
    ),
    all = TRUE
  )
  
  top1_check[
    ,
    `:=`(
      relative_area_difference = abs(
        new_area_km2 - existing_area_km2
      ) / pmax(
        1,
        abs(existing_area_km2)
      ),
      absolute_share_difference = abs(
        new_area_share - existing_area_share
      )
    )
  ]
  
  fwrite(
    top1_check,
    file.path(
      output_dir,
      "matched_top1_consistency_check.csv"
    )
  )
  
  if (
    nrow(top1_check) !=
    length(method_order) *
    length(future_order) *
    nrow(change_lookup) ||
    any(!is.finite(top1_check$relative_area_difference)) ||
    max(top1_check$relative_area_difference) > 1e-8 ||
    max(top1_check$absolute_share_difference) > 1e-8
  ) {
    stop(
      "Top-1 matched results do not reproduce the existing transition table."
    )
  }
}


# 6. Final assertions ===========================================================

expected_rows <- length(method_order) *
  length(future_order) *
  length(topk_values) *
  nrow(change_lookup)

stopifnot(
  nrow(matched_summary) == expected_rows,
  !anyNA(matched_summary$area_km2),
  !anyNA(matched_summary$area_share),
  max(
    abs(
      matched_summary[
        ,
        sum(area_share),
        by = .(
          method,
          scenario,
          k
        )
      ]$V1 - 1
    )
  ) < 1e-10
)

novel_topk_check <- matched_summary[
  change_code == 3L,
  .(
    area_range_km2 = max(area_km2) -
      min(area_km2),
    relative_area_range = (
      max(area_km2) -
        min(area_km2)
    ) / pmax(
      1,
      max(area_km2)
    ),
    pixel_range = max(pixel_count) -
      min(pixel_count),
    share_range = max(area_share) -
      min(area_share)
  ),
  by = .(
    method,
    scenario
  )
]

fwrite(
  novel_topk_check,
  file.path(
    output_dir,
    "matched_topk_novel_ecosystem_check.csv"
  )
)

stopifnot(
  max(novel_topk_check$pixel_range) == 0,
  max(novel_topk_check$relative_area_range) < 1e-10,
  max(novel_topk_check$share_range) < 1e-10
)

cat(
  "Matched Top-k future comparison completed.\n",
  "Reference and future maps use the same Top-k rule.\n",
  "Novel ecosystem remains the k-independent Zone 99 definition.\n",
  "Output: ",
  file.path(
    output_dir,
    "matched_topk_change_summary.csv"
  ),
  "\n",
  sep = ""
)
