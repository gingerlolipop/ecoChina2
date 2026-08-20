# 11.4 Top-k downstream figures
# ==============================================================================
# Run after scripts 6.1, 8.1, 8.2 and 11.2.
#
# This script adds rank-aware downstream figures without refitting any model.
# It uses the saved ranked-zone and ranked-suitability rasters and preserves the
# dual-suitability threshold of 0.4.
#
# Definitions
# -----------
# Figure 6b:
#   Top-k non-retention is the area where the workflow's normal-period Top-1
#   ecotype is not retained within the future first k positive-suitability
#   ranks. Novel cells are included. This is a relative-rank diagnostic, not
#   novel niche space or a threshold-based suitable-area estimate, and must
#   decrease or remain unchanged from Top-1 to Top-3 to Top-5.
#
# Supplementary analogue scarcity:
#   Scarcity at minimum count k is the area with fewer than k current ecotype
#   analogues reaching dual suitability 0.4. The k = 1 case is exactly novel
#   niche space. Because the minimum required candidate count increases with k,
#   scarcity must increase or remain unchanged from k = 1 to 3 to 5.
#
# Figure 8:
#   A species is represented within Top-k when at least one of its source
#   ecotypes satisfies the Top-k diagnostic criterion at a cell. Each cell is
#   counted once per species even when several source ecotypes are represented.
#
# Figure 10a:
#   A population proxy is represented within Top-k when its source ecotype
#   satisfies the Top-k diagnostic criterion at a cell.
#
# Species/population Top-k criterion
# ----------------------------------
#   Top-1 uses the existing assigned map, preserving the 1e-4 tie rule.
#   Top-3 and Top-5 include Top-1 plus source ecotypes occurring within the
#   first k saved ranks with dual suitability >= 0.4. This matches the nested
#   diagnostic logic used for reference agreement while retaining the absolute
#   suitability threshold required for species/population niche area.
#
# Regression guards verify that Top-1 areas reproduce the existing novel,
# assigned-species and assigned-population results from script 6.1.
# ==============================================================================

library(terra)
library(data.table)
library(ggplot2)

rm(list = ls())
gc()


# 0. Paths and settings =========================================================

base_dir <- "H:/Jing/ecoChina2"

assessment_dir <- file.path(
  base_dir,
  "assessment_var"
)

assigned_result_root <- file.path(
  base_dir,
  "future tree niche var"
)

assigned_table_dir <- file.path(
  assigned_result_root,
  "tables"
)

ranking_root <- file.path(
  base_dir,
  "dual suit ranking var"
)

result_map_root <- file.path(
  base_dir,
  "result maps"
)

reference_file <- file.path(
  base_dir,
  "raster/ecosys_ori.tif"
)

population_lookup_file <- file.path(
  assigned_table_dir,
  "population_projection_lookup_var.csv"
)

future_topk_retention_file <- file.path(
  assessment_dir,
  "future_topk_analysis",
  "future_topk_retention_overall.csv"
)

output_root <- file.path(
  base_dir,
  "visualization var threshold0.4"
)

figure_dir <- file.path(
  output_root,
  "figures"
)

table_dir <- file.path(
  output_root,
  "tables"
)

population_detail_dir <- file.path(
  figure_dir,
  "population detail"
)

cache_dir <- file.path(
  table_dir,
  "topk downstream cache"
)

for (directory in c(
  figure_dir,
  table_dir,
  population_detail_dir,
  cache_dir
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

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

period_levels <- c(
  "2011-2040",
  "2041-2070",
  "2071-2100"
)

ssp_levels <- c(
  "SSP245",
  "SSP585"
)

rank_cutoffs <- c(
  1L,
  3L,
  5L
)

rank_labels <- paste0(
  "Top-",
  rank_cutoffs
)

model_zoneID <- c(
  1:7,
  9:50,
  52:55
)

dual_threshold <- 0.4
tie_tol <- 1e-4
novel_value <- 99L
chunk_rows <- 30L
reuse_cache <- TRUE
cache_version <- "v1"

# Relative tolerance for independent cell-area summations.
area_check_tolerance <- 1e-7


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


safe_name <- function(x) {
  x <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    x
  )
  
  x <- gsub(
    "^_+|_+$",
    "",
    x
  )
  
  ifelse(
    nchar(x) == 0,
    "unnamed",
    x
  )
}


save_plot <- function(
    plot_object,
    filename,
    width,
    height) {
  ggsave(
    filename = filename,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
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


ranked_suitability_file <- function(
    method,
    scenario) {
  file.path(
    ranking_root,
    method,
    scenario,
    "ranked_suitability.tif"
  )
}


cache_is_current <- function(
    cache_files,
    input_files) {
  if (!reuse_cache) {
    return(FALSE)
  }
  
  cache_files <- unlist(
    cache_files
  )
  
  input_files <- unlist(
    input_files
  )
  
  if (!all(file.exists(cache_files))) {
    return(FALSE)
  }
  
  if (!all(file.exists(input_files))) {
    return(FALSE)
  }
  
  min(
    file.info(cache_files)$mtime
  ) >= max(
    file.info(input_files)$mtime
  )
}


check_geometry <- function(
    raster,
    template,
    description) {
  if (!compareGeom(
    raster,
    template,
    stopOnError = FALSE
  )) {
    stop(
      "Raster geometry mismatch: ",
      description
    )
  }
  
  invisible(TRUE)
}


update_zone_totals <- function(
    area_total,
    pixel_total,
    zone,
    area,
    keep) {
  if (!any(keep)) {
    return(list(
      area = area_total,
      pixels = pixel_total
    ))
  }
  
  zone_index <- match(
    zone[keep],
    model_zoneID
  )
  
  pixel_total <- pixel_total +
    tabulate(
      zone_index,
      nbins = length(model_zoneID)
    )
  
  area_sum <- rowsum(
    area[keep],
    zone_index,
    reorder = FALSE
  )
  
  area_index <- as.integer(
    rownames(area_sum)
  )
  
  area_total[area_index] <-
    area_total[area_index] +
    as.numeric(
      area_sum[, 1]
    )
  
  list(
    area = area_total,
    pixels = pixel_total
  )
}


relative_area_error <- function(
    new_value,
    old_value) {
  abs(
    new_value - old_value
  ) / pmax(
    1,
    abs(old_value)
  )
}


# 2. Read lookup data and existing Top-1 results ================================

population <- fread(
  require_file(
    population_lookup_file
  )
)

required_population_columns <- c(
  "PopulationID",
  "Species",
  "source_zone"
)

missing_population_columns <- setdiff(
  required_population_columns,
  names(population)
)

if (length(missing_population_columns) > 0L) {
  stop(
    "Population lookup is missing columns: ",
    paste(
      missing_population_columns,
      collapse = ", "
    )
  )
}

population[
  ,
  `:=`(
    PopulationID = as.character(PopulationID),
    Species = as.character(Species),
    source_zone = as.integer(source_zone)
  )
]

if ("projected" %in% names(population)) {
  population[
    ,
    projected := as.logical(projected)
  ]
  
  population <- population[
    projected == TRUE
  ]
}

population <- population[
  source_zone %in% model_zoneID
]

if (!nrow(population)) {
  stop(
    "No modeled source populations remain in the lookup table."
  )
}

if (anyDuplicated(
  population[
    ,
    .(
      Species,
      source_zone
    )
  ]
)) {
  stop(
    "Duplicated species x source-zone populations were found."
  )
}

species_names <- sort(
  unique(population$Species)
)

if (length(species_names) > 30L) {
  stop(
    "The bit-mask implementation supports at most 30 species."
  )
}

species_index <- setNames(
  seq_along(species_names),
  species_names
)

species_bits <- bitwShiftL(
  1L,
  seq_along(species_names) - 1L
)

zone_mask_lookup <- integer(
  max(
    c(
      model_zoneID,
      novel_value
    )
  ) + 1L
)

for (row_index in seq_len(nrow(population))) {
  zone_index <- population$source_zone[row_index] + 1L
  bit_value <- species_bits[
    species_index[
      population$Species[row_index]
    ]
  ]
  
  zone_mask_lookup[zone_index] <- bitwOr(
    zone_mask_lookup[zone_index],
    bit_value
  )
}

existing_ecosystem <- fread(
  require_file(
    file.path(
      assigned_table_dir,
      "future_ecosystem_area_var.csv"
    )
  )
)[
  method %in% method_order
]

existing_species <- fread(
  require_file(
    file.path(
      assigned_table_dir,
      "future_species_niche_area_var.csv"
    )
  )
)[
  method %in% method_order
]

existing_population <- fread(
  require_file(
    file.path(
      assigned_table_dir,
      "future_population_niche_area_var.csv"
    )
  )
)[
  method %in% method_order
]

future_topk_retention <- fread(
  require_file(
    future_topk_retention_file
  )
)[
  method %in% method_order &
    scenario %in% future_order
]

required_future_retention_columns <- c(
  "method",
  "method_label",
  "scenario",
  "period",
  "ssp",
  "common_valid_pixels",
  "common_valid_area_km2",
  "stable_top1_pixels",
  "stable_top1_area_km2",
  "changed_existing_pixels",
  "changed_existing_area_km2",
  "novel_pixels",
  "novel_area_km2",
  "former_top3_among_changed_pixel",
  "former_top3_among_changed_area",
  "former_top5_among_changed_pixel",
  "former_top5_among_changed_area"
)

missing_future_retention_columns <- setdiff(
  required_future_retention_columns,
  names(future_topk_retention)
)

if (length(missing_future_retention_columns) > 0L) {
  stop(
    "Future Top-k retention table is missing columns: ",
    paste(
      missing_future_retention_columns,
      collapse = ", "
    )
  )
}

expected_future_retention_jobs <- CJ(
  method = method_order,
  scenario = future_order,
  unique = TRUE
)

observed_future_retention_jobs <- unique(
  future_topk_retention[
    ,
    .(
      method,
      scenario
    )
  ]
)

missing_future_retention_jobs <- fsetdiff(
  expected_future_retention_jobs,
  observed_future_retention_jobs
)

duplicated_future_retention_jobs <- future_topk_retention[
  duplicated(
    future_topk_retention,
    by = c(
      "method",
      "scenario"
    )
  ) |
    duplicated(
      future_topk_retention,
      by = c(
        "method",
        "scenario"
      ),
      fromLast = TRUE
    )
]

if (nrow(missing_future_retention_jobs) > 0L) {
  print(
    missing_future_retention_jobs
  )
  stop(
    "Future Top-k retention table is incomplete. Re-run script 8.2."
  )
}

if (nrow(duplicated_future_retention_jobs) > 0L) {
  print(
    duplicated_future_retention_jobs[
      order(
        method,
        scenario
      )
    ]
  )
  stop(
    "Future Top-k retention table has duplicated method x scenario rows."
  )
}

if (nrow(future_topk_retention) !=
    nrow(expected_future_retention_jobs)) {
  stop(
    "Future Top-k retention table has an unexpected number of rows."
  )
}

future_retention_fraction_columns <- c(
  "former_top3_among_changed_pixel",
  "former_top3_among_changed_area",
  "former_top5_among_changed_pixel",
  "former_top5_among_changed_area"
)

if (anyNA(
  future_topk_retention[
    ,
    ..future_retention_fraction_columns
  ]
) || any(
  future_topk_retention[
    ,
    unlist(
      lapply(
        .SD,
        function(x) {
          !is.finite(x) |
            x < 0 |
            x > 1
        }
      )
    ),
    .SDcols = future_retention_fraction_columns
  ]
)) {
  stop(
    "Future Top-k retention fractions must be finite values within [0, 1]."
  )
}

if (any(
  future_topk_retention$former_top3_among_changed_pixel >
  future_topk_retention$former_top5_among_changed_pixel +
  1e-12
) || any(
  future_topk_retention$former_top3_among_changed_area >
  future_topk_retention$former_top5_among_changed_area +
  1e-12
)) {
  stop(
    "Future retention is not nested: expected Top-3 <= Top-5."
  )
}

reference_map <- rast(
  require_file(
    reference_file
  )
)

cell_area <- cellSize(
  reference_map,
  unit = "km"
)

names(cell_area) <- "cell_area_km2"


# 3. Scan each method x scenario once ===========================================

analogue_results <- list()
zone_results <- list()
species_results <- list()

result_index <- 0L

for (method in method_order) {
  for (scenario in future_order) {
    result_index <- result_index + 1L
    
    cat(
      "\n[TOP-K DOWNSTREAM] ",
      method,
      " | ",
      scenario,
      "\n",
      sep = ""
    )
    
    assigned_file <- assigned_map_file(
      method,
      scenario
    )
    
    rank_zone_file <- ranked_zone_file(
      method,
      scenario
    )
    
    rank_suit_file <- ranked_suitability_file(
      method,
      scenario
    )
    
    input_files <- c(
      reference_file,
      population_lookup_file,
      assigned_file,
      rank_zone_file,
      rank_suit_file
    )
    
    lapply(
      input_files,
      require_file
    )
    
    cache_prefix <- paste0(
      cache_version,
      "_",
      method,
      "_",
      scenario
    )
    
    cache_analogue <- file.path(
      cache_dir,
      paste0(
        cache_prefix,
        "_analogue.csv"
      )
    )
    
    cache_zone <- file.path(
      cache_dir,
      paste0(
        cache_prefix,
        "_zone.csv"
      )
    )
    
    cache_species <- file.path(
      cache_dir,
      paste0(
        cache_prefix,
        "_species.csv"
      )
    )
    
    cache_files <- c(
      cache_analogue,
      cache_zone,
      cache_species
    )
    
    if (cache_is_current(
      cache_files,
      input_files
    )) {
      cat(
        "[REUSE CACHE]\n"
      )
      
      analogue_job <- fread(
        cache_analogue
      )
      
      zone_job <- fread(
        cache_zone
      )
      
      species_job <- fread(
        cache_species
      )
      
    } else {
      assigned_map <- rast(
        assigned_file
      )[[1]]
      
      rank_zone <- rast(
        rank_zone_file
      )[[1:5]]
      
      rank_suitability <- rast(
        rank_suit_file
      )[[1:5]]
      
      check_geometry(
        assigned_map,
        reference_map,
        assigned_file
      )
      
      check_geometry(
        rank_zone,
        reference_map,
        rank_zone_file
      )
      
      check_geometry(
        rank_suitability,
        reference_map,
        rank_suit_file
      )
      
      scarcity_area <- numeric(
        length(rank_cutoffs)
      )
      
      scarcity_pixels <- numeric(
        length(rank_cutoffs)
      )
      
      common_area <- 0
      common_pixels <- 0
      
      zone_area <- matrix(
        0,
        nrow = length(model_zoneID),
        ncol = length(rank_cutoffs)
      )
      
      zone_pixels <- matrix(
        0,
        nrow = length(model_zoneID),
        ncol = length(rank_cutoffs)
      )
      
      species_area <- matrix(
        0,
        nrow = length(species_names),
        ncol = length(rank_cutoffs)
      )
      
      species_pixels <- matrix(
        0,
        nrow = length(species_names),
        ncol = length(rank_cutoffs)
      )
      
      readStart(assigned_map)
      readStart(rank_zone)
      readStart(rank_suitability)
      readStart(cell_area)
      
      for (row_start in seq(
        1L,
        nrow(reference_map),
        by = chunk_rows
      )) {
        nrows_now <- min(
          chunk_rows,
          nrow(reference_map) -
            row_start + 1L
        )
        
        assigned <- as.integer(
          readValues(
            assigned_map,
            row = row_start,
            nrows = nrows_now,
            mat = FALSE
          )
        )
        
        ranked_zones <- readValues(
          rank_zone,
          row = row_start,
          nrows = nrows_now,
          mat = TRUE
        )
        
        ranked_suitability <- readValues(
          rank_suitability,
          row = row_start,
          nrows = nrows_now,
          mat = TRUE
        )
        
        area <- as.numeric(
          readValues(
            cell_area,
            row = row_start,
            nrows = nrows_now,
            mat = FALSE
          )
        )
        
        if (is.null(dim(ranked_zones))) {
          ranked_zones <- matrix(
            ranked_zones,
            ncol = 5L
          )
        }
        
        if (is.null(dim(ranked_suitability))) {
          ranked_suitability <- matrix(
            ranked_suitability,
            ncol = 5L
          )
        }
        
        storage.mode(ranked_zones) <- "integer"
        storage.mode(ranked_suitability) <- "double"
        
        valid <- (
          assigned %in% c(
            model_zoneID,
            novel_value
          ) &
            is.finite(area)
        )
        
        if (!any(valid)) {
          next
        }
        
        common_pixels <- common_pixels +
          sum(valid)
        
        common_area <- common_area +
          sum(
            area[valid],
            na.rm = TRUE
          )
        
        above_threshold <- (
          is.finite(ranked_suitability) &
            ranked_suitability >= dual_threshold
        )
        
        n_above_threshold <- rowSums(
          above_threshold,
          na.rm = TRUE
        )
        
        for (rank_index in seq_along(rank_cutoffs)) {
          k <- rank_cutoffs[rank_index]
          
          scarce <- (
            valid &
              n_above_threshold < k
          )
          
          scarcity_pixels[rank_index] <-
            scarcity_pixels[rank_index] +
            sum(scarce)
          
          scarcity_area[rank_index] <-
            scarcity_area[rank_index] +
            sum(
              area[scarce],
              na.rm = TRUE
            )
          
          assigned_candidate <- assigned
          assigned_candidate[
            !(assigned_candidate %in% model_zoneID)
          ] <- NA_integer_
          
          if (k == 1L) {
            candidates <- matrix(
              assigned_candidate,
              ncol = 1L
            )
          } else {
            ranked_candidates <- ranked_zones[
              ,
              seq_len(k),
              drop = FALSE
            ]
            
            ranked_candidates[
              !above_threshold[
                ,
                seq_len(k),
                drop = FALSE
              ]
            ] <- NA_integer_
            
            candidates <- cbind(
              assigned_candidate,
              ranked_candidates
            )
          }
          
          # Accumulate each candidate zone once per cell. Ranked zones are
          # unique; this also removes duplication with the assigned Top-1 zone.
          for (candidate_index in seq_len(ncol(candidates))) {
            candidate <- candidates[
              ,
              candidate_index
            ]
            
            duplicated_here <- rep(
              FALSE,
              length(candidate)
            )
            
            if (candidate_index > 1L) {
              duplicated_here <- rowSums(
                candidates[
                  ,
                  seq_len(candidate_index - 1L),
                  drop = FALSE
                ] == candidate,
                na.rm = TRUE
              ) > 0
            }
            
            keep_candidate <- (
              valid &
                candidate %in% model_zoneID &
                !duplicated_here
            )
            
            updated_totals <- update_zone_totals(
              area_total = zone_area[
                ,
                rank_index
              ],
              pixel_total = zone_pixels[
                ,
                rank_index
              ],
              zone = candidate,
              area = area,
              keep = keep_candidate
            )
            
            zone_area[
              ,
              rank_index
            ] <- updated_totals$area
            
            zone_pixels[
              ,
              rank_index
            ] <- updated_totals$pixels
          }
          
          # Convert candidate ecotypes to a bit mask so each cell is counted
          # only once for a species represented by multiple source ecotypes.
          cell_species_mask <- integer(
            nrow(candidates)
          )
          
          for (candidate_index in seq_len(ncol(candidates))) {
            candidate <- candidates[
              ,
              candidate_index
            ]
            
            candidate_mask <- integer(
              length(candidate)
            )
            
            keep_lookup <- (
              !is.na(candidate) &
                candidate >= 0L &
                candidate < length(zone_mask_lookup)
            )
            
            candidate_mask[keep_lookup] <-
              zone_mask_lookup[
                candidate[keep_lookup] + 1L
              ]
            
            cell_species_mask <- bitwOr(
              cell_species_mask,
              candidate_mask
            )
          }
          
          for (species_number in seq_along(species_names)) {
            species_present <- (
              valid &
                bitwAnd(
                  cell_species_mask,
                  species_bits[species_number]
                ) != 0L
            )
            
            species_pixels[
              species_number,
              rank_index
            ] <- species_pixels[
              species_number,
              rank_index
            ] + sum(species_present)
            
            species_area[
              species_number,
              rank_index
            ] <- species_area[
              species_number,
              rank_index
            ] + sum(
              area[species_present],
              na.rm = TRUE
            )
          }
        }
      }
      
      readStop(assigned_map)
      readStop(rank_zone)
      readStop(rank_suitability)
      readStop(cell_area)
      
      fields <- scenario_fields(
        scenario
      )
      
      analogue_job <- data.table(
        method = method,
        method_label = method_labels[[method]],
        scenario = scenario,
        period = fields$period,
        ssp = fields$ssp,
        rank_cutoff = rank_cutoffs,
        rank_label = rank_labels,
        criterion = paste0(
          "fewer_than_",
          rank_cutoffs,
          "_ecotypes_at_or_above_",
          dual_threshold
        ),
        common_valid_pixels = common_pixels,
        common_valid_area_km2 = common_area,
        pixel_count = scarcity_pixels,
        pixel_share = scarcity_pixels /
          common_pixels,
        area_km2 = scarcity_area,
        area_share = scarcity_area /
          common_area
      )
      
      zone_job <- CJ(
        source_zone = model_zoneID,
        rank_index = seq_along(rank_cutoffs),
        sorted = FALSE
      )
      
      zone_job[
        ,
        `:=`(
          method = method,
          method_label = method_labels[[method]],
          scenario = scenario,
          period = fields$period,
          ssp = fields$ssp,
          common_valid_pixels = common_pixels,
          common_valid_area_km2 = common_area,
          rank_cutoff = rank_cutoffs[rank_index],
          rank_label = rank_labels[rank_index],
          pixel_count = as.numeric(
            zone_pixels[
              cbind(
                match(source_zone, model_zoneID),
                rank_index
              )
            ]
          ),
          area_km2 = as.numeric(
            zone_area[
              cbind(
                match(source_zone, model_zoneID),
                rank_index
              )
            ]
          )
        )
      ]
      
      zone_job[
        ,
        rank_index := NULL
      ]
      
      zone_job[
        ,
        `:=`(
          pixel_share = pixel_count /
            common_valid_pixels,
          area_share = area_km2 /
            common_valid_area_km2
        )
      ]
      
      species_job <- CJ(
        Species = species_names,
        rank_index = seq_along(rank_cutoffs),
        sorted = FALSE
      )
      
      species_job[
        ,
        `:=`(
          method = method,
          method_label = method_labels[[method]],
          scenario = scenario,
          period = fields$period,
          ssp = fields$ssp,
          common_valid_pixels = common_pixels,
          common_valid_area_km2 = common_area,
          rank_cutoff = rank_cutoffs[rank_index],
          rank_label = rank_labels[rank_index],
          pixel_count = as.numeric(
            species_pixels[
              cbind(
                match(Species, species_names),
                rank_index
              )
            ]
          ),
          area_km2 = as.numeric(
            species_area[
              cbind(
                match(Species, species_names),
                rank_index
              )
            ]
          )
        )
      ]
      
      species_job[
        ,
        rank_index := NULL
      ]
      
      species_job[
        ,
        `:=`(
          pixel_share = pixel_count /
            common_valid_pixels,
          area_share = area_km2 /
            common_valid_area_km2
        )
      ]
      
      fwrite(
        analogue_job,
        cache_analogue
      )
      
      fwrite(
        zone_job,
        cache_zone
      )
      
      fwrite(
        species_job,
        cache_species
      )
      
      rm(
        assigned_map,
        rank_zone,
        rank_suitability
      )
      
      gc()
    }
    
    analogue_results[[result_index]] <-
      analogue_job
    
    zone_results[[result_index]] <-
      zone_job
    
    species_results[[result_index]] <-
      species_job
  }
}

analogue_topk <- rbindlist(
  analogue_results,
  fill = TRUE
)

zone_topk <- rbindlist(
  zone_results,
  fill = TRUE
)

species_topk <- rbindlist(
  species_results,
  fill = TRUE
)

population_topk <- merge(
  population,
  zone_topk,
  by = "source_zone",
  allow.cartesian = TRUE,
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(
  population_topk$area_km2
)) {
  stop(
    "Top-k area is missing for at least one source population."
  )
}

# The scarcity metric asks whether the number of qualifying analogues is below
# an increasingly demanding minimum. Its direction is therefore the opposite
# of rank-based non-retention: scarcity must satisfy count-1 <= count-3 <=
# count-5 for every method and scenario.
scarcity_area_wide <- dcast(
  analogue_topk,
  method + scenario ~ rank_cutoff,
  value.var = "area_km2"
)

scarcity_pixel_wide <- dcast(
  analogue_topk,
  method + scenario ~ rank_cutoff,
  value.var = "pixel_count"
)

required_rank_columns <- as.character(
  rank_cutoffs
)

if (!all(
  required_rank_columns %in%
  names(scarcity_area_wide)
) || !all(
  required_rank_columns %in%
  names(scarcity_pixel_wide)
)) {
  stop(
    "Failed to construct all minimum-1/minimum-3/minimum-5 scarcity columns."
  )
}

if (any(
  scarcity_area_wide[["1"]] >
  scarcity_area_wide[["3"]] +
  1e-6 |
  scarcity_area_wide[["3"]] >
  scarcity_area_wide[["5"]] +
  1e-6
) || any(
  scarcity_pixel_wide[["1"]] >
  scarcity_pixel_wide[["3"]] |
  scarcity_pixel_wide[["3"]] >
  scarcity_pixel_wide[["5"]]
)) {
  stop(
    "Analogue scarcity is not monotonic: expected minimum-1 <= minimum-3 <= minimum-5."
  )
}


# 4. Future rank non-retention of the normal-period Top-1 ecotype ==============

# This is distinct from novel niche space. Top-1 non-retention includes every
# changed or novel cell. At Top-3 and Top-5, changed cells are removed from the
# non-retained set when their normal-period Top-1 ecotype remains within the
# corresponding future rank set. Novel cells remain non-retained at every k.
# For non-novel cells, this is a relative-rank diagnostic; the former ecotype
# itself is not additionally required to reach 0.4.
rank_nonretention <- rbindlist(
  lapply(
    rank_cutoffs,
    function(k) {
      result <- copy(
        future_topk_retention
      )
      
      if (k == 1L) {
        retained_changed_pixel_fraction <- 0
        retained_changed_area_fraction <- 0
      } else if (k == 3L) {
        retained_changed_pixel_fraction <-
          result$former_top3_among_changed_pixel
        retained_changed_area_fraction <-
          result$former_top3_among_changed_area
      } else if (k == 5L) {
        retained_changed_pixel_fraction <-
          result$former_top5_among_changed_pixel
        retained_changed_area_fraction <-
          result$former_top5_among_changed_area
      } else {
        stop(
          "Unsupported rank cutoff: ",
          k
        )
      }
      
      result[
        ,
        `:=`(
          rank_cutoff = k,
          rank_label = paste0(
            "Top-",
            k
          ),
          changed_nonretained_pixels = round(
            changed_existing_pixels *
              (
                1 -
                  retained_changed_pixel_fraction
              )
          ),
          changed_nonretained_area_km2 =
            changed_existing_area_km2 *
            (
              1 -
                retained_changed_area_fraction
            )
        )
      ]
      
      result[
        ,
        `:=`(
          pixel_count =
            changed_nonretained_pixels +
            novel_pixels,
          area_km2 =
            changed_nonretained_area_km2 +
            novel_area_km2,
          criterion = paste0(
            "normal_period_top1_not_retained_within_future_top",
            k,
            "_including_novel"
          ),
          baseline =
            "workflow-specific normal-period Top-1 assigned map"
        )
      ]
      
      result[
        ,
        `:=`(
          pixel_share =
            pixel_count /
            common_valid_pixels,
          area_share =
            area_km2 /
            common_valid_area_km2,
          novel_pixel_share =
            novel_pixels /
            common_valid_pixels,
          novel_area_share =
            novel_area_km2 /
            common_valid_area_km2
        )
      ]
      
      result[
        ,
        .(
          method,
          method_label,
          scenario,
          period,
          ssp,
          baseline,
          rank_cutoff,
          rank_label,
          criterion,
          common_valid_pixels,
          common_valid_area_km2,
          stable_top1_pixels,
          stable_top1_area_km2,
          changed_nonretained_pixels,
          changed_nonretained_area_km2,
          novel_pixels,
          novel_area_km2,
          pixel_count,
          pixel_share,
          area_km2,
          area_share,
          novel_pixel_share,
          novel_area_share
        )
      ]
    }
  ),
  fill = TRUE
)

nonretention_area_wide <- dcast(
  rank_nonretention,
  method + scenario ~ rank_cutoff,
  value.var = "area_km2"
)

nonretention_pixel_wide <- dcast(
  rank_nonretention,
  method + scenario ~ rank_cutoff,
  value.var = "pixel_count"
)

if (!all(
  required_rank_columns %in%
  names(nonretention_area_wide)
) || !all(
  required_rank_columns %in%
  names(nonretention_pixel_wide)
)) {
  stop(
    "Failed to construct all Top-1/Top-3/Top-5 non-retention columns."
  )
}

# Relaxing the retention criterion cannot increase non-retained area or pixels.
if (any(
  nonretention_area_wide[["1"]] +
  1e-6 <
  nonretention_area_wide[["3"]] |
  nonretention_area_wide[["3"]] +
  1e-6 <
  nonretention_area_wide[["5"]]
) || any(
  nonretention_pixel_wide[["1"]] <
  nonretention_pixel_wide[["3"]] |
  nonretention_pixel_wide[["3"]] <
  nonretention_pixel_wide[["5"]]
)) {
  stop(
    "Top-k non-retention is not monotonic: expected Top-1 >= Top-3 >= Top-5."
  )
}

top1_nonretention_check <- rank_nonretention[
  rank_cutoff == 1L
]

stopifnot(
  max(
    abs(
      top1_nonretention_check$area_km2 -
        (
          top1_nonretention_check$common_valid_area_km2 -
            top1_nonretention_check$stable_top1_area_km2
        )
    )
  ) < 1e-6,
  max(
    abs(
      top1_nonretention_check$pixel_count -
        (
          top1_nonretention_check$common_valid_pixels -
            top1_nonretention_check$stable_top1_pixels
        )
    )
  ) == 0
)

# Zone 99 is independent of k. Confirm that the novel component imported from
# script 8.2 exactly matches the direct count of cells with no ecotype reaching
# dual suitability 0.4 in the ranked-suitability rasters.
novel_definition_check <- merge(
  rank_nonretention[
    rank_cutoff == 1L,
    .(
      method,
      scenario,
      retention_table_novel_pixels = novel_pixels,
      retention_table_novel_area_km2 = novel_area_km2
    )
  ],
  analogue_topk[
    rank_cutoff == 1L,
    .(
      method,
      scenario,
      threshold_scan_novel_pixels = pixel_count,
      threshold_scan_novel_area_km2 = area_km2
    )
  ],
  by = c(
    "method",
    "scenario"
  ),
  all = TRUE,
  sort = TRUE
)

novel_definition_check[
  ,
  `:=`(
    pixel_difference =
      retention_table_novel_pixels -
      threshold_scan_novel_pixels,
    area_km2_difference =
      retention_table_novel_area_km2 -
      threshold_scan_novel_area_km2
  )
]

stopifnot(
  nrow(novel_definition_check) ==
    nrow(expected_future_retention_jobs),
  !anyNA(novel_definition_check),
  max(
    abs(
      novel_definition_check$pixel_difference
    )
  ) == 0,
  max(
    abs(
      novel_definition_check$area_km2_difference
    )
  ) < 1e-6
)


# 5. Top-1 consistency checks ===================================================

novel_check <- merge(
  analogue_topk[
    rank_cutoff == 1L,
    .(
      method,
      scenario,
      new_area_km2 = area_km2
    )
  ],
  existing_ecosystem[
    zoneID == novel_value,
    .(
      method,
      scenario,
      old_area_km2 = area_km2
    )
  ],
  by = c(
    "method",
    "scenario"
  ),
  all.x = TRUE
)

# terra::zonal omits a category that has zero cells. In that case the absent
# Zone-99 row in the existing table represents an area of exactly zero.
novel_check[
  is.na(old_area_km2),
  old_area_km2 := 0
]

species_check <- merge(
  species_topk[
    rank_cutoff == 1L,
    .(
      method,
      scenario,
      Species,
      new_area_km2 = area_km2
    )
  ],
  existing_species[
    ,
    .(
      method,
      scenario,
      Species,
      old_area_km2 = future_area_km2
    )
  ],
  by = c(
    "method",
    "scenario",
    "Species"
  ),
  all = TRUE
)

population_check <- merge(
  population_topk[
    rank_cutoff == 1L,
    .(
      method,
      scenario,
      Species,
      PopulationID,
      source_zone,
      new_area_km2 = area_km2
    )
  ],
  existing_population[
    ,
    .(
      method,
      scenario,
      Species,
      PopulationID = as.character(PopulationID),
      source_zone = as.integer(source_zone),
      old_area_km2 = future_area_km2
    )
  ],
  by = c(
    "method",
    "scenario",
    "Species",
    "PopulationID",
    "source_zone"
  ),
  all = TRUE
)

for (check_name in c(
  "novel_check",
  "species_check",
  "population_check"
)) {
  check_table <- get(
    check_name
  )
  
  if (anyNA(check_table)) {
    stop(
      "Top-1 consistency table contains missing values: ",
      check_name
    )
  }
  
  check_table[
    ,
    relative_error := relative_area_error(
      new_area_km2,
      old_area_km2
    )
  ]
  
  if (max(
    check_table$relative_error,
    na.rm = TRUE
  ) > area_check_tolerance) {
    print(
      check_table[
        order(-relative_error)
      ][
        seq_len(
          min(
            nrow(check_table),
            20L
          )
        )
      ]
    )
    
    stop(
      "Top-1 area does not reproduce the existing result: ",
      check_name
    )
  }
  
  fwrite(
    check_table,
    file.path(
      table_dir,
      paste0(
        "Figure_topk_",
        check_name,
        ".csv"
      )
    )
  )
}


# 6. Save final source tables ===================================================

fwrite(
  rank_nonretention,
  file.path(
    table_dir,
    "Figure_var_6b_normal_top1_rank_nonretention.csv"
  )
)

fwrite(
  novel_definition_check,
  file.path(
    table_dir,
    "Figure_var_6_vs_6b_novel_definition_check.csv"
  )
)

fwrite(
  analogue_topk,
  file.path(
    table_dir,
    "Figure_S_analogue_scarcity_by_minimum_count.csv"
  )
)

fwrite(
  zone_topk,
  file.path(
    table_dir,
    "Figure_var_topk_source_ecotype_area.csv"
  )
)

fwrite(
  species_topk,
  file.path(
    table_dir,
    "Figure_var_8_species_topk_area.csv"
  )
)

fwrite(
  population_topk,
  file.path(
    table_dir,
    "Figure_var_10a_population_topk_area.csv"
  )
)


# 7. Shared plotting fields =====================================================

for (table_name in c(
  "rank_nonretention",
  "analogue_topk",
  "species_topk",
  "population_topk"
)) {
  table_object <- get(
    table_name
  )
  
  table_object[
    ,
    `:=`(
      period = factor(
        period,
        levels = period_levels
      ),
      ssp = factor(
        ssp,
        levels = ssp_levels
      ),
      rank_label = factor(
        rank_label,
        levels = rank_labels
      ),
      method_label = factor(
        method_label,
        levels = method_labels[method_order]
      )
    )
  ]
  
  assign(
    table_name,
    table_object
  )
}


# 8. Figure 6b: rank non-retention of the normal-period Top-1 ecotype ===========

rank_nonretention[
  ,
  area_million_km2 := area_km2 /
    1e6
]

figure_6b_topk <- ggplot(
  rank_nonretention,
  aes(
    x = period,
    y = area_million_km2,
    group = method_label,
    shape = method_label,
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.78,
    color = "black"
  ) +
  geom_point(
    size = 2.0,
    fill = "white",
    color = "black",
    stroke = 0.5
  ) +
  facet_grid(
    rank_label ~ ssp
  ) +
  scale_shape_manual(
    values = c(
      21,
      24
    )
  ) +
  scale_linetype_manual(
    values = c(
      "solid",
      "22"
    )
  ) +
  labs(
    x = "Future period",
    y = expression("Area without rank retention (million km"^2*")"),
    shape = NULL,
    linetype = NULL,
    title = "Future rank non-retention of the normal-period Top-1 ecotype",
    subtitle = paste(
      "A cell is counted when its normal-period Top-1 ecotype is absent from the first k future ranks.",
      "Novel cells are included at every k; relaxing k can only reduce area without rank retention.",
      sep = "\n"
    )
  ) +
  theme_bw(
    base_size = 10.2
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    panel.spacing = grid::unit(
      0.8,
      "lines"
    )
  )

save_plot(
  figure_6b_topk,
  file.path(
    figure_dir,
    "Figure_var_6b_normal_top1_rank_nonretention_topk.png"
  ),
  9.2,
  8.6
)

common_y_max_6b <- max(
  rank_nonretention$area_million_km2,
  na.rm = TRUE
) * 1.04

for (k in rank_cutoffs) {
  figure_6b_k <- ggplot(
    rank_nonretention[
      rank_cutoff == k
    ],
    aes(
      x = period,
      y = area_million_km2,
      group = method_label,
      shape = method_label,
      linetype = method_label
    )
  ) +
    geom_line(
      linewidth = 0.82,
      color = "black"
    ) +
    geom_point(
      size = 2.0,
      fill = "white",
      color = "black",
      stroke = 0.5
    ) +
    facet_wrap(
      ~ ssp,
      nrow = 1
    ) +
    scale_shape_manual(
      values = c(
        21,
        24
      )
    ) +
    scale_linetype_manual(
      values = c(
        "solid",
        "22"
      )
    ) +
    coord_cartesian(
      ylim = c(
        0,
        common_y_max_6b
      )
    ) +
    labs(
      x = "Future period",
      y = expression("Area without rank retention (million km"^2*")"),
      shape = NULL,
      linetype = NULL,
      title = paste0(
        "Normal-period Top-1 ecotype absent from future Top-",
        k
      ),
      subtitle =
        "Relative-rank diagnostic; each workflow's normal-period Top-1 map is the baseline and novel cells are included."
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  save_plot(
    figure_6b_k,
    file.path(
      figure_dir,
      paste0(
        "Figure_var_6b_normal_top1_rank_nonretention_top",
        k,
        ".png"
      )
    ),
    8.5,
    4.8
  )
}


# 8b. Supplement: absolute scarcity of suitable current analogues ==============

analogue_topk[
  ,
  `:=`(
    area_million_km2 = area_km2 /
      1e6,
    minimum_count_label = factor(
      fifelse(
        rank_cutoff == 1L,
        "No suitable analogue",
        paste0(
          "Fewer than ",
          rank_cutoff,
          " suitable analogues"
        )
      ),
      levels = c(
        "No suitable analogue",
        "Fewer than 3 suitable analogues",
        "Fewer than 5 suitable analogues"
      )
    )
  )
]

figure_s_analogue_scarcity <- ggplot(
  analogue_topk,
  aes(
    x = period,
    y = area_million_km2,
    group = method_label,
    shape = method_label,
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.78,
    color = "black"
  ) +
  geom_point(
    size = 2.0,
    fill = "white",
    color = "black",
    stroke = 0.5
  ) +
  facet_grid(
    minimum_count_label ~ ssp
  ) +
  scale_shape_manual(
    values = c(
      21,
      24
    )
  ) +
  scale_linetype_manual(
    values = c(
      "solid",
      "22"
    )
  ) +
  labs(
    x = "Future period",
    y = expression("Area below the minimum analogue count (million km"^2*")"),
    shape = NULL,
    linetype = NULL,
    title = "Scarcity of suitable current-ecotype analogues",
    subtitle = paste(
      "A suitable analogue has dual suitability >= 0.4.",
      "Rows impose increasingly demanding minimum candidate counts, so scarcity increases with the count.",
      sep = "\n"
    )
  ) +
  theme_bw(
    base_size = 10.2
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    panel.spacing = grid::unit(
      0.8,
      "lines"
    )
  )

save_plot(
  figure_s_analogue_scarcity,
  file.path(
    figure_dir,
    "Figure_S_analogue_scarcity_by_minimum_count.png"
  ),
  9.2,
  8.6
)


# 9. Figure 8: species area within Top-k ecotype analogues ======================

species_topk[
  ,
  area_million_km2 := area_km2 /
    1e6
]

species_topk[
  ,
  Species := factor(
    Species,
    levels = sort(
      unique(
        as.character(Species)
      )
    )
  )
]

figure_8_topk <- ggplot(
  species_topk,
  aes(
    x = period,
    y = area_million_km2,
    group = method_label,
    shape = method_label,
    linetype = method_label
  )
) +
  geom_line(
    linewidth = 0.70,
    color = "black"
  ) +
  geom_point(
    size = 1.65,
    fill = "white",
    color = "black",
    stroke = 0.45
  ) +
  facet_grid(
    Species ~ rank_label + ssp,
    scales = "free_y"
  ) +
  scale_shape_manual(
    values = c(
      21,
      24
    )
  ) +
  scale_linetype_manual(
    values = c(
      "solid",
      "22"
    )
  ) +
  labs(
    x = "Future period",
    y = expression("Species area within Top-k analogues (million km"^2*")"),
    shape = NULL,
    linetype = NULL,
    title = "Projected species-niche area under Top-k ecotype-analogue criteria",
    subtitle = paste(
      "A cell is counted once when at least one source ecotype of the species",
      "is retained by the Top-k criterion at dual suitability >= 0.4."
    )
  ) +
  theme_bw(
    base_size = 8.6
  ) +
  theme(
    legend.position = "top",
    strip.text.y = element_text(
      angle = 0,
      face = "italic",
      size = 7.6
    ),
    strip.text.x = element_text(
      face = "bold",
      size = 7.6
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey93",
      linewidth = 0.25
    ),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1,
      size = 6.7
    ),
    panel.spacing = grid::unit(
      0.45,
      "lines"
    )
  )

save_plot(
  figure_8_topk,
  file.path(
    figure_dir,
    "Figure_var_8_assigned_species_area_topk.png"
  ),
  16.5,
  max(
    8.5,
    1.02 * length(species_names)
  )
)

# Separate Top-1, Top-3 and Top-5 versions retain the original Figure 8 layout.
for (k in rank_cutoffs) {
  figure_8_k <- ggplot(
    species_topk[
      rank_cutoff == k
    ],
    aes(
      x = period,
      y = area_million_km2,
      group = method_label,
      shape = method_label,
      linetype = method_label
    )
  ) +
    geom_line(
      linewidth = 0.75,
      color = "black"
    ) +
    geom_point(
      size = 1.8,
      fill = "white",
      color = "black",
      stroke = 0.48
    ) +
    facet_grid(
      Species ~ ssp,
      scales = "free_y"
    ) +
    scale_shape_manual(
      values = c(
        21,
        24
      )
    ) +
    scale_linetype_manual(
      values = c(
        "solid",
        "22"
      )
    ) +
    labs(
      x = "Future period",
      y = expression("Species area (million km"^2*")"),
      shape = NULL,
      linetype = NULL,
      title = paste0(
        "Projected species-niche area under the Top-",
        k,
        " ecotype-analogue criterion"
      ),
      subtitle = "Only ecotype analogues with dual suitability >= 0.4 are retained."
    ) +
    theme_bw(
      base_size = 9.2
    ) +
    theme(
      legend.position = "top",
      strip.text.y = element_text(
        angle = 0,
        face = "italic"
      ),
      strip.text.x = element_text(
        face = "bold"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(
        color = "grey94",
        linewidth = 0.25
      ),
      panel.spacing = grid::unit(
        0.80,
        "lines"
      )
    )
  
  save_plot(
    figure_8_k,
    file.path(
      figure_dir,
      paste0(
        "Figure_var_8_assigned_species_area_top",
        k,
        ".png"
      )
    ),
    10.4,
    max(
      8,
      1.15 * length(species_names)
    )
  )
}


# 10. Figure 10a: population source ecotypes within Top-k =======================

population_topk[
  ,
  population_label := if (
    "zone_name" %in% names(population_topk)
  ) {
    fifelse(
      !is.na(zone_name) &
        nzchar(zone_name),
      paste0(
        "P",
        PopulationID,
        " | Z",
        source_zone,
        " ",
        zone_name
      ),
      paste0(
        "P",
        PopulationID,
        " | Z",
        source_zone
      )
    )
  } else {
    paste0(
      "P",
      PopulationID,
      " | Z",
      source_zone
    )
  }
]

population_topk[
  ,
  area_million_km2 := area_km2 /
    1e6
]

population_species_summary <- population_topk[
  ,
  .(
    n_populations = uniqueN(
      PopulationID
    ),
    median_area_million_km2 = median(
      area_million_km2,
      na.rm = TRUE
    ),
    q25_area_million_km2 = quantile(
      area_million_km2,
      0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    q75_area_million_km2 = quantile(
      area_million_km2,
      0.75,
      na.rm = TRUE,
      names = FALSE
    )
  ),
  by = .(
    Species,
    method_label,
    scenario,
    period,
    ssp,
    rank_cutoff,
    rank_label
  )
]

population_species_summary[
  ,
  Species := factor(
    Species,
    levels = rev(
      sort(
        unique(Species)
      )
    )
  )
]

figure_10a_topk <- ggplot(
  population_species_summary,
  aes(
    x = median_area_million_km2,
    y = Species,
    shape = method_label
  )
) +
  geom_errorbarh(
    aes(
      xmin = q25_area_million_km2,
      xmax = q75_area_million_km2
    ),
    height = 0.16,
    position = position_dodge(
      width = 0.50
    ),
    linewidth = 0.48,
    color = "grey45"
  ) +
  geom_point(
    position = position_dodge(
      width = 0.50
    ),
    size = 2.0,
    fill = "white",
    color = "black",
    stroke = 0.50
  ) +
  facet_grid(
    rank_label ~ ssp + period,
    scales = "free_x"
  ) +
  scale_shape_manual(
    values = c(
      21,
      24
    )
  ) +
  labs(
    x = expression("Median population area within Top-k analogues (million km"^2*")"),
    y = "Species",
    shape = NULL,
    title = "Projected population-niche area under Top-k ecotype-analogue criteria",
    subtitle = paste(
      "Rows are Top-1, Top-3 and Top-5. Points are species medians across source",
      "populations; horizontal bars show interquartile ranges."
    )
  ) +
  theme_bw(
    base_size = 8.7
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(
      color = "grey92",
      linewidth = 0.28
    ),
    strip.text = element_text(
      face = "bold",
      size = 7.5
    ),
    axis.text.y = element_text(
      size = 7.2
    ),
    axis.text.x = element_text(
      size = 6.8
    ),
    panel.spacing = grid::unit(
      0.55,
      "lines"
    )
  )

save_plot(
  figure_10a_topk,
  file.path(
    figure_dir,
    "Figure_var_10a_population_species_summary_topk.png"
  ),
  16.0,
  8.5
)

fwrite(
  population_species_summary,
  file.path(
    table_dir,
    "Figure_var_10a_population_species_summary_topk.csv"
  )
)

# One detailed three-row figure per species. Rank rows share the x scale within
# each scenario column, so expansion from Top-1 to Top-3 to Top-5 is comparable.
population_figure_index <- list()

for (species_name in sort(
  unique(
    as.character(population_topk$Species)
  )
)) {
  species_dt <- copy(
    population_topk[
      Species == species_name
    ]
  )
  
  population_order <- species_dt[
    ,
    .(
      source_zone = first(source_zone)
    ),
    by = .(
      PopulationID,
      population_label
    )
  ][
    order(
      source_zone,
      PopulationID
    )
  ]$population_label
  
  species_dt[
    ,
    population_label := factor(
      population_label,
      levels = rev(population_order)
    )
  ]
  
  population_figure <- ggplot(
    species_dt,
    aes(
      x = area_million_km2,
      y = population_label,
      shape = method_label
    )
  ) +
    geom_point(
      position = position_dodge(
        width = 0.46
      ),
      size = 1.85,
      fill = "white",
      color = "black",
      stroke = 0.48
    ) +
    facet_grid(
      rank_label ~ ssp + period,
      scales = "free_x"
    ) +
    scale_shape_manual(
      values = c(
        21,
        24
      )
    ) +
    labs(
      x = expression("Population area within Top-k analogues (million km"^2*")"),
      y = "Population and source ecotype",
      shape = NULL,
      title = paste0(
        species_name,
        ": population-niche area under Top-k ecotype-analogue criteria"
      ),
      subtitle = paste(
        "Rows are Top-1, Top-3 and Top-5; only source ecotypes with dual",
        "suitability >= 0.4 are retained."
      )
    ) +
    theme_bw(
      base_size = 8.8
    ) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(
        color = "grey92",
        linewidth = 0.28
      ),
      strip.text = element_text(
        face = "bold",
        size = 7.4
      ),
      axis.text.y = element_text(
        size = 6.6
      ),
      axis.text.x = element_text(
        size = 6.7
      ),
      panel.spacing = grid::unit(
        0.52,
        "lines"
      )
    )
  
  population_figure_file <- file.path(
    population_detail_dir,
    paste0(
      "Figure_var_10a_population_topk_area_",
      safe_name(species_name),
      ".png"
    )
  )
  
  save_plot(
    population_figure,
    population_figure_file,
    16.0,
    max(
      7.5,
      4.9 +
        0.27 * length(population_order)
    )
  )
  
  population_figure_index[[
    length(population_figure_index) + 1L
  ]] <- data.table(
    Species = species_name,
    n_populations = length(population_order),
    figure_file = population_figure_file
  )
}

population_figure_index <- rbindlist(
  population_figure_index,
  fill = TRUE
)

fwrite(
  population_figure_index,
  file.path(
    table_dir,
    "Figure_var_10a_population_topk_figure_index.csv"
  )
)


# 11. Completion report =========================================================

obsolete_figure6_outputs <- c(
  file.path(
    figure_dir,
    "Figure_var_6_topk_analogue_availability.png"
  ),
  file.path(
    figure_dir,
    paste0(
      "Figure_var_6_analogue_availability_top",
      rank_cutoffs,
      ".png"
    )
  ),
  file.path(
    table_dir,
    "Figure_var_6_topk_analogue_availability.csv"
  )
)

obsolete_figure6_outputs <- obsolete_figure6_outputs[
  file.exists(
    obsolete_figure6_outputs
  )
]

if (length(obsolete_figure6_outputs) > 0L) {
  warning(
    paste(
      "Superseded Figure 6 files remain on disk and should not be used:",
      paste(
        obsolete_figure6_outputs,
        collapse = "\n"
      ),
      sep = "\n"
    ),
    call. = FALSE
  )
}

cat(
  "\nCOMPLETE\n",
  "Threshold: ",
  dual_threshold,
  "\n",
  "Ranks: ",
  paste(
    rank_labels,
    collapse = ", "
  ),
  "\n",
  "Top-k non-retention monotonicity: Top-1 >= Top-3 >= Top-5 PASS\n",
  "Analogue-scarcity monotonicity: minimum-1 <= minimum-3 <= minimum-5 PASS\n",
  "Top-1 consistency checks: PASS\n",
  "Figures: ",
  figure_dir,
  "\n",
  "Tables: ",
  table_dir,
  "\n",
  sep = ""
)
