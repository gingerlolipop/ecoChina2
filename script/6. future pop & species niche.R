# Future population and species niches from final Multi-Forest dual suitability
# ==============================================================================
# A population proxy is one species x source-ecotype combination. Its projected
# niche is the dual-suitability surface of that source ecotype. A species niche
# is the pixel-wise maximum across its projection-eligible populations.
#
# Run after scripts 1.2, 3.2 and 4. Legacy mf_var paths are retained so valid
# raster products from the original workflow are reused rather than rebuilt.

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

modeled_zones <- c(1:7, 9:50, 52:55)
scenarios <- c(
  "normal",
  "2011-2040SSP245", "2041-2070SSP245", "2071-2100SSP245",
  "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585"
)
future_scenarios <- setdiff(scenarios, "normal")

dual_threshold <- 0.4
novel_value <- 99L
write_species_rasters <- TRUE
reuse_existing <- !identical(
  tolower(Sys.getenv("ECOCHINA2_REUSE_POPULATION", "true")),
  "false"
)

reference_candidates <- c(
  file.path(base_dir, "raster", "ecosys_ori.tif"),
  file.path(base_dir, "data", "ecosys_ori.tif")
)
reference_file <- reference_candidates[file.exists(reference_candidates)][1]
if (is.na(reference_file)) stop("Missing data/ecosys_ori.tif.")
dual_root <- file.path(base_dir, "dual suit", "mf_var")
map_root <- file.path(base_dir, "result maps", "mf_var")
output_root <- file.path(base_dir, "future tree niche dual suitability var")
table_dir <- file.path(output_root, "tables")
map_output_root <- file.path(base_dir, "future tree niche var")
map_table_dir <- file.path(map_output_root, "tables")
raster_root <- file.path(output_root, "mf_var")

for (directory in c(table_dir, map_table_dir, raster_root)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}


# 1. Helpers ===================================================================

first_existing <- function(paths, label) {
  found <- paths[file.exists(paths)]
  if (!length(found)) stop("Missing ", label, ":\n", paste(paths, collapse = "\n"))
  found[[1]]
}

population_file <- first_existing(
  c(
    file.path(base_dir, "species_zone_population_long.csv"),
    file.path(base_dir, "data", "species_zone_population_long.csv"),
    file.path(base_dir, "data", "processed", "species_zone_population_long.csv"),
    # Read the normalized legacy lookup only when the script-1.2 source table
    # is unavailable. Never let this script's own output hide a newer source.
    file.path(base_dir, "future tree niche var", "tables", "population_projection_lookup_var.csv")
  ),
  "population lookup"
)

dual_file <- function(scenario, zone) {
  first_existing(
    c(
      file.path(dual_root, scenario, paste0("dual_suitability_zone", zone, ".tif")),
      file.path(dual_root, scenario, paste0("dual_zone", zone, ".tif"))
    ),
    paste("dual suitability", scenario, "zone", zone)
  )
}

assigned_map_file <- function(scenario, required = TRUE) {
  candidates <- c(
    file.path(
      map_root,
      paste0(
        "assigned_zone_", scenario,
        "_threshold0.4_tol1e-04_novel99_maskNA8_noNovelNormal.tif"
      )
    ),
    file.path(map_root, paste0("assigned_zone_", scenario, ".tif"))
  )
  found <- candidates[file.exists(candidates)]
  if (length(found)) return(found[[1]])
  if (required) stop("Missing assigned map for ", scenario)
  NA_character_
}

scenario_fields <- function(scenario) {
  if (scenario == "normal") {
    return(list(period = "1961-1990", ssp = "Reference"))
  }
  list(
    period = sub("SSP.*$", "", scenario),
    ssp = sub("^.*(SSP[0-9]+)$", "\\1", scenario)
  )
}

global_sum <- function(r) {
  value <- global(r, "sum", na.rm = TRUE)[1, 1]
  if (is.na(value)) 0 else as.numeric(value)
}

global_mean <- function(r) {
  value <- global(r, "mean", na.rm = TRUE)[1, 1]
  if (is.nan(value)) NA_real_ else as.numeric(value)
}

check_geometry <- function(x, template, label) {
  if (!compareGeom(x, template, stopOnError = FALSE)) {
    stop(label, " does not match data/ecosys_ori.tif geometry.")
  }
  invisible(TRUE)
}

raster_valid <- function(file, template, expected_layers = NULL) {
  if (!file.exists(file)) return(FALSE)
  tryCatch({
    x <- rast(file)
    layer_ok <- is.null(expected_layers) || nlyr(x) == expected_layers
    layer_ok && compareGeom(x, template, stopOnError = FALSE)
  }, error = function(e) FALSE)
}

files_current <- function(outputs, inputs) {
  outputs <- outputs[!is.na(outputs) & nzchar(outputs)]
  if (!length(outputs) || !all(file.exists(outputs))) return(FALSE)
  inputs <- inputs[!is.na(inputs) & nzchar(inputs) & file.exists(inputs)]
  if (!length(inputs)) return(TRUE)
  isTRUE(min(file.info(outputs)$mtime) >= max(file.info(inputs)$mtime))
}

write_table_if_changed <- function(table, file) {
  unchanged <- if (file.exists(file)) {
    previous <- tryCatch(fread(file), error = function(e) NULL)
    if (is.null(previous) || !setequal(names(previous), names(table)) ||
        nrow(previous) != nrow(table)) {
      FALSE
    } else {
      current <- copy(table)
      setcolorder(previous, names(current))
      if ("PopulationID" %in% names(current)) {
        setorder(previous, PopulationID)
        setorder(current, PopulationID)
      }
      isTRUE(all.equal(
        as.data.frame(previous), as.data.frame(current),
        check.attributes = FALSE
      ))
    }
  } else {
    FALSE
  }
  if (!unchanged) fwrite(table, file)
  invisible(!unchanged)
}

population_rank_cell <- function(values, source_zones) {
  n_population <- length(source_zones)
  ranked_zone <- rep(NA_real_, n_population)
  ranked_suitability <- rep(NA_real_, n_population)
  valid <- which(is.finite(values) & values > 0)
  
  if (length(valid)) {
    ordered <- valid[order(-values[valid], source_zones[valid])]
    ranked_zone[seq_along(ordered)] <- source_zones[ordered]
    ranked_suitability[seq_along(ordered)] <- values[ordered]
  }
  
  margin <- if (length(valid) >= 2L) {
    ranked_suitability[1] - ranked_suitability[2]
  } else {
    NA_real_
  }
  
  c(
    ranked_zone,
    ranked_suitability,
    n_pop_ranked = length(valid),
    n_pop_suitable = sum(is.finite(values) & values >= dual_threshold),
    top1_minus_top2 = margin
  )
}

species_output_paths <- function(scenario, species_name) {
  safe_species <- gsub("[^A-Za-z0-9_]", "_", species_name)
  safe_species <- gsub("^_+|_+$", "", safe_species)
  scenario_dir <- file.path(raster_root, scenario)
  
  list(
    species = file.path(
      scenario_dir, "species niche",
      paste0(safe_species, "_species_dual_suitability.tif")
    ),
    binary = file.path(
      scenario_dir, "species niche binary",
      paste0(safe_species, "_species_dual_binary_threshold", dual_threshold, ".tif")
    ),
    rank_zone = file.path(
      scenario_dir, "population ranking",
      paste0(safe_species, "_ranked_population_zones.tif")
    ),
    rank_suitability = file.path(
      scenario_dir, "population ranking",
      paste0(safe_species, "_ranked_population_suitability.tif")
    ),
    rank_summary = file.path(
      scenario_dir, "population ranking",
      paste0(safe_species, "_ranked_population_summary.tif")
    ),
    rank_index = file.path(
      scenario_dir, "population ranking",
      paste0(safe_species, "_population_rank_index.csv")
    )
  )
}

zone_area_table <- function(zone_map, cell_area, scenario) {
  area <- as.data.table(zonal(cell_area, zone_map, fun = "sum", na.rm = TRUE))
  if (!nrow(area)) return(data.table())
  setnames(area, names(area)[1:2], c("zoneID", "area_km2"))
  fields <- scenario_fields(scenario)
  area[, `:=`(
    zoneID = as.integer(zoneID),
    scenario = scenario,
    period = fields$period,
    ssp = fields$ssp
  )]
  area[]
}


# 2. Population lookup ==========================================================

population <- fread(population_file)

required_columns <- c("PopulationID", "Species")
if (!all(required_columns %in% names(population))) {
  stop("Population lookup must contain: ", paste(required_columns, collapse = ", "))
}

if (!"Population" %in% names(population) &&
    !"reference_abundance" %in% names(population)) {
  stop("Population lookup needs Population or reference_abundance.")
}

if (!"source_zone" %in% names(population)) {
  if (!"Zone" %in% names(population)) stop("Population lookup needs Zone or source_zone.")
  population[, source_zone := as.integer(gsub("[^0-9]", "", as.character(Zone)))]
}

# Recreate the original lookup schema. In particular, keep `projected` as the
# established column name so a first clean-branch run does not rewrite a valid
# legacy checkpoint merely because an interim draft used `projection_eligible`.
status_columns <- intersect(c("projected", "projection_eligible"), names(population))
if (length(status_columns)) population[, (status_columns) := NULL]

population[, `:=`(
  PopulationID = as.character(PopulationID),
  Species = as.character(Species),
  source_zone = as.integer(source_zone),
  reference_abundance = if ("Population" %in% names(population)) {
    as.numeric(Population)
  } else {
    as.numeric(reference_abundance)
  }
)]

# Keep the metadata columns used by the original future-niche tables, without
# making the expensive raster work depend on a palette file.
palette_candidates <- c(
  file.path(base_dir, "color_palette_China.csv"),
  file.path(base_dir, "data", "zone_palette.csv")
)
palette_file <- palette_candidates[file.exists(palette_candidates)][1]
if (!is.na(palette_file)) {
  palette <- fread(palette_file)
  if (all(c("zoneID", "zone", "COLOR") %in% names(palette))) {
    existing_metadata <- intersect(c("zone_name", "COLOR"), names(population))
    if (length(existing_metadata)) population[, (existing_metadata) := NULL]
    population <- merge(
      population,
      palette[, .(
        source_zone = as.integer(zoneID),
        zone_name = as.character(zone),
        COLOR = as.character(COLOR)
      )],
      by = "source_zone",
      all.x = TRUE,
      sort = FALSE
    )
  }
}
if (!"zone_name" %in% names(population)) population[, zone_name := NA_character_]
if (!"COLOR" %in% names(population)) population[, COLOR := NA_character_]
population[, projected := source_zone %in% modeled_zones]

if (anyDuplicated(population[, .(Species, source_zone)])) {
  stop("Duplicated species x source-zone population proxies were found.")
}

projected_population <- population[projected == TRUE]
if (!nrow(projected_population)) stop("No projection-eligible populations were found.")

population_lookup_output <- file.path(
  map_table_dir, "population_projection_lookup_var.csv"
)
population_exclusion_output <- file.path(
  map_table_dir, "population_projection_exclusions_var.csv"
)
dual_population_lookup_output <- file.path(
  table_dir, "dual_population_projection_lookup_var.csv"
)
dual_population_exclusion_output <- file.path(
  table_dir, "dual_population_projection_exclusions_var.csv"
)

write_table_if_changed(
  population,
  population_lookup_output
)
write_table_if_changed(
  population[projected == FALSE],
  population_exclusion_output
)
write_table_if_changed(
  projected_population,
  dual_population_lookup_output
)
write_table_if_changed(
  population[projected == FALSE],
  dual_population_exclusion_output
)

# Each legacy raster family depends on the lookup written immediately before it
# in the original workflow. Keeping these two checkpoints separate preserves
# their historical mtimes on the first clean-branch run.
assigned_population_dependency <- population_lookup_output
dual_population_dependency <- dual_population_lookup_output


# 3. Reference raster and cell area ============================================

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
assigned_files <- setNames(
  vapply(scenarios, assigned_map_file, character(1), required = FALSE),
  scenarios
)


# 4. Ecosystem assigned-area summaries =========================================

ecosystem_area_file <- file.path(map_table_dir, "future_ecosystem_area_var.csv")

if (reuse_existing && files_current(ecosystem_area_file, assigned_files)) {
  ecosystem_area <- fread(ecosystem_area_file)
  if ("method" %in% names(ecosystem_area)) {
    ecosystem_area <- ecosystem_area[method == "mf_var"]
  } else {
    ecosystem_area[, method := "mf_var"]
  }
  cat("[REUSE TABLE]", ecosystem_area_file, "\n")
} else {
  ecosystem_area_results <- list()
  
  for (scenario in scenarios) {
    map_file <- assigned_map_file(scenario, required = FALSE)
    if (is.na(map_file)) next
    
    assigned <- rast(map_file)
    check_geometry(assigned, reference, paste("Assigned map", scenario))
    ecosystem_area_results[[length(ecosystem_area_results) + 1L]] <-
      zone_area_table(assigned, get_cell_area(), scenario)
    rm(assigned)
  }
  
  ecosystem_area <- rbindlist(ecosystem_area_results, fill = TRUE)
  ecosystem_area[, method := "mf_var"]
  fwrite(ecosystem_area, ecosystem_area_file)
}


# 5. Normal-to-future assigned-zone changes ====================================

transition_file <- file.path(map_table_dir, "future_ecosystem_transition_var.csv")

if (reuse_existing && files_current(transition_file, assigned_files)) {
  transitions <- fread(transition_file)
  if ("method" %in% names(transitions)) {
    transitions <- transitions[method == "mf_var"]
  } else {
    transitions[, method := "mf_var"]
  }
  cat("[REUSE TABLE]", transition_file, "\n")
} else {
  transition_results <- list()
  normal_file <- assigned_map_file("normal", required = FALSE)
  
  if (!is.na(normal_file)) {
    normal_map <- rast(normal_file)
    check_geometry(normal_map, reference, "Normal assigned map")
    
    for (scenario in future_scenarios) {
      future_file <- assigned_map_file(scenario, required = FALSE)
      if (is.na(future_file)) next
      
      future_map <- rast(future_file)
      check_geometry(future_map, reference, paste("Future assigned map", scenario))
      transition_code <- ifel(
        is.na(normal_map) | is.na(future_map),
        NA,
        normal_map * 1000L + future_map
      )
      transition <- as.data.table(
        zonal(get_cell_area(), transition_code, fun = "sum", na.rm = TRUE)
      )
      
      if (nrow(transition)) {
        setnames(transition, names(transition)[1:2], c("transition_code", "area_km2"))
        fields <- scenario_fields(scenario)
        transition[, `:=`(
          normal_zone = as.integer(transition_code %/% 1000L),
          future_zone = as.integer(transition_code %% 1000L),
          scenario = scenario,
          period = fields$period,
          ssp = fields$ssp
        )]
        transition[, change_class := fifelse(
          future_zone == novel_value,
          "novel",
          fifelse(normal_zone == future_zone, "stable", "changed")
        )]
        transition_results[[length(transition_results) + 1L]] <- transition
      }
      
      rm(future_map, transition_code, transition)
    }
    rm(normal_map)
  }
  
  transitions <- rbindlist(transition_results, fill = TRUE)
  transitions[, method := "mf_var"]
  fwrite(transitions, transition_file)
}


# 6. Assigned-map population and species niches ===============================

# Preserve the original, inexpensive categorical products in their established
# paths. Complete tables plus existing raster files form a checkpoint, so an
# ordinary rerun does not even open the six future assigned maps.
assigned_population_file <- file.path(
  map_table_dir, "future_population_niche_area_var.csv"
)
assigned_species_file <- file.path(
  map_table_dir, "future_species_niche_area_var.csv"
)

read_mf_cache <- function(file) {
  empty <- data.table(
    method = character(), scenario = character(), Species = character(),
    PopulationID = character(), population_raster = character(),
    species_raster = character()
  )
  if (!reuse_existing || !file.exists(file)) return(empty)
  result <- fread(file)
  if (!all(c("scenario", "Species") %in% names(result))) return(empty)
  if ("method" %in% names(result)) {
    result <- result[method == "mf_var"]
  } else {
    result[, method := "mf_var"]
  }
  result
}

# Use explicit column vectors for filters involving loop/function arguments.
# The data.table `..name` shortcut belongs in j; in i it can be looked up as
# a literal missing object instead of the surrounding scenario/species value.
dual_cache_rows <- function(cache, scenario_value, species_value = NULL) {
  if (!all(c("method", "scenario", "Species") %in% names(cache))) {
    return(cache[0L])
  }
  keep <- cache$method == "mf_var" & cache$scenario == scenario_value
  if (!is.null(species_value)) keep <- keep & cache$Species == species_value
  row_index <- which(keep)
  cache[row_index]
}

dual_cache_identity_valid <- function(population_rows, species_rows,
                                      expected_population) {
  population_columns <- c(
    "PopulationID", "source_zone", "suitable_area_km2",
    "suitability_weighted_area_km2", "mean_dual_suitability"
  )
  species_columns <- c(
    "suitable_area_km2", "suitability_weighted_area_km2", "mean_dual_suitability"
  )
  if (!all(population_columns %in% names(population_rows)) ||
      !all(species_columns %in% names(species_rows)) ||
      nrow(population_rows) != nrow(expected_population) ||
      nrow(species_rows) != 1L ||
      anyDuplicated(population_rows$PopulationID) ||
      !setequal(population_rows$PopulationID, expected_population$PopulationID)) {
    return(FALSE)
  }
  expected_order <- match(population_rows$PopulationID, expected_population$PopulationID)
  isTRUE(all(population_rows$source_zone ==
               expected_population$source_zone[expected_order]))
}

cached_assigned_population <- read_mf_cache(assigned_population_file)
cached_assigned_species <- read_mf_cache(assigned_species_file)
expected_assigned_population_rows <-
  nrow(projected_population) * length(future_scenarios)
expected_assigned_species_rows <-
  uniqueN(projected_population$Species) * length(future_scenarios)

assigned_scenario_cache_valid <- function(population_rows, species_rows,
                                          assigned_input) {
  species_names <- sort(unique(projected_population$Species))
  required_population <- c("PopulationID", "Species", "population_raster")
  required_species <- c("Species", "species_raster")
  if (!file.exists(assigned_input) ||
      !all(required_population %in% names(population_rows)) ||
      !all(required_species %in% names(species_rows)) ||
      nrow(population_rows) != nrow(projected_population) ||
      nrow(species_rows) != length(species_names) ||
      !setequal(population_rows$PopulationID, projected_population$PopulationID) ||
      !setequal(species_rows$Species, species_names) ||
      anyDuplicated(population_rows$PopulationID) ||
      anyDuplicated(species_rows$Species)) {
    return(FALSE)
  }
  
  population_rasters <- unique(as.character(population_rows$population_raster))
  species_rasters <- unique(as.character(species_rows$species_raster))
  if (anyNA(c(population_rasters, species_rasters)) ||
      any(!nzchar(c(population_rasters, species_rasters))) ||
      !all(vapply(
        population_rasters, raster_valid, logical(1),
        template = reference, expected_layers = 1L
      )) ||
      !all(vapply(
        species_rasters, raster_valid, logical(1),
        template = reference, expected_layers = 1L
      )) ||
      !all(vapply(
        population_rasters, files_current, logical(1),
        inputs = c(assigned_input, assigned_population_dependency)
      ))) {
    return(FALSE)
  }
  
  all(vapply(species_names, function(species_name) {
    species_raster <- unique(as.character(
      species_rows[Species == species_name]$species_raster
    ))
    population_raster <- unique(as.character(
      population_rows[Species == species_name]$population_raster
    ))
    length(species_raster) == 1L && length(population_raster) == 1L &&
      files_current(
        species_raster,
        c(assigned_input, assigned_population_dependency, population_raster)
      )
  }, logical(1)))
}

assigned_tables_complete <-
  nrow(cached_assigned_population) == expected_assigned_population_rows &&
  nrow(cached_assigned_species) == expected_assigned_species_rows &&
  all(future_scenarios %in% cached_assigned_population$scenario) &&
  all(future_scenarios %in% cached_assigned_species$scenario) &&
  files_current(
    c(assigned_population_file, assigned_species_file),
    c(assigned_files[future_scenarios], assigned_population_dependency)
  ) &&
  all(vapply(future_scenarios, function(scenario_value) {
    assigned_scenario_cache_valid(
      cached_assigned_population[scenario == scenario_value],
      cached_assigned_species[scenario == scenario_value],
      assigned_files[[scenario_value]]
    )
  }, logical(1)))

if (assigned_tables_complete) {
  assigned_population_area <- cached_assigned_population
  assigned_species_area <- cached_assigned_species
  cat("[REUSE TABLES] Assigned-map population/species niches are complete.\n")
} else {
  assigned_population_results <- list()
  assigned_species_results <- list()
  species_names <- sort(unique(projected_population$Species))
  
  for (scenario in future_scenarios) {
    scenario_key <- scenario
    cached_scenario_population <- cached_assigned_population[scenario == scenario_key]
    cached_scenario_species <- cached_assigned_species[scenario == scenario_key]
    scenario_complete <-
      files_current(
        c(assigned_population_file, assigned_species_file),
        c(assigned_files[[scenario]], assigned_population_dependency)
      ) &&
      assigned_scenario_cache_valid(
        cached_scenario_population,
        cached_scenario_species,
        assigned_files[[scenario]]
      )
    
    if (reuse_existing && scenario_complete) {
      assigned_population_results[[length(assigned_population_results) + 1L]] <-
        cached_scenario_population
      assigned_species_results[[length(assigned_species_results) + 1L]] <-
        cached_scenario_species
      cat("[REUSE ASSIGNED SCENARIO]", scenario, "\n")
      next
    }
    
    future_map <- rast(assigned_map_file(scenario))
    check_geometry(future_map, reference, paste("Assigned map", scenario))
    fields <- scenario_fields(scenario)
    population_dir <- file.path(map_output_root, "mf_var", scenario, "population niche")
    species_dir <- file.path(map_output_root, "mf_var", scenario, "species niche")
    dir.create(population_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (species_name in species_names) {
      species_population <- projected_population[
        Species == species_name
      ][order(source_zone)]
      population_output <- file.path(
        population_dir, paste0(species_name, "_population_niche.tif")
      )
      species_output <- file.path(
        species_dir, paste0(species_name, "_species_niche.tif")
      )
      
      cached_population_rows <- cached_scenario_population[Species == species_name]
      cached_species_rows <- cached_scenario_species[Species == species_name]
      species_checkpoint <-
        nrow(cached_population_rows) == nrow(species_population) &&
        nrow(cached_species_rows) == 1L &&
        setequal(
          cached_population_rows$PopulationID,
          species_population$PopulationID
        ) &&
        files_current(
          c(assigned_population_file, assigned_species_file),
          c(assigned_files[[scenario]], assigned_population_dependency)
        ) &&
        raster_valid(population_output, reference, 1L) &&
        raster_valid(species_output, reference, 1L) &&
        files_current(
          population_output,
          c(assigned_files[[scenario]], assigned_population_dependency)
        ) &&
        files_current(
          species_output,
          c(
            assigned_files[[scenario]], assigned_population_dependency,
            population_output
          )
        )
      if (reuse_existing && species_checkpoint) {
        assigned_population_results[[length(assigned_population_results) + 1L]] <-
          cached_population_rows
        assigned_species_results[[length(assigned_species_results) + 1L]] <-
          cached_species_rows
        cat("[REUSE ASSIGNED]", species_name, "|", scenario, "\n")
        next
      }
      
      population_output_current <-
        raster_valid(population_output, reference, 1L) &&
        files_current(
          population_output,
          c(assigned_files[[scenario]], assigned_population_dependency)
        )
      if (reuse_existing && population_output_current) {
        population_map <- rast(population_output)
      } else {
        population_map <- subst(
          future_map,
          from = species_population$source_zone,
          to = species_population$source_zone,
          others = NA
        )
        names(population_map) <- "source_zone"
        writeRaster(
          population_map, population_output, overwrite = TRUE,
          datatype = "INT2S", NAflag = -32768,
          wopt = list(gdal = "COMPRESS=LZW")
        )
      }
      
      species_output_current <-
        raster_valid(species_output, reference, 1L) &&
        files_current(
          species_output,
          c(
            assigned_files[[scenario]], assigned_population_dependency,
            population_output
          )
        )
      if (reuse_existing && species_output_current) {
        species_map <- rast(species_output)
      } else {
        species_map <- ifel(!is.na(population_map), 1L, NA)
        names(species_map) <- "species_suitable"
        writeRaster(
          species_map, species_output, overwrite = TRUE,
          datatype = "INT1U", wopt = list(gdal = "COMPRESS=LZW")
        )
      }
      
      population_area <- as.data.table(
        zonal(get_cell_area(), population_map, fun = "sum", na.rm = TRUE)
      )
      if (nrow(population_area)) {
        setnames(
          population_area,
          names(population_area)[1:2],
          c("source_zone", "future_area_km2")
        )
        population_area[, source_zone := as.integer(source_zone)]
      } else {
        population_area <- data.table(
          source_zone = integer(), future_area_km2 = numeric()
        )
      }
      population_area <- merge(
        species_population[, .(
          PopulationID, Species, source_zone, zone_name, COLOR,
          reference_abundance
        )],
        population_area,
        by = "source_zone",
        all.x = TRUE,
        sort = FALSE
      )
      population_area[is.na(future_area_km2), future_area_km2 := 0]
      population_area[, `:=`(
        method = "mf_var", scenario = scenario,
        period = fields$period, ssp = fields$ssp,
        population_raster = population_output
      )]
      assigned_population_results[[length(assigned_population_results) + 1L]] <-
        population_area
      
      assigned_species_results[[length(assigned_species_results) + 1L]] <- data.table(
        Species = species_name,
        method = "mf_var",
        scenario = scenario,
        period = fields$period,
        ssp = fields$ssp,
        populations_projected = nrow(species_population),
        future_area_km2 = global_sum(get_cell_area() * species_map),
        species_raster = species_output
      )
      
      rm(population_map, species_map, population_area)
      gc()
    }
    
    rm(future_map)
    gc()
  }
  
  assigned_population_area <- rbindlist(assigned_population_results, fill = TRUE)
  assigned_species_area <- rbindlist(assigned_species_results, fill = TRUE)
  if (nrow(assigned_population_area)) {
    setorder(assigned_population_area, Species, source_zone, ssp, period)
  }
  if (nrow(assigned_species_area)) {
    setorder(assigned_species_area, Species, ssp, period)
  }
  fwrite(assigned_population_area, assigned_population_file)
  fwrite(assigned_species_area, assigned_species_file)
}


# 7. Population and species dual-suitability niches ============================

population_results <- list()
species_results <- list()
population_layer_index <- list()
species_raster_index <- list()
population_rank_index <- list()

population_table_file <- file.path(
  table_dir,
  "dual_population_niche_area_var.csv"
)
species_table_file <- file.path(
  table_dir,
  "dual_species_niche_area_var.csv"
)

cached_population <- if (reuse_existing && file.exists(population_table_file)) {
  fread(population_table_file)
} else {
  data.table(method = character(), scenario = character(), Species = character())
}
cached_species <- if (reuse_existing && file.exists(species_table_file)) {
  fread(species_table_file)
} else {
  data.table(method = character(), scenario = character(), Species = character())
}

if (nrow(cached_population) && !"method" %in% names(cached_population)) {
  cached_population[, method := "mf_var"]
}
if (nrow(cached_species) && !"method" %in% names(cached_species)) {
  cached_species[, method := "mf_var"]
}

dual_checkpoint_dir <- file.path(table_dir, "cache", "dual_scenarios")
dir.create(dual_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

for (scenario in scenarios) {
  cat("\n[POPULATION/SPECIES]", scenario, "\n")
  
  dual_files <- vapply(modeled_zones, function(z) dual_file(scenario, z), character(1))
  species_names <- sort(unique(projected_population$Species))
  scenario_cached_population <- cached_population
  scenario_cached_species <- cached_species
  summary_cache_files <- c(population_table_file, species_table_file)
  checkpoint_file <- file.path(dual_checkpoint_dir, paste0(scenario, ".rds"))
  
  # Preserve completed scenario summaries even if a later scenario fails before
  # the canonical combined CSV files can be written. Only this small checkpoint
  # is new; all model and raster output paths remain the established paths.
  if (reuse_existing &&
      files_current(checkpoint_file, c(dual_files, dual_population_dependency)) &&
      !files_current(summary_cache_files, checkpoint_file)) {
    checkpoint <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
    if (is.list(checkpoint) &&
        identical(checkpoint$scenario, scenario) &&
        identical(checkpoint$dual_threshold, dual_threshold) &&
        is.data.frame(checkpoint$population) &&
        is.data.frame(checkpoint$species)) {
      scenario_cached_population <- as.data.table(checkpoint$population)
      scenario_cached_species <- as.data.table(checkpoint$species)
      summary_cache_files <- checkpoint_file
      cat("[RESUME SUMMARY CHECKPOINT]", scenario, "\n")
    }
  }
  
  scenario_reusable <- reuse_existing && all(vapply(
    species_names,
    function(species_name) {
      species_population <- projected_population[Species == species_name]
      paths <- species_output_paths(scenario, species_name)
      cached_population_rows <- dual_cache_rows(
        scenario_cached_population, scenario, species_name
      )
      cached_species_rows <- dual_cache_rows(
        scenario_cached_species, scenario, species_name
      )
      
      dual_cache_identity_valid(
        cached_population_rows, cached_species_rows, species_population
      ) &&
        raster_valid(paths$species, reference, 1L) &&
        raster_valid(paths$binary, reference, 1L) &&
        raster_valid(paths$rank_zone, reference, nrow(species_population)) &&
        raster_valid(paths$rank_suitability, reference, nrow(species_population)) &&
        raster_valid(paths$rank_summary, reference, 3L) &&
        files_current(
          summary_cache_files,
          c(
            dual_files[match(species_population$source_zone, modeled_zones)],
            dual_population_dependency, paths$species
          )
        ) &&
        files_current(
          paths$binary,
          c(paths$species, dual_files[match(species_population$source_zone, modeled_zones)],
            dual_population_dependency)
        ) &&
        files_current(
          unlist(paths[c(
            "species", "binary", "rank_zone", "rank_suitability", "rank_summary"
          )]),
          c(
            dual_files[match(species_population$source_zone, modeled_zones)],
            dual_population_dependency
          )
        )
    },
    logical(1)
  ))
  
  if (scenario_reusable) {
    cat("[REUSE COMPLETE SCENARIO]", scenario, "\n")
    for (species_name in species_names) {
      species_population <- projected_population[Species == species_name]
      paths <- species_output_paths(scenario, species_name)
      source_zones <- species_population$source_zone
      population_results[[length(population_results) + 1L]] <- dual_cache_rows(
        scenario_cached_population, scenario, species_name
      )
      species_results[[length(species_results) + 1L]] <- dual_cache_rows(
        scenario_cached_species, scenario, species_name
      )
      population_layer_index[[length(population_layer_index) + 1L]] <- data.table(
        Species = species_name,
        PopulationID = species_population$PopulationID,
        source_zone = source_zones,
        method = "mf_var",
        scenario = scenario,
        layer = seq_along(source_zones),
        layer_name = species_population$PopulationID,
        source_dual_raster = dual_files[match(source_zones, modeled_zones)]
      )
      population_rank_index[[length(population_rank_index) + 1L]] <- data.table(
        Species = species_name,
        PopulationID = species_population$PopulationID,
        source_zone = source_zones,
        method = "mf_var",
        scenario = scenario,
        ranked_zone_raster = paths$rank_zone,
        ranked_suitability_raster = paths$rank_suitability,
        ranked_summary_raster = paths$rank_summary
      )
      species_raster_index[[length(species_raster_index) + 1L]] <- data.table(
        Species = species_name,
        method = "mf_var",
        scenario = scenario,
        species_raster = paths$species,
        species_binary_raster = paths$binary,
        ranked_zone_raster = paths$rank_zone,
        ranked_suitability_raster = paths$rank_suitability,
        ranked_summary_raster = paths$rank_summary
      )
    }
    next
  }
  
  dual_stack <- rast(dual_files)
  names(dual_stack) <- paste0("zone", modeled_zones)
  check_geometry(dual_stack, reference, paste("Dual suitability", scenario))
  
  fields <- scenario_fields(scenario)
  scenario_dir <- file.path(raster_root, scenario)
  species_dir <- file.path(scenario_dir, "species niche")
  species_binary_dir <- file.path(scenario_dir, "species niche binary")
  rank_dir <- file.path(scenario_dir, "population ranking")
  for (directory in c(species_dir, species_binary_dir, rank_dir)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  }
  
  for (species_name in species_names) {
    species_population <- projected_population[Species == species_name]
    source_zones <- species_population$source_zone
    layer_index <- match(source_zones, modeled_zones)
    population_stack <- dual_stack[[layer_index]]
    names(population_stack) <- species_population$PopulationID
    
    paths <- species_output_paths(scenario, species_name)
    species_file <- paths$species
    species_binary_file <- paths$binary
    rank_zone_file <- paths$rank_zone
    rank_suitability_file <- paths$rank_suitability
    rank_summary_file <- paths$rank_summary
    rank_index_file <- paths$rank_index
    
    cached_population_rows <- dual_cache_rows(
      scenario_cached_population, scenario, species_name
    )
    cached_species_rows <- dual_cache_rows(
      scenario_cached_species, scenario, species_name
    )
    cached_tables_valid <-
      dual_cache_identity_valid(
        cached_population_rows, cached_species_rows, species_population
      ) &&
      files_current(
        summary_cache_files,
        c(dual_files[layer_index], dual_population_dependency)
      )
    
    species_file_current <-
      raster_valid(species_file, reference, 1L) &&
      files_current(species_file, c(dual_files[layer_index], dual_population_dependency))
    binary_file_current <-
      raster_valid(species_binary_file, reference, 1L) &&
      files_current(
        species_binary_file,
        c(species_file, dual_files[layer_index], dual_population_dependency)
      )
    rank_valid <-
      raster_valid(rank_zone_file, reference, length(source_zones)) &&
      raster_valid(rank_suitability_file, reference, length(source_zones)) &&
      raster_valid(rank_summary_file, reference, 3L) &&
      files_current(
        c(rank_zone_file, rank_suitability_file, rank_summary_file),
        c(dual_files[layer_index], dual_population_dependency)
      )
    species_rasters_valid <-
      species_file_current && binary_file_current && rank_valid
    
    if (reuse_existing && cached_tables_valid && species_rasters_valid &&
        files_current(summary_cache_files, species_file)) {
      cat("[REUSE]", species_name, "|", scenario, "\n")
      population_results[[length(population_results) + 1L]] <- cached_population_rows
      species_results[[length(species_results) + 1L]] <- cached_species_rows
      population_layer_index[[length(population_layer_index) + 1L]] <- data.table(
        Species = species_name,
        PopulationID = species_population$PopulationID,
        source_zone = source_zones,
        method = "mf_var",
        scenario = scenario,
        layer = seq_along(source_zones),
        layer_name = names(population_stack),
        source_dual_raster = dual_files[layer_index]
      )
      population_rank_index[[length(population_rank_index) + 1L]] <- data.table(
        Species = species_name,
        PopulationID = species_population$PopulationID,
        source_zone = source_zones,
        method = "mf_var",
        scenario = scenario,
        ranked_zone_raster = rank_zone_file,
        ranked_suitability_raster = rank_suitability_file,
        ranked_summary_raster = rank_summary_file
      )
      species_raster_index[[length(species_raster_index) + 1L]] <- data.table(
        Species = species_name,
        method = "mf_var",
        scenario = scenario,
        species_raster = species_file,
        species_binary_raster = species_binary_file,
        ranked_zone_raster = rank_zone_file,
        ranked_suitability_raster = rank_suitability_file,
        ranked_summary_raster = rank_summary_file
      )
      rm(population_stack)
      next
    }
    
    # A missing rank/binary file does not invalidate unchanged population areas.
    # Preserve those summaries while repairing only the missing raster product.
    if (reuse_existing && cached_tables_valid) {
      population_results[[length(population_results) + 1L]] <- cached_population_rows
      cat("[REUSE POPULATION SUMMARY]", species_name, "|", scenario, "\n")
    }
    population_indices <- if (reuse_existing && cached_tables_valid) {
      integer()
    } else {
      seq_len(nrow(species_population))
    }
    for (population_index in population_indices) {
      suitability <- population_stack[[population_index]]
      suitable <- ifel(is.na(suitability), NA, suitability >= dual_threshold)
      
      area_km2 <- global_sum(ifel(suitable, get_cell_area(), 0))
      weighted_area_km2 <- global_sum(ifel(
        is.na(suitability),
        NA,
        get_cell_area() * suitability
      ))
      
      population_results[[length(population_results) + 1L]] <- data.table(
        PopulationID = species_population$PopulationID[[population_index]],
        Species = species_name,
        source_zone = source_zones[[population_index]],
        reference_abundance = species_population$reference_abundance[[population_index]],
        method = "mf_var",
        scenario = scenario,
        period = fields$period,
        ssp = fields$ssp,
        suitable_area_km2 = area_km2,
        suitability_weighted_area_km2 = weighted_area_km2,
        mean_dual_suitability = global_mean(suitability),
        dual_raster = dual_files[[layer_index[[population_index]]]]
      )
    }
    
    if (reuse_existing && species_file_current) {
      species_suitability <- rast(species_file)
      names(species_suitability) <- "species_dual_suitability"
    } else {
      species_suitability <- app(
        population_stack,
        fun = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE),
        filename = species_file,
        overwrite = TRUE,
        wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW")
      )
      names(species_suitability) <- "species_dual_suitability"
    }
    
    species_summary_current <- reuse_existing && cached_tables_valid &&
      species_file_current && files_current(summary_cache_files, species_file)
    binary_file_current <-
      raster_valid(species_binary_file, reference, 1L) &&
      files_current(
        species_binary_file,
        c(species_file, dual_files[layer_index], dual_population_dependency)
      )
    write_species_binary <- write_species_rasters &&
      !(reuse_existing && binary_file_current)
    species_suitable <- if (!species_summary_current || write_species_binary) {
      ifel(is.na(species_suitability), NA, species_suitability >= dual_threshold)
    } else {
      NULL
    }
    
    if (species_summary_current) {
      species_area_km2 <- cached_species_rows$suitable_area_km2[[1]]
      species_weighted_area_km2 <- cached_species_rows$suitability_weighted_area_km2[[1]]
      species_mean_suitability <- cached_species_rows$mean_dual_suitability[[1]]
      cat("[REUSE SPECIES SUMMARY]", species_name, "|", scenario, "\n")
    } else {
      species_area_km2 <- global_sum(ifel(species_suitable, get_cell_area(), 0))
      species_weighted_area_km2 <- global_sum(ifel(
        is.na(species_suitability),
        NA,
        get_cell_area() * species_suitability
      ))
      species_mean_suitability <- global_mean(species_suitability)
    }
    
    if (write_species_binary) {
      writeRaster(
        species_suitable,
        species_binary_file,
        overwrite = TRUE,
        datatype = "INT1U",
        wopt = list(gdal = "COMPRESS=LZW")
      )
    }
    
    if (!(reuse_existing && rank_valid)) {
      rank_all <- app(
        population_stack,
        fun = function(x) population_rank_cell(x, source_zones)
      )
      n_population <- length(source_zones)
      rank_zone <- rank_all[[seq_len(n_population)]]
      rank_suitability <- rank_all[[n_population + seq_len(n_population)]]
      rank_summary <- rank_all[[2L * n_population + 1:3]]
      names(rank_zone) <- paste0("rank", seq_len(n_population), "_zone")
      names(rank_suitability) <- paste0("rank", seq_len(n_population), "_suit")
      names(rank_summary) <- c("n_pop_ranked", "n_pop_suitable", "top1_minus_top2")
      writeRaster(
        rank_zone,
        rank_zone_file,
        overwrite = TRUE,
        datatype = "INT2S",
        wopt = list(gdal = "COMPRESS=LZW")
      )
      writeRaster(
        rank_suitability,
        rank_suitability_file,
        overwrite = TRUE,
        datatype = "FLT4S",
        wopt = list(gdal = "COMPRESS=LZW")
      )
      writeRaster(
        rank_summary,
        rank_summary_file,
        overwrite = TRUE,
        datatype = "FLT4S",
        wopt = list(gdal = "COMPRESS=LZW")
      )
      rm(rank_all, rank_zone, rank_suitability, rank_summary)
    }
    
    rank_index <- data.table(
      Species = species_name,
      PopulationID = species_population$PopulationID,
      source_zone = source_zones,
      method = "mf_var",
      scenario = scenario,
      ranked_zone_raster = rank_zone_file,
      ranked_suitability_raster = rank_suitability_file,
      ranked_summary_raster = rank_summary_file
    )
    fwrite(rank_index, rank_index_file)
    population_rank_index[[length(population_rank_index) + 1L]] <- rank_index
    population_layer_index[[length(population_layer_index) + 1L]] <- data.table(
      Species = species_name,
      PopulationID = species_population$PopulationID,
      source_zone = source_zones,
      method = "mf_var",
      scenario = scenario,
      layer = seq_along(source_zones),
      layer_name = names(population_stack),
      source_dual_raster = dual_files[layer_index]
    )
    species_raster_index[[length(species_raster_index) + 1L]] <- data.table(
      Species = species_name,
      method = "mf_var",
      scenario = scenario,
      species_raster = species_file,
      species_binary_raster = species_binary_file,
      ranked_zone_raster = rank_zone_file,
      ranked_suitability_raster = rank_suitability_file,
      ranked_summary_raster = rank_summary_file
    )
    
    species_results[[length(species_results) + 1L]] <- data.table(
      Species = species_name,
      n_populations = nrow(species_population),
      method = "mf_var",
      scenario = scenario,
      period = fields$period,
      ssp = fields$ssp,
      suitable_area_km2 = species_area_km2,
      suitability_weighted_area_km2 = species_weighted_area_km2,
      mean_dual_suitability = species_mean_suitability,
      species_raster = if (write_species_rasters) species_file else NA_character_,
      species_binary_raster = if (write_species_rasters) species_binary_file else NA_character_,
      ranked_zone_raster = rank_zone_file,
      ranked_suitability_raster = rank_suitability_file,
      ranked_summary_raster = rank_summary_file
    )
    
    rm(population_stack, species_suitability, species_suitable)
    gc()
  }
  
  saveRDS(
    list(
      scenario = scenario,
      dual_threshold = dual_threshold,
      population = dual_cache_rows(rbindlist(population_results, fill = TRUE), scenario),
      species = dual_cache_rows(rbindlist(species_results, fill = TRUE), scenario)
    ),
    checkpoint_file
  )
  cat("[SAVED SUMMARY CHECKPOINT]", scenario, "\n")
  
  rm(dual_stack)
  gc()
}

population_area <- rbindlist(population_results, fill = TRUE)
species_area <- rbindlist(species_results, fill = TRUE)
population_layers <- rbindlist(population_layer_index, fill = TRUE)
species_index <- rbindlist(species_raster_index, fill = TRUE)
population_rank_index_all <- rbindlist(population_rank_index, fill = TRUE)

setorder(population_area, Species, source_zone, scenario)
setorder(species_area, Species, scenario)

fwrite(population_area, population_table_file)
fwrite(species_area, species_table_file)
fwrite(
  population_layers,
  file.path(table_dir, "dual_population_layer_index_var.csv")
)
fwrite(
  species_index,
  file.path(table_dir, "dual_species_raster_index_var.csv")
)
fwrite(
  population_rank_index_all,
  file.path(table_dir, "dual_population_rank_index_var.csv")
)

settings <- data.table(
  setting = c(
    "workflow", "dual_threshold", "novel_value", "modeled_zones",
    "population_definition", "species_definition"
  ),
  value = c(
    "selected-variable Multi-Forest",
    dual_threshold,
    novel_value,
    paste(modeled_zones, collapse = ","),
    "species x source ecotype; >=10 occupied cells upstream",
    "pixel-wise maximum across projection-eligible populations"
  )
)
fwrite(settings, file.path(table_dir, "dual_niche_settings_var.csv"))

cat(
  "\nCOMPLETE\n",
  "Population rows:", nrow(population_area), "\n",
  "Species rows:", nrow(species_area), "\n",
  "Output:", output_root, "\n"
)
