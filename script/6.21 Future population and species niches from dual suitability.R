# 6.21 Future population and species niches from dual suitability
# Corrected climate models: rf_var and mf_var
# ==============================================================================
#
# Population niche:
#   The source-zone dual suitability raster for one species x source-zone
#   population.
#
# Species niche:
#   Pixel-wise maximum of all source-population dual suitability rasters.
#
# Threshold rule:
#   dual suitability >= 0.4 is suitable;
#   a pixel is novel in the assigned-zone workflow only when every zone has
#   dual suitability < 0.4.
#
# Ranking:
#   All positive dual suitability values are ranked. The 0.4 threshold is used
#   only for n_pop_suitable and binary/area summaries.
#
# Inputs:
#   dual suit/rf_var/{scenario}/dual_suitability_zone*.tif
#   dual suit/mf_var/{scenario}/dual_suitability_zone*.tif
#   future tree niche var/tables/population_projection_lookup_var.csv
#
# Outputs:
#   future tree niche dual suitability var/
#
# Existing script 6.2 outputs are not overwritten.
# ==============================================================================

library(terra)
library(data.table)

rm(list = ls())
gc()


# 0. Paths and settings ==========================================================

base_dir <- "H:/Jing/ecoChina2"

dual_root <- file.path(
  base_dir,
  "dual suit"
)

population_lookup_file <- file.path(
  base_dir,
  "future tree niche var",
  "tables",
  "population_projection_lookup_var.csv"
)

reference_file <- file.path(
  base_dir,
  "raster/ecosys_ori.tif"
)

output_root <- file.path(
  base_dir,
  "future tree niche dual suitability var"
)

table_dir <- file.path(
  output_root,
  "tables"
)

tmp_dir <- file.path(
  base_dir,
  "tmp_dual_population_niche_var"
)

dir.create(
  output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tmp_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

terraOptions(
  tempdir = tmp_dir,
  memfrac = 0.15
)

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

model_zoneID <- c(
  1:7,
  9:50,
  52:55
)

# Consistent with assign-zone: maximum < 0.4 is novel.
suitability_threshold <- 0.4

# Rank every positive dual suitability value.
rank_min_suitability <- 0

reuse_existing_outputs <- TRUE

# Normally FALSE because each population layer already exists as a source-zone
# dual suitability raster.
write_population_stack <- FALSE


# 1. Helpers ====================================================================

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


remove_raster_files <- function(filepath) {
  
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
  
  if (length(files) > 0) {
    gc()
    suppressWarnings(
      unlink(
        files,
        force = TRUE
      )
    )
  }
  
  invisible(TRUE)
}


dual_file <- function(
    method,
    scenario,
    zone) {
  
  file.path(
    dual_root,
    method,
    scenario,
    paste0(
      "dual_suitability_zone",
      zone,
      ".tif"
    )
  )
}


raster_valid <- function(
    raster_file,
    template,
    expected_layers = NULL) {
  
  if (!file.exists(raster_file)) {
    return(FALSE)
  }
  
  tryCatch({
    
    x <- rast(raster_file)
    
    if (!is.null(expected_layers)) {
      if (nlyr(x) != expected_layers) {
        return(FALSE)
      }
    }
    
    compareGeom(
      x,
      template,
      stopOnError = FALSE
    )
    
  }, error = function(e) {
    FALSE
  })
}


global_scalar <- function(
    raster,
    fun,
    zero_if_missing = FALSE) {
  
  value <- tryCatch(
    global(
      raster,
      fun = fun,
      na.rm = TRUE
    )[1, 1],
    error = function(e) {
      NA_real_
    }
  )
  
  value <- as.numeric(value)
  
  if (
    length(value) == 0 ||
    !is.finite(value)
  ) {
    
    if (zero_if_missing) {
      return(0)
    }
    
    return(NA_real_)
  }
  
  value
}


suitable_area <- function(
    suitability,
    cell_area,
    threshold) {
  
  # Equality is retained because only values strictly below 0.4 are novel.
  area_raster <- ifel(
    !is.na(suitability) &
      suitability >= threshold,
    cell_area,
    NA
  )
  
  global_scalar(
    area_raster,
    "sum",
    zero_if_missing = TRUE
  )
}


weighted_area <- function(
    suitability,
    cell_area) {
  
  global_scalar(
    suitability * cell_area,
    "sum",
    zero_if_missing = TRUE
  )
}


mean_suitability <- function(suitability) {
  
  global_scalar(
    suitability,
    "mean",
    zero_if_missing = FALSE
  )
}


max_suitability <- function(suitability) {
  
  global_scalar(
    suitability,
    "max",
    zero_if_missing = FALSE
  )
}


max_stack_na <- function(x) {
  
  if (all(is.na(x))) {
    return(NA_real_)
  }
  
  max(
    x,
    na.rm = TRUE
  )
}


rank_layer_names <- function(
    number_of_layers,
    suffix) {
  
  paste0(
    "rank",
    seq_len(number_of_layers),
    "_",
    suffix
  )
}


rank_dual_cell <- function(
    x,
    population_zones,
    rank_min,
    suitable_threshold) {
  
  number_of_populations <- length(
    population_zones
  )
  
  output_zone <- rep(
    NA_real_,
    number_of_populations
  )
  
  output_suitability <- rep(
    NA_real_,
    number_of_populations
  )
  
  valid_rank <- which(
    !is.na(x) &
      is.finite(x) &
      x > rank_min
  )
  
  valid_suitable <- which(
    !is.na(x) &
      is.finite(x) &
      x >= suitable_threshold
  )
  
  if (length(valid_rank) > 0) {
    
    # Higher suitability first. Exact ties use smaller source-zone ID.
    order_index <- valid_rank[
      order(
        -x[valid_rank],
        population_zones[valid_rank]
      )
    ]
    
    output_zone[
      seq_along(order_index)
    ] <- population_zones[
      order_index
    ]
    
    output_suitability[
      seq_along(order_index)
    ] <- x[
      order_index
    ]
  }
  
  top1_minus_top2 <- NA_real_
  
  if (length(valid_rank) >= 2) {
    top1_minus_top2 <-
      output_suitability[1] -
      output_suitability[2]
  }
  
  c(
    output_zone,
    output_suitability,
    n_pop_ranked = length(valid_rank),
    n_pop_suitable = length(valid_suitable),
    top1_minus_top2 = top1_minus_top2
  )
}


read_dual_stack <- function(
    files,
    layer_names,
    template) {
  
  stack <- rast(files)
  names(stack) <- layer_names
  
  if (!compareGeom(
    stack,
    template,
    stopOnError = FALSE
  )) {
    
    cat(
      "[RESAMPLE] dual rasters to reference geometry\n"
    )
    
    stack <- resample(
      stack,
      template,
      method = "bilinear"
    )
    
    names(stack) <- layer_names
  }
  
  stack <- mask(
    stack,
    template
  )
  
  names(stack) <- layer_names
  
  stack
}


# 2. Population lookup ===========================================================

if (!file.exists(population_lookup_file)) {
  stop(
    "Missing corrected population lookup:\n",
    population_lookup_file,
    "\nRun script 6.1 first."
  )
}

if (!file.exists(reference_file)) {
  stop(
    "Missing reference raster: ",
    reference_file
  )
}

population <- fread(
  population_lookup_file
)

required_columns <- c(
  "PopulationID",
  "Species",
  "source_zone"
)

if (!all(
  required_columns %in%
  names(population)
)) {
  stop(
    "Population lookup must contain: ",
    paste(
      required_columns,
      collapse = ", "
    )
  )
}

population[
  ,
  `:=`(
    PopulationID = as.character(
      PopulationID
    ),
    Species = as.character(
      Species
    ),
    source_zone = as.integer(
      source_zone
    )
  )
]

if (!"projected" %in% names(population)) {
  population[
    ,
    projected :=
      source_zone %in% model_zoneID
  ]
} else {
  population[
    ,
    projected := as.logical(projected)
  ]
}

if (!"reference_abundance" %in% names(population)) {
  
  if ("Population" %in% names(population)) {
    
    population[
      ,
      reference_abundance :=
        as.numeric(Population)
    ]
    
  } else {
    
    population[
      ,
      reference_abundance :=
        NA_real_
    ]
  }
}

if (!"zone_name" %in% names(population)) {
  population[
    ,
    zone_name := NA_character_
  ]
}

if (!"COLOR" %in% names(population)) {
  population[
    ,
    COLOR := NA_character_
  ]
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

projected_population <- population[
  projected == TRUE &
    source_zone %in% model_zoneID
]

if (!nrow(projected_population)) {
  stop(
    "No projected populations remain."
  )
}

required_zones <- sort(
  unique(
    projected_population$source_zone
  )
)

species_names <- sort(
  unique(
    projected_population$Species
  )
)

fwrite(
  projected_population,
  file.path(
    table_dir,
    "dual_population_projection_lookup_var.csv"
  )
)

fwrite(
  population[
    projected == FALSE |
      !(source_zone %in% model_zoneID)
  ],
  file.path(
    table_dir,
    "dual_population_projection_exclusions_var.csv"
  )
)


# 3. Available dual-suitability jobs ============================================

job_list <- list()

for (method in method_order) {
  for (scenario in future_order) {
    
    files <- dual_file(
      method,
      scenario,
      required_zones
    )
    
    missing_zones <- required_zones[
      !file.exists(files)
    ]
    
    job_list[[length(job_list) + 1L]] <- data.table(
      method = method,
      scenario = scenario,
      n_required_zones =
        length(required_zones),
      n_missing_zones =
        length(missing_zones),
      missing_zones = paste(
        missing_zones,
        collapse = ","
      )
    )
  }
}

jobs <- rbindlist(
  job_list
)

available_jobs <- jobs[
  n_missing_zones == 0
]

missing_jobs <- jobs[
  n_missing_zones > 0
]

fwrite(
  available_jobs,
  file.path(
    table_dir,
    "dual_suitability_jobs_available_var.csv"
  )
)

fwrite(
  missing_jobs,
  file.path(
    table_dir,
    "dual_suitability_jobs_missing_var.csv"
  )
)

if (!nrow(available_jobs)) {
  stop(
    "No complete rf_var or mf_var dual-suitability jobs were found."
  )
}

cat(
  "\n[AVAILABLE JOBS]\n",
  "Available: ",
  nrow(available_jobs),
  " of ",
  nrow(jobs),
  "\n",
  sep = ""
)

print(
  available_jobs[
    ,
    .(
      method,
      scenario
    )
  ]
)


# 4. Template and result containers =============================================

reference_map <- rast(
  reference_file
)

names(reference_map) <- "zoneID"

cell_area <- cellSize(
  reference_map,
  unit = "km"
)

population_area_results <- list()
species_area_results <- list()
population_layer_index <- list()
species_raster_index <- list()
population_rank_index <- list()


# 5. Generate dual-based niches =================================================

for (job_index in seq_len(
  nrow(available_jobs)
)) {
  
  method <- available_jobs$method[
    job_index
  ]
  
  scenario <- available_jobs$scenario[
    job_index
  ]
  
  cat(
    "\n==============================\n",
    "DUAL-BASED TREE NICHE: ",
    method,
    " | ",
    scenario,
    "\n",
    "==============================\n",
    sep = ""
  )
  
  population_dir <- file.path(
    output_root,
    method,
    scenario,
    "population niche"
  )
  
  species_dir <- file.path(
    output_root,
    method,
    scenario,
    "species niche"
  )
  
  species_binary_dir <- file.path(
    output_root,
    method,
    scenario,
    "species niche binary"
  )
  
  rank_dir <- file.path(
    output_root,
    method,
    scenario,
    "population ranking"
  )
  
  dir.create(
    population_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    species_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    species_binary_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    rank_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  for (species_name in species_names) {
    
    species_population <- projected_population[
      Species == species_name
    ][order(source_zone)]
    
    population_zones <-
      species_population$source_zone
    
    if (!length(population_zones)) {
      next
    }
    
    safe_species <- safe_name(
      species_name
    )
    
    population_files <- dual_file(
      method,
      scenario,
      population_zones
    )
    
    if (!all(
      file.exists(population_files)
    )) {
      
      cat(
        "[SKIP SPECIES] missing dual raster:",
        species_name,
        "\n"
      )
      
      next
    }
    
    layer_names <- make.unique(
      paste0(
        "pop_zone",
        population_zones
      )
    )
    
    population_stack_file <- file.path(
      population_dir,
      paste0(
        safe_species,
        "_population_dual_suitability.tif"
      )
    )
    
    species_file <- file.path(
      species_dir,
      paste0(
        safe_species,
        "_species_dual_suitability.tif"
      )
    )
    
    species_binary_file <- file.path(
      species_binary_dir,
      paste0(
        safe_species,
        "_species_dual_binary_threshold",
        suitability_threshold,
        ".tif"
      )
    )
    
    rank_zone_file <- file.path(
      rank_dir,
      paste0(
        safe_species,
        "_ranked_population_zones.tif"
      )
    )
    
    rank_suitability_file <- file.path(
      rank_dir,
      paste0(
        safe_species,
        "_ranked_population_suitability.tif"
      )
    )
    
    rank_summary_file <- file.path(
      rank_dir,
      paste0(
        safe_species,
        "_ranked_population_summary.tif"
      )
    )
    
    rank_index_file <- file.path(
      rank_dir,
      paste0(
        safe_species,
        "_population_rank_index.csv"
      )
    )
    
    population_stack <- read_dual_stack(
      population_files,
      layer_names,
      reference_map
    )
    
    if (write_population_stack) {
      
      if (
        !reuse_existing_outputs ||
        !raster_valid(
          population_stack_file,
          reference_map,
          length(population_zones)
        )
      ) {
        
        remove_raster_files(
          population_stack_file
        )
        
        writeRaster(
          population_stack,
          population_stack_file,
          overwrite = TRUE,
          wopt = list(
            datatype = "FLT4S",
            gdal = "COMPRESS=LZW"
          )
        )
      }
    }
    
    
    # 5.1 Species continuous dual suitability -----------------------------------
    
    species_map <- NULL
    
    if (
      reuse_existing_outputs &&
      raster_valid(
        species_file,
        reference_map,
        1L
      )
    ) {
      
      species_map <- rast(
        species_file
      )
      
      names(species_map) <-
        "species_dual_suitability"
      
      cat(
        "[REUSE SPECIES]",
        species_name,
        "|",
        method,
        "|",
        scenario,
        "\n"
      )
    }
    
    if (is.null(species_map)) {
      
      if (nlyr(population_stack) == 1L) {
        
        species_map <-
          population_stack[[1]]
        
      } else {
        
        remove_raster_files(
          species_file
        )
        
        species_map <- app(
          population_stack,
          fun = max_stack_na,
          filename = species_file,
          overwrite = TRUE,
          wopt = list(
            datatype = "FLT4S",
            gdal = "COMPRESS=LZW"
          )
        )
      }
      
      names(species_map) <-
        "species_dual_suitability"
      
      if (nlyr(population_stack) == 1L) {
        
        remove_raster_files(
          species_file
        )
        
        writeRaster(
          species_map,
          species_file,
          overwrite = TRUE,
          wopt = list(
            datatype = "FLT4S",
            gdal = "COMPRESS=LZW"
          )
        )
      }
    }
    
    
    # 5.2 Binary species niche at threshold 0.4 ---------------------------------
    
    if (
      !reuse_existing_outputs ||
      !raster_valid(
        species_binary_file,
        reference_map,
        1L
      )
    ) {
      
      species_binary <- ifel(
        !is.na(species_map) &
          species_map >=
          suitability_threshold,
        1,
        NA
      )
      
      names(species_binary) <-
        "species_suitable"
      
      remove_raster_files(
        species_binary_file
      )
      
      writeRaster(
        species_binary,
        species_binary_file,
        overwrite = TRUE,
        wopt = list(
          datatype = "INT1U",
          gdal = "COMPRESS=LZW"
        )
      )
      
      rm(species_binary)
      gc()
    }
    
    
    # 5.3 Population ranking within species -------------------------------------
    
    number_of_populations <-
      length(population_zones)
    
    reuse_rank <- (
      reuse_existing_outputs &&
        raster_valid(
          rank_zone_file,
          reference_map,
          number_of_populations
        ) &&
        raster_valid(
          rank_suitability_file,
          reference_map,
          number_of_populations
        ) &&
        raster_valid(
          rank_summary_file,
          reference_map,
          3L
        )
    )
    
    if (!reuse_rank) {
      
      rank_all <- app(
        population_stack,
        fun = function(x) {
          rank_dual_cell(
            x = x,
            population_zones =
              population_zones,
            rank_min =
              rank_min_suitability,
            suitable_threshold =
              suitability_threshold
          )
        }
      )
      
      zone_layers <- seq_len(
        number_of_populations
      )
      
      suitability_layers <-
        number_of_populations +
        seq_len(number_of_populations)
      
      summary_layers <-
        2L * number_of_populations +
        c(1L, 2L, 3L)
      
      rank_zone <- rank_all[[zone_layers]]
      
      names(rank_zone) <- rank_layer_names(
        number_of_populations,
        "zone"
      )
      
      rank_suitability <- rank_all[[suitability_layers]]
      
      names(rank_suitability) <-
        rank_layer_names(
          number_of_populations,
          "suit"
        )
      
      rank_summary <- rank_all[[summary_layers]]
      
      names(rank_summary) <- c(
        "n_pop_ranked",
        "n_pop_suitable",
        "top1_minus_top2"
      )
      
      remove_raster_files(
        rank_zone_file
      )
      
      remove_raster_files(
        rank_suitability_file
      )
      
      remove_raster_files(
        rank_summary_file
      )
      
      writeRaster(
        rank_zone,
        rank_zone_file,
        overwrite = TRUE,
        wopt = list(
          datatype = "INT2S",
          gdal = "COMPRESS=LZW"
        )
      )
      
      writeRaster(
        rank_suitability,
        rank_suitability_file,
        overwrite = TRUE,
        wopt = list(
          datatype = "FLT4S",
          gdal = "COMPRESS=LZW"
        )
      )
      
      writeRaster(
        rank_summary,
        rank_summary_file,
        overwrite = TRUE,
        wopt = list(
          datatype = "FLT4S",
          gdal = "COMPRESS=LZW"
        )
      )
      
      rm(
        rank_all,
        rank_zone,
        rank_suitability,
        rank_summary
      )
      
      gc()
    }
    
    rank_index <- data.table(
      Species = species_name,
      PopulationID =
        species_population$PopulationID,
      source_zone =
        species_population$source_zone,
      zone_name =
        species_population$zone_name,
      method = method,
      scenario = scenario,
      original_layer_in_population_stack =
        seq_along(population_zones),
      rank_min_suitability =
        rank_min_suitability,
      dual_threshold =
        suitability_threshold,
      source_dual_raster =
        population_files,
      ranked_zone_raster =
        rank_zone_file,
      ranked_suitability_raster =
        rank_suitability_file,
      ranked_summary_raster =
        rank_summary_file
    )
    
    fwrite(
      rank_index,
      rank_index_file
    )
    
    population_rank_index[[length(population_rank_index) + 1L]] <- rank_index
    
    
    # 5.4 Population-level summaries --------------------------------------------
    
    population_area <- rbindlist(
      lapply(
        seq_along(population_zones),
        function(layer_index) {
          
          layer <-
            population_stack[[layer_index]]
          
          data.table(
            Species = species_name,
            PopulationID =
              species_population$PopulationID[
                layer_index
              ],
            source_zone =
              species_population$source_zone[
                layer_index
              ],
            zone_name =
              species_population$zone_name[
                layer_index
              ],
            reference_abundance =
              species_population$reference_abundance[
                layer_index
              ],
            method = method,
            scenario = scenario,
            dual_threshold =
              suitability_threshold,
            rank_min_suitability =
              rank_min_suitability,
            suitable_area_km2 =
              suitable_area(
                layer,
                cell_area,
                suitability_threshold
              ),
            suitability_weighted_area_km2 =
              weighted_area(
                layer,
                cell_area
              ),
            mean_dual_suitability =
              mean_suitability(layer),
            max_dual_suitability =
              max_suitability(layer),
            population_niche_definition =
              "source-zone dual suitability raster",
            source_dual_raster =
              population_files[layer_index],
            population_stack_file =
              if (write_population_stack) {
                population_stack_file
              } else {
                NA_character_
              },
            layer = layer_index,
            layer_name =
              layer_names[layer_index],
            COLOR =
              species_population$COLOR[
                layer_index
              ]
          )
        }
      ),
      fill = TRUE
    )
    
    population_area_results[[length(population_area_results) + 1L]] <- population_area
    
    
    # 5.5 Species-level summary --------------------------------------------------
    
    species_area_results[[length(species_area_results) + 1L]] <- data.table(
      Species = species_name,
      method = method,
      scenario = scenario,
      populations_projected =
        length(population_zones),
      dual_threshold =
        suitability_threshold,
      rank_min_suitability =
        rank_min_suitability,
      suitable_area_km2 =
        suitable_area(
          species_map,
          cell_area,
          suitability_threshold
        ),
      suitability_weighted_area_km2 =
        weighted_area(
          species_map,
          cell_area
        ),
      mean_dual_suitability =
        mean_suitability(species_map),
      max_dual_suitability =
        max_suitability(species_map),
      species_niche_definition =
        "pixel-wise maximum of source population dual suitability rasters",
      species_raster = species_file,
      species_binary_raster =
        species_binary_file,
      ranked_zone_raster =
        rank_zone_file,
      ranked_suitability_raster =
        rank_suitability_file,
      ranked_summary_raster =
        rank_summary_file
    )
    
    
    # 5.6 Raster indexes ---------------------------------------------------------
    
    population_layer_index[[length(population_layer_index) + 1L]] <- data.table(
      Species = species_name,
      PopulationID =
        species_population$PopulationID,
      source_zone =
        species_population$source_zone,
      zone_name =
        species_population$zone_name,
      method = method,
      scenario = scenario,
      layer =
        seq_along(population_zones),
      layer_name =
        layer_names,
      source_dual_raster =
        population_files,
      population_stack_file =
        if (write_population_stack) {
          population_stack_file
        } else {
          NA_character_
        }
    )
    
    species_raster_index[[length(species_raster_index) + 1L]] <- data.table(
      Species = species_name,
      method = method,
      scenario = scenario,
      populations_projected =
        length(population_zones),
      species_raster =
        species_file,
      species_binary_raster =
        species_binary_file,
      ranked_zone_raster =
        rank_zone_file,
      ranked_suitability_raster =
        rank_suitability_file,
      ranked_summary_raster =
        rank_summary_file
    )
    
    latest_species_result <-
      species_area_results[[length(species_area_results)]]
    
    cat(
      "[SAVED]",
      species_name,
      "| populations:",
      length(population_zones),
      "| suitable area km2:",
      round(
        latest_species_result$suitable_area_km2,
        2
      ),
      "\n"
    )
    
    rm(
      population_stack,
      species_map,
      population_area,
      rank_index
    )
    
    gc()
  }
}


# 6. Save summary tables =========================================================

if (!length(population_area_results)) {
  stop(
    "No population-level results were generated."
  )
}

population_area <- rbindlist(
  population_area_results,
  fill = TRUE
)

species_area <- rbindlist(
  species_area_results,
  fill = TRUE
)

population_layers <- rbindlist(
  population_layer_index,
  fill = TRUE
)

species_index <- rbindlist(
  species_raster_index,
  fill = TRUE
)

population_rank_index_all <- rbindlist(
  population_rank_index,
  fill = TRUE
)

setorder(
  population_area,
  Species,
  source_zone,
  method,
  scenario
)

setorder(
  species_area,
  Species,
  method,
  scenario
)

setorder(
  population_layers,
  Species,
  source_zone,
  method,
  scenario
)

setorder(
  species_index,
  Species,
  method,
  scenario
)

setorder(
  population_rank_index_all,
  Species,
  source_zone,
  method,
  scenario
)

fwrite(
  population_area,
  file.path(
    table_dir,
    "dual_population_niche_area_var.csv"
  )
)

fwrite(
  species_area,
  file.path(
    table_dir,
    "dual_species_niche_area_var.csv"
  )
)

fwrite(
  population_layers,
  file.path(
    table_dir,
    "dual_population_layer_index_var.csv"
  )
)

fwrite(
  species_index,
  file.path(
    table_dir,
    "dual_species_raster_index_var.csv"
  )
)

fwrite(
  population_rank_index_all,
  file.path(
    table_dir,
    "dual_population_rank_index_var.csv"
  )
)

cat(
  "\nCOMPLETE\n",
  "Methods: ",
  paste(
    method_order,
    collapse = ", "
  ),
  "\n",
  "Dual suitable threshold: ",
  suitability_threshold,
  "\n",
  "Output root: ",
  output_root,
  "\n",
  "Available jobs processed: ",
  nrow(available_jobs),
  " of ",
  nrow(jobs),
  "\n",
  "Population niche: source-zone dual suitability raster\n",
  "Species niche: pixel-wise maximum of population dual suitability rasters\n",
  sep = ""
)
