# 8.1 Rank dual suitability by pixel
# Corrected climate models: rf_var and mf_var
# ==============================================================================
#
# For each method x scenario:
#   1. Read all 53 zone-level dual suitability rasters.
#   2. Rank zones by dual suitability for every pixel.
#   3. Save top-10 zones, top-10 suitability values and uncertainty summaries.
#
# Threshold rule:
#   dual suitability >= 0.4 is above the assignment threshold;
#   novel_by_threshold = 1 only when every zone is < 0.4.
#
# Ranking itself includes every positive dual suitability value.
#
# Outputs:
#   dual suit ranking var/{method}/{scenario}/
#
# Existing script 8 outputs are not overwritten.
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

output_root <- file.path(
  base_dir,
  "dual suit ranking var"
)

table_dir <- file.path(
  output_root,
  "tables"
)

tmp_dir <- file.path(
  base_dir,
  "tmp_dual_suit_ranking_var"
)

reference_file <- file.path(
  base_dir,
  "raster/ecosys_ori.tif"
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

zoneID <- c(
  1:7,
  9:50,
  52:55
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

scenario_order <- c(
  "normal",
  future_order
)

methods_to_run <- method_order
scenarios_to_run <- scenario_order

dual_threshold <- 0.4
rank_min_suitability <- 0
top_k_rank <- 10L

rank_cores <- min(
  4L,
  max(
    1L,
    parallel::detectCores() - 1L
  )
)

reuse_existing_outputs <- TRUE
overwrite_outputs <- TRUE


# 1. Helpers ====================================================================

cat0 <- function(...) {
  cat(
    ...,
    "\n",
    sep = ""
  )
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


regex_escape <- function(x) {
  
  gsub(
    "([][{}()+*^$|\\\\?.])",
    "\\\\\\1",
    x,
    perl = TRUE
  )
}


remove_raster_files <- function(
    filepath,
    must_remove = FALSE) {
  
  base <- tools::file_path_sans_ext(
    filepath
  )
  
  folder <- dirname(filepath)
  base_name <- basename(base)
  
  side_files <- character()
  
  if (dir.exists(folder)) {
    
    side_files <- list.files(
      folder,
      pattern = paste0(
        "^",
        regex_escape(base_name),
        "\\."
      ),
      full.names = TRUE
    )
  }
  
  files <- unique(c(
    filepath,
    paste0(filepath, ".aux.xml"),
    paste0(filepath, ".ovr"),
    paste0(filepath, ".msk"),
    paste0(base, ".aux.xml"),
    paste0(base, ".ovr"),
    paste0(base, ".tfw"),
    side_files
  ))
  
  files <- files[
    file.exists(files)
  ]
  
  if (length(files) > 0) {
    
    gc()
    
    suppressWarnings(
      unlink(
        files,
        force = TRUE,
        recursive = TRUE
      )
    )
    
    Sys.sleep(0.2)
    gc()
  }
  
  remaining <- files[
    file.exists(files)
  ]
  
  if (
    must_remove &&
    length(remaining) > 0
  ) {
    stop(
      "Cannot remove existing raster files:\n",
      paste(
        remaining,
        collapse = "\n"
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


rank_output_files <- function(
    method,
    scenario) {
  
  output_dir <- file.path(
    output_root,
    method,
    scenario
  )
  
  list(
    output_dir = output_dir,
    ranked_zone_file = file.path(
      output_dir,
      "ranked_zone.tif"
    ),
    ranked_suitability_file = file.path(
      output_dir,
      "ranked_suitability.tif"
    ),
    ranked_summary_file = file.path(
      output_dir,
      "ranked_summary.tif"
    )
  )
}


valid_output <- function(
    raster_file,
    expected_layers,
    template) {
  
  if (!file.exists(raster_file)) {
    return(FALSE)
  }
  
  raster <- tryCatch(
    rast(raster_file),
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(raster)) {
    return(FALSE)
  }
  
  if (nlyr(raster) != expected_layers) {
    return(FALSE)
  }
  
  isTRUE(
    compareGeom(
      raster,
      template,
      stopOnError = FALSE
    )
  )
}


valid_rank_outputs <- function(
    method,
    scenario,
    template,
    number_of_ranks,
    number_of_summary_layers) {
  
  output <- rank_output_files(
    method,
    scenario
  )
  
  valid_zone <- valid_output(
    output$ranked_zone_file,
    number_of_ranks,
    template
  )
  
  valid_suitability <- valid_output(
    output$ranked_suitability_file,
    number_of_ranks,
    template
  )
  
  valid_summary <- valid_output(
    output$ranked_summary_file,
    number_of_summary_layers,
    template
  )
  
  list(
    all_valid =
      valid_zone &&
      valid_suitability &&
      valid_summary,
    valid_zone = valid_zone,
    valid_suitability =
      valid_suitability,
    valid_summary =
      valid_summary,
    ranked_zone_file =
      output$ranked_zone_file,
    ranked_suitability_file =
      output$ranked_suitability_file,
    ranked_summary_file =
      output$ranked_summary_file
  )
}


make_layer_index <- function(
    method,
    scenario,
    number_of_ranks,
    summary_layer_names) {
  
  rbindlist(
    list(
      data.table(
        method = method,
        scenario = scenario,
        raster_type = "ranked_zone",
        layer = seq_len(
          number_of_ranks
        ),
        layer_name = rank_layer_names(
          number_of_ranks,
          "zone"
        ),
        rank = seq_len(
          number_of_ranks
        ),
        meaning =
          "zone ID at this rank"
      ),
      data.table(
        method = method,
        scenario = scenario,
        raster_type =
          "ranked_suitability",
        layer = seq_len(
          number_of_ranks
        ),
        layer_name = rank_layer_names(
          number_of_ranks,
          "suit"
        ),
        rank = seq_len(
          number_of_ranks
        ),
        meaning =
          "dual suitability at this rank"
      ),
      data.table(
        method = method,
        scenario = scenario,
        raster_type =
          "ranked_summary",
        layer = seq_along(
          summary_layer_names
        ),
        layer_name =
          summary_layer_names,
        rank = NA_integer_,
        meaning = c(
          "number of zones with positive dual suitability",
          "number of zones with dual suitability >= 0.4",
          "rank1 suitability minus rank2 suitability",
          "highest positive dual suitability",
          "second-highest positive dual suitability",
          "1 when every zone is < 0.4, otherwise 0"
        )
      )
    )
  )
}


read_dual_stack <- function(
    files,
    template,
    valid_mask) {
  
  stack <- rast(files)
  names(stack) <- as.character(zoneID)
  
  if (!compareGeom(
    stack,
    template,
    stopOnError = FALSE
  )) {
    
    cat0(
      "[RESAMPLE] dual suitability rasters to reference geometry"
    )
    
    stack <- resample(
      stack,
      template,
      method = "bilinear"
    )
    
    names(stack) <-
      as.character(zoneID)
  }
  
  stack <- mask(
    stack,
    valid_mask
  )
  
  names(stack) <-
    as.character(zoneID)
  
  stack
}


make_rank_function <- function(
    zone_ids,
    rank_min,
    threshold,
    top_k) {
  
  force(zone_ids)
  force(rank_min)
  force(threshold)
  force(top_k)
  
  function(x) {
    
    number_of_ranks <- min(
      top_k,
      length(zone_ids)
    )
    
    output_zone <- rep(
      NA_real_,
      number_of_ranks
    )
    
    output_suitability <- rep(
      NA_real_,
      number_of_ranks
    )
    
    valid_rank <- which(
      !is.na(x) &
        is.finite(x) &
        x > rank_min
    )
    
    valid_above_threshold <- which(
      !is.na(x) &
        is.finite(x) &
        x >= threshold
    )
    
    if (length(valid_rank) > 0) {
      
      order_index <- valid_rank[
        order(
          -x[valid_rank],
          zone_ids[valid_rank]
        )
      ]
      
      order_index <- order_index[
        seq_len(
          min(
            length(order_index),
            number_of_ranks
          )
        )
      ]
      
      output_zone[
        seq_along(order_index)
      ] <- zone_ids[
        order_index
      ]
      
      output_suitability[
        seq_along(order_index)
      ] <- x[
        order_index
      ]
    }
    
    top1_suitability <- if (
      !is.na(output_suitability[1])
    ) {
      output_suitability[1]
    } else {
      NA_real_
    }
    
    top2_suitability <- if (
      number_of_ranks >= 2L &&
      !is.na(output_suitability[2])
    ) {
      output_suitability[2]
    } else {
      NA_real_
    }
    
    top1_minus_top2 <- if (
      is.finite(top1_suitability) &&
      is.finite(top2_suitability)
    ) {
      top1_suitability -
        top2_suitability
    } else {
      NA_real_
    }
    
    novel_by_threshold <- as.numeric(
      length(valid_above_threshold) == 0
    )
    
    c(
      output_zone,
      output_suitability,
      n_zone_ranked =
        length(valid_rank),
      n_zone_above_threshold =
        length(valid_above_threshold),
      top1_minus_top2 =
        top1_minus_top2,
      top1_suit =
        top1_suitability,
      top2_suit =
        top2_suitability,
      novel_by_threshold =
        novel_by_threshold
    )
  }
}


# 2. Reference map ==============================================================

if (!file.exists(reference_file)) {
  stop(
    "Missing reference raster: ",
    reference_file
  )
}

reference_map <- rast(
  reference_file
)

names(reference_map) <-
  "original_zone"

valid_mask <- ifel(
  !is.na(reference_map) &
    reference_map != 8,
  1,
  NA
)

number_of_ranks <- min(
  top_k_rank,
  length(zoneID)
)

number_of_summary_layers <- 6L

summary_layer_names <- c(
  "n_zone_ranked",
  "n_zone_above_threshold",
  "top1_minus_top2",
  "top1_suit",
  "top2_suit",
  "novel_by_threshold"
)

cat(
  "\n[RANKING SETTINGS]\n",
  "Methods: ",
  paste(
    methods_to_run,
    collapse = ", "
  ),
  "\n",
  "Top K: ",
  number_of_ranks,
  "\n",
  "Rank minimum: ",
  rank_min_suitability,
  "\n",
  "Dual threshold: ",
  dual_threshold,
  "\n",
  "Cores: ",
  rank_cores,
  "\n",
  sep = ""
)


# 3. Build jobs =================================================================

job_list <- list()

for (method in methods_to_run) {
  for (scenario in scenarios_to_run) {
    
    files <- dual_file(
      method,
      scenario,
      zoneID
    )
    
    missing_zones <- zoneID[
      !file.exists(files)
    ]
    
    output_validity <- valid_rank_outputs(
      method,
      scenario,
      reference_map,
      number_of_ranks,
      number_of_summary_layers
    )
    
    job_list[[length(job_list) + 1L]] <- data.table(
      method = method,
      scenario = scenario,
      n_required_zones =
        length(zoneID),
      n_missing_zones =
        length(missing_zones),
      missing_zones = paste(
        missing_zones,
        collapse = ","
      ),
      input_complete =
        length(missing_zones) == 0,
      output_valid =
        output_validity$all_valid,
      valid_ranked_zone =
        output_validity$valid_zone,
      valid_ranked_suitability =
        output_validity$valid_suitability,
      valid_ranked_summary =
        output_validity$valid_summary,
      ranked_zone_file =
        output_validity$ranked_zone_file,
      ranked_suitability_file =
        output_validity$ranked_suitability_file,
      ranked_summary_file =
        output_validity$ranked_summary_file
    )
  }
}

jobs <- rbindlist(
  job_list
)

setorder(
  jobs,
  method,
  scenario
)

fwrite(
  jobs,
  file.path(
    table_dir,
    "ranking_jobs_all_var.csv"
  )
)

fwrite(
  jobs[
    input_complete == TRUE
  ],
  file.path(
    table_dir,
    "ranking_jobs_input_complete_var.csv"
  )
)

fwrite(
  jobs[
    input_complete == FALSE &
      output_valid == FALSE
  ],
  file.path(
    table_dir,
    "ranking_jobs_missing_input_and_output_var.csv"
  )
)


# 4. Rank dual suitability =======================================================

output_index <- list()
layer_index <- list()

for (job_index in seq_len(
  nrow(jobs)
)) {
  
  method <- jobs$method[
    job_index
  ]
  
  scenario <- jobs$scenario[
    job_index
  ]
  
  cat(
    "\n==============================\n",
    "RANK DUAL SUITABILITY: ",
    method,
    " | ",
    scenario,
    "\n",
    "==============================\n",
    sep = ""
  )
  
  output <- rank_output_files(
    method,
    scenario
  )
  
  dir.create(
    output$output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (
    reuse_existing_outputs &&
    isTRUE(
      jobs$output_valid[
        job_index
      ]
    )
  ) {
    
    cat0(
      "[REUSE] ",
      method,
      " | ",
      scenario
    )
    
    output_index[[length(output_index) + 1L]] <- data.table(
      method = method,
      scenario = scenario,
      status = "reused",
      message = NA_character_,
      n_input_zones =
        length(zoneID),
      n_missing_zones =
        jobs$n_missing_zones[
          job_index
        ],
      missing_zones =
        jobs$missing_zones[
          job_index
        ],
      top_k_rank =
        number_of_ranks,
      rank_min_suitability =
        rank_min_suitability,
      dual_threshold =
        dual_threshold,
      ranked_zone_file =
        output$ranked_zone_file,
      ranked_suitability_file =
        output$ranked_suitability_file,
      ranked_summary_file =
        output$ranked_summary_file
    )
    
    layer_index[[length(layer_index) + 1L]] <- make_layer_index(
      method,
      scenario,
      number_of_ranks,
      summary_layer_names
    )
    
    next
  }
  
  if (!isTRUE(
    jobs$input_complete[
      job_index
    ]
  )) {
    
    message <- paste0(
      "Missing dual suitability rasters for zones: ",
      jobs$missing_zones[
        job_index
      ]
    )
    
    cat0(
      "[SKIP] ",
      message
    )
    
    output_index[[length(output_index) + 1L]] <- data.table(
      method = method,
      scenario = scenario,
      status =
        "skipped_missing_input",
      message = message,
      n_input_zones =
        length(zoneID) -
        jobs$n_missing_zones[
          job_index
        ],
      n_missing_zones =
        jobs$n_missing_zones[
          job_index
        ],
      missing_zones =
        jobs$missing_zones[
          job_index
        ],
      top_k_rank =
        number_of_ranks,
      rank_min_suitability =
        rank_min_suitability,
      dual_threshold =
        dual_threshold,
      ranked_zone_file =
        output$ranked_zone_file,
      ranked_suitability_file =
        output$ranked_suitability_file,
      ranked_summary_file =
        output$ranked_summary_file
    )
    
    next
  }
  
  tryCatch({
    
    files <- dual_file(
      method,
      scenario,
      zoneID
    )
    
    names(files) <- zoneID
    
    temporary_rank_file <- tempfile(
      pattern = paste0(
        "tmp_rank_",
        safe_name(method),
        "_",
        safe_name(scenario),
        "_"
      ),
      tmpdir = tmp_dir,
      fileext = ".tif"
    )
    
    dual_stack <- read_dual_stack(
      files,
      reference_map,
      valid_mask
    )
    
    if (overwrite_outputs) {
      
      remove_raster_files(
        output$ranked_zone_file,
        must_remove = TRUE
      )
      
      remove_raster_files(
        output$ranked_suitability_file,
        must_remove = TRUE
      )
      
      remove_raster_files(
        output$ranked_summary_file,
        must_remove = TRUE
      )
    }
    
    remove_raster_files(
      temporary_rank_file,
      must_remove = TRUE
    )
    
    rank_function <- make_rank_function(
      zoneID,
      rank_min_suitability,
      dual_threshold,
      number_of_ranks
    )
    
    rank_all <- app(
      dual_stack,
      fun = rank_function,
      filename = temporary_rank_file,
      overwrite = TRUE,
      cores = rank_cores,
      wopt = list(
        datatype = "FLT4S",
        gdal = "COMPRESS=NONE"
      )
    )
    
    zone_layers <- seq_len(
      number_of_ranks
    )
    
    suitability_layers <-
      number_of_ranks +
      seq_len(number_of_ranks)
    
    summary_layers <-
      2L * number_of_ranks +
      seq_len(
        number_of_summary_layers
      )
    
    ranked_zone <- rank_all[[zone_layers]]
    
    names(ranked_zone) <- rank_layer_names(
      number_of_ranks,
      "zone"
    )
    
    ranked_suitability <- rank_all[[suitability_layers]]
    
    names(ranked_suitability) <-
      rank_layer_names(
        number_of_ranks,
        "suit"
      )
    
    ranked_summary <- rank_all[[summary_layers]]
    
    names(ranked_summary) <-
      summary_layer_names
    
    writeRaster(
      ranked_zone,
      output$ranked_zone_file,
      overwrite = TRUE,
      wopt = list(
        datatype = "INT2S",
        gdal = "COMPRESS=LZW"
      )
    )
    
    writeRaster(
      ranked_suitability,
      output$ranked_suitability_file,
      overwrite = TRUE,
      wopt = list(
        datatype = "FLT4S",
        gdal = "COMPRESS=LZW"
      )
    )
    
    writeRaster(
      ranked_summary,
      output$ranked_summary_file,
      overwrite = TRUE,
      wopt = list(
        datatype = "FLT4S",
        gdal = "COMPRESS=LZW"
      )
    )
    
    output_index[[length(output_index) + 1L]] <- data.table(
      method = method,
      scenario = scenario,
      status = "created",
      message = NA_character_,
      n_input_zones =
        length(zoneID),
      n_missing_zones = 0L,
      missing_zones = "",
      top_k_rank =
        number_of_ranks,
      rank_min_suitability =
        rank_min_suitability,
      dual_threshold =
        dual_threshold,
      ranked_zone_file =
        output$ranked_zone_file,
      ranked_suitability_file =
        output$ranked_suitability_file,
      ranked_summary_file =
        output$ranked_summary_file
    )
    
    layer_index[[length(layer_index) + 1L]] <- make_layer_index(
      method,
      scenario,
      number_of_ranks,
      summary_layer_names
    )
    
    cat(
      "[SAVED]\n",
      "  ",
      output$ranked_zone_file,
      "\n",
      "  ",
      output$ranked_suitability_file,
      "\n",
      "  ",
      output$ranked_summary_file,
      "\n",
      sep = ""
    )
    
    rm(
      dual_stack,
      rank_function,
      rank_all,
      ranked_zone,
      ranked_suitability,
      ranked_summary
    )
    
    gc()
    
    remove_raster_files(
      temporary_rank_file
    )
    
  }, error = function(e) {
    
    error_message <- conditionMessage(e)
    
    cat0(
      "[ERROR] ",
      method,
      " | ",
      scenario,
      ": ",
      error_message
    )
    
    output_index[[length(output_index) + 1L]] <<- data.table(
      method = method,
      scenario = scenario,
      status = "error",
      message = error_message,
      n_input_zones =
        length(zoneID),
      n_missing_zones = 0L,
      missing_zones = "",
      top_k_rank =
        number_of_ranks,
      rank_min_suitability =
        rank_min_suitability,
      dual_threshold =
        dual_threshold,
      ranked_zone_file =
        output$ranked_zone_file,
      ranked_suitability_file =
        output$ranked_suitability_file,
      ranked_summary_file =
        output$ranked_summary_file
    )
    
    gc()
  })
}


# 5. Save indexes ================================================================

output_index <- rbindlist(
  output_index,
  fill = TRUE
)

if (length(layer_index) > 0) {
  
  layer_index <- rbindlist(
    layer_index,
    fill = TRUE
  )
  
} else {
  
  layer_index <- data.table(
    method = character(),
    scenario = character(),
    raster_type = character(),
    layer = integer(),
    layer_name = character(),
    rank = integer(),
    meaning = character()
  )
}

setorder(
  output_index,
  method,
  scenario
)

setorder(
  layer_index,
  method,
  scenario,
  raster_type,
  layer
)

fwrite(
  output_index,
  file.path(
    table_dir,
    "ranking_output_index_var.csv"
  )
)

fwrite(
  layer_index,
  file.path(
    table_dir,
    "ranking_layer_index_var.csv"
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
  "Output root: ",
  output_root,
  "\n",
  "Created: ",
  output_index[
    ,
    sum(status == "created")
  ],
  "\n",
  "Reused: ",
  output_index[
    ,
    sum(status == "reused")
  ],
  "\n",
  "Errors: ",
  output_index[
    ,
    sum(status == "error")
  ],
  "\n",
  "Dual threshold: ",
  dual_threshold,
  "\n",
  sep = ""
)
