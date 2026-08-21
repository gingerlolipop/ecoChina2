# 11.7 Manuscript Top-k diagnostics and population-coloured species maps
# ==============================================================================
# Run after scripts 8.2, 11.2, 11.3, and 11.4.
#
# This script does not refit models, repeat projections, or recalculate the
# Top-k area curves. It reads existing assigned maps, ranked-zone rasters,
# ranked-suitability rasters, and visualization tables to create two missing
# reader-facing products:
#
#   1. A three-row Top-1/Top-3/Top-5 zone-level F1 bubble comparison against
#      the Multiclass RF, evaluated on one common reference-period mask.
#   2. Future species-niche maps in which each contributing source-ecotype
#      population is drawn with its reference ecotype colour. Transparent
#      overlays show cells represented by more than one source population.
#
# Top-k definitions
# -----------------
# Reference-period bubble plot:
#   Top-1 is the existing assigned map, including the 1e-4 tie rule. For
#   Top-3 and Top-5, the observed reference ecotype is restored when it occurs
#   within the first k ranked analogues. These are reference-conditioned
#   agreement diagnostics, not independent predictions.
#
# Future species maps:
#   Top-1 uses the existing assigned map. Top-3 and Top-5 include the assigned
#   Top-1 ecotype plus ecotypes occurring within the first k ranked positions
#   when their dual suitability is >= 0.4. A species is mapped wherever at
#   least one of its source ecotypes meets that criterion. This is the same
#   rank-conditioned definition used by Figure_var_8_assigned_species_area_topk;
#   it is distinct from the continuous species envelope formed by the cellwise
#   maximum of all source-population suitability surfaces.
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

visualization_dir <- file.path(
  base_dir,
  "visualization var threshold0.4"
)

figure_dir <- file.path(
  visualization_dir,
  "figures"
)

table_dir <- file.path(
  visualization_dir,
  "tables"
)

chord_dir <- file.path(
  figure_dir,
  "chord diagrams"
)

species_map_dir <- file.path(
  figure_dir,
  "species population overlay maps"
)

dir.create(
  species_map_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

reference_file <- file.path(
  base_dir,
  "raster",
  "ecosys_ori.tif"
)

result_map_root <- file.path(
  base_dir,
  "result maps"
)

ranking_root <- file.path(
  base_dir,
  "dual suit ranking var"
)

multiclass_reference_file <- file.path(
  result_map_root,
  "multiclass_rf",
  "assigned_zone_normal_multiclass_rf.tif"
)

existing_top1_bubble_table <- file.path(
  table_dir,
  "Figure10a_all_binary_workflows_vs_multiclass_reference_map_F1.csv"
)

population_source_file <- file.path(
  table_dir,
  "Figure8_species_population_source_zones.csv"
)

zone_metadata_file <- file.path(
  assessment_dir,
  "feature_importance_analysis",
  "tables",
  "FI_03_zone_metadata_and_reference_area.csv"
)

output_bubble <- file.path(
  figure_dir,
  "Figure_var_10d_reference_map_F1_vs_multiclass_bubble_topk.png"
)

output_bubble_table <- file.path(
  table_dir,
  "Figure_var_10d_reference_map_F1_vs_multiclass_bubble_topk.csv"
)

output_bubble_summary <- file.path(
  table_dir,
  "Figure_var_10d_reference_map_F1_vs_multiclass_bubble_topk_summary.csv"
)

output_species_map_index <- file.path(
  table_dir,
  "Figure_var_9_population_coloured_species_niche_topk_map_index.csv"
)

output_population_colour_key <- file.path(
  species_map_dir,
  "Figure_S_source_population_zone_colour_key.png"
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

rank_cutoffs <- c(
  1L,
  3L,
  5L
)

model_zoneID <- c(
  1:7,
  9:50,
  52:55
)

dual_threshold <- 0.4
tie_tol <- 1e-4
novel_value <- 99L

# The rasters are reduced only for plotting. All areas in Figure 8 remain the
# exact, full-resolution values already produced by script 11.4.
target_plot_cells <- 120000L
population_alpha <- 0.48
chunk_rows <- 48L
reuse_existing_species_maps <- TRUE


# 1. General helpers ============================================================

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path)
  }
  path
}

assigned_map_file <- function(method, scenario) {
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

ranked_zone_file <- function(method, scenario) {
  file.path(
    ranking_root,
    method,
    scenario,
    "ranked_zone.tif"
  )
}

ranked_suitability_file <- function(method, scenario) {
  file.path(
    ranking_root,
    method,
    scenario,
    "ranked_suitability.tif"
  )
}

check_geometry <- function(x, reference, label) {
  if (!compareGeom(
    x,
    reference,
    stopOnError = FALSE
  )) {
    stop("Raster geometry mismatch: ", label)
  }
  invisible(TRUE)
}

matrix_add_pairs <- function(
    target,
    observed,
    predicted,
    keep,
    zone_values) {
  observed_index <- match(
    observed[keep],
    zone_values
  )
  predicted_index <- match(
    predicted[keep],
    zone_values
  )
  
  pair_index <- (
    (observed_index - 1L) *
      length(zone_values) +
      predicted_index
  )
  
  target + matrix(
    tabulate(
      pair_index,
      nbins = length(zone_values)^2
    ),
    nrow = length(zone_values),
    byrow = TRUE
  )
}

confusion_to_metrics <- function(
    confusion_matrix,
    method,
    rank_cutoff) {
  true_positive <- diag(
    confusion_matrix
  )
  
  observed_total <- rowSums(
    confusion_matrix
  )
  
  predicted_total <- colSums(
    confusion_matrix
  )
  
  false_positive <- predicted_total -
    true_positive
  
  false_negative <- observed_total -
    true_positive
  
  denominator <- 2 * true_positive +
    false_positive +
    false_negative
  
  f1 <- rep(
    NA_real_,
    length(model_zoneID)
  )
  
  valid <- denominator > 0
  
  f1[valid] <- 2 * true_positive[valid] /
    denominator[valid]
  
  data.table(
    method = method,
    rank_cutoff = as.integer(rank_cutoff),
    zone = model_zoneID,
    true_positive = as.numeric(true_positive),
    false_positive = as.numeric(false_positive),
    false_negative = as.numeric(false_negative),
    f1 = as.numeric(f1)
  )
}


# 2. Exact common-mask Top-k bubble data ========================================

calculate_reference_topk_metrics <- function() {
  reference_map <- rast(
    require_file(reference_file)
  )[[1]]
  
  multiclass_map <- rast(
    require_file(multiclass_reference_file)
  )[[1]]
  
  assigned_maps <- lapply(
    method_order,
    function(method) {
      rast(
        require_file(
          assigned_map_file(
            method,
            "normal"
          )
        )
      )[[1]]
    }
  )
  
  names(assigned_maps) <- method_order
  
  ranked_maps <- lapply(
    method_order,
    function(method) {
      rast(
        require_file(
          ranked_zone_file(
            method,
            "normal"
          )
        )
      )[[1:5]]
    }
  )
  
  names(ranked_maps) <- method_order
  
  check_geometry(
    multiclass_map,
    reference_map,
    "Multiclass RF reference map"
  )
  
  for (method in method_order) {
    check_geometry(
      assigned_maps[[method]],
      reference_map,
      paste(method, "assigned reference map")
    )
    
    check_geometry(
      ranked_maps[[method]],
      reference_map,
      paste(method, "ranked reference map")
    )
  }
  
  cell_area <- cellSize(
    reference_map,
    unit = "km"
  )
  
  all_read_rasters <- c(
    list(
      reference_map,
      multiclass_map,
      cell_area
    ),
    assigned_maps,
    ranked_maps
  )
  
  invisible(
    lapply(
      all_read_rasters,
      readStart
    )
  )
  
  on.exit(
    invisible(
      lapply(
        all_read_rasters,
        readStop
      )
    ),
    add = TRUE
  )
  
  n_zones <- length(
    model_zoneID
  )
  
  multiclass_confusion <- matrix(
    0,
    nrow = n_zones,
    ncol = n_zones
  )
  
  topk_confusion <- list()
  
  for (method in method_order) {
    for (k in rank_cutoffs) {
      topk_confusion[[
        paste(
          method,
          k,
          sep = "_top"
        )
      ]] <- matrix(
        0,
        nrow = n_zones,
        ncol = n_zones
      )
    }
  }
  
  reference_area_km2 <- numeric(
    n_zones
  )
  
  common_pixels <- 0
  
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
    
    observed <- as.integer(
      readValues(
        reference_map,
        row = row_start,
        nrows = nrows_now,
        mat = FALSE
      )
    )
    
    multiclass <- as.integer(
      readValues(
        multiclass_map,
        row = row_start,
        nrows = nrows_now,
        mat = FALSE
      )
    )
    
    area <- as.numeric(
      readValues(
        cell_area,
        row = row_start,
        nrows = nrows_now,
        mat = FALSE
      )
    )
    
    assigned_values <- lapply(
      method_order,
      function(method) {
        as.integer(
          readValues(
            assigned_maps[[method]],
            row = row_start,
            nrows = nrows_now,
            mat = FALSE
          )
        )
      }
    )
    
    names(assigned_values) <- method_order
    
    ranked_values <- lapply(
      method_order,
      function(method) {
        value <- readValues(
          ranked_maps[[method]],
          row = row_start,
          nrows = nrows_now,
          mat = TRUE
        )
        
        if (is.null(dim(value))) {
          value <- matrix(
            value,
            ncol = 5L
          )
        }
        
        storage.mode(value) <- "integer"
        value
      }
    )
    
    names(ranked_values) <- method_order
    
    common <- (
      observed %in% model_zoneID &
        multiclass %in% model_zoneID &
        is.finite(area)
    )
    
    for (method in method_order) {
      common <- common &
        assigned_values[[method]] %in%
        model_zoneID
    }
    
    if (!any(common)) {
      next
    }
    
    common_pixels <- common_pixels +
      sum(common)
    
    observed_index <- match(
      observed[common],
      model_zoneID
    )
    
    area_by_zone <- rowsum(
      area[common],
      observed_index,
      reorder = FALSE
    )
    
    area_rows <- as.integer(
      rownames(area_by_zone)
    )
    
    reference_area_km2[area_rows] <-
      reference_area_km2[area_rows] +
      as.numeric(
        area_by_zone[, 1]
      )
    
    multiclass_confusion <- matrix_add_pairs(
      multiclass_confusion,
      observed,
      multiclass,
      common,
      model_zoneID
    )
    
    for (method in method_order) {
      assigned <- assigned_values[[method]]
      ranked <- ranked_values[[method]]
      
      for (k in rank_cutoffs) {
        predicted <- assigned
        
        if (k > 1L) {
          restored <- rowSums(
            ranked[
              ,
              seq_len(k),
              drop = FALSE
            ] == observed,
            na.rm = TRUE
          ) > 0
          
          predicted[restored] <-
            observed[restored]
        }
        
        key <- paste(
          method,
          k,
          sep = "_top"
        )
        
        topk_confusion[[key]] <-
          matrix_add_pairs(
            topk_confusion[[key]],
            observed,
            predicted,
            common,
            model_zoneID
          )
      }
    }
  }
  
  if (common_pixels <= 0) {
    stop(
      "No common valid reference-period pixels were found."
    )
  }
  
  multiclass_metrics <- confusion_to_metrics(
    multiclass_confusion,
    "multiclass_rf",
    1L
  )[
    ,
    .(
      zone,
      multiclass_f1 = f1
    )
  ]
  
  topk_metrics <- rbindlist(
    lapply(
      method_order,
      function(method) {
        rbindlist(
          lapply(
            rank_cutoffs,
            function(k) {
              key <- paste(
                method,
                k,
                sep = "_top"
              )
              
              confusion_to_metrics(
                topk_confusion[[key]],
                method,
                k
              )
            }
          )
        )
      }
    )
  )
  
  topk_metrics <- merge(
    topk_metrics,
    multiclass_metrics,
    by = "zone",
    all.x = TRUE,
    sort = FALSE
  )
  
  topk_metrics[
    ,
    `:=`(
      method_label = method_labels[method],
      rank_label = paste0(
        "Top-",
        rank_cutoff
      ),
      common_reference_area_km2 = reference_area_km2[
        match(
          zone,
          model_zoneID
        )
      ],
      common_valid_pixels = common_pixels
    )
  ]
  
  topk_metrics[]
}


topk_bubble <- calculate_reference_topk_metrics()

zone_metadata <- fread(
  require_file(
    zone_metadata_file
  )
)[
  zoneID %in% model_zoneID,
  .(
    zone = as.integer(zoneID),
    zone_name = as.character(zone_name),
    category_label = as.character(category_label),
    zone_color = as.character(zone_color)
  )
]

topk_bubble <- merge(
  topk_bubble,
  zone_metadata,
  by = "zone",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(
  topk_bubble$zone_color
)) {
  stop(
    "At least one modeled zone lacks a plotting colour."
  )
}

# Regression check: the new common-mask Top-1 row must reproduce the existing
# Figure 10d source table exactly, apart from floating-point representation.
existing_top1 <- fread(
  require_file(
    existing_top1_bubble_table
  )
)[
  method %in% method_order,
  .(
    method,
    zone = as.integer(zone),
    old_binary_f1 = as.numeric(binary_f1),
    old_multiclass_f1 = as.numeric(multiclass_f1),
    reference_area_km2 = as.numeric(area_km2)
  )
]

area_lookup <- existing_top1[
  ,
  .(
    reference_area_km2 = mean(
      reference_area_km2,
      na.rm = TRUE
    )
  ),
  by = zone
]

topk_bubble <- merge(
  topk_bubble,
  area_lookup,
  by = "zone",
  all.x = TRUE,
  sort = FALSE
)

top1_check <- merge(
  topk_bubble[
    rank_cutoff == 1L,
    .(
      method,
      zone,
      new_binary_f1 = f1,
      new_multiclass_f1 = multiclass_f1
    )
  ],
  existing_top1,
  by = c(
    "method",
    "zone"
  ),
  all = TRUE
)

if (anyNA(
  top1_check
)) {
  stop(
    "Top-1 bubble regression check contains missing values."
  )
}

if (
  max(
    abs(
      top1_check$new_binary_f1 -
      top1_check$old_binary_f1
    )
  ) > 1e-10 ||
  max(
    abs(
      top1_check$new_multiclass_f1 -
      top1_check$old_multiclass_f1
    )
  ) > 1e-10
) {
  stop(
    "The Top-1 row does not reproduce the existing Figure 10d data."
  )
}

topk_wide <- dcast(
  topk_bubble,
  method + zone ~ rank_cutoff,
  value.var = "f1"
)

if (any(
  topk_wide[["1"]] >
  topk_wide[["3"]] + 1e-12 |
  topk_wide[["3"]] >
  topk_wide[["5"]] + 1e-12
)) {
  stop(
    "Zone-level Top-k F1 is not monotonic."
  )
}

topk_bubble_summary <- topk_bubble[
  ,
  .(
    mean_f1 = mean(
      f1,
      na.rm = TRUE
    ),
    se_f1 = sd(
      f1,
      na.rm = TRUE
    ) /
      sqrt(
        sum(
          is.finite(f1)
        )
      ),
    median_f1 = median(
      f1,
      na.rm = TRUE
    ),
    min_f1 = min(
      f1,
      na.rm = TRUE
    ),
    max_f1 = max(
      f1,
      na.rm = TRUE
    )
  ),
  by = .(
    method,
    method_label,
    rank_cutoff,
    rank_label,
    common_valid_pixels
  )
]

fwrite(
  topk_bubble,
  output_bubble_table
)

fwrite(
  topk_bubble_summary,
  output_bubble_summary
)

topk_bubble[
  ,
  `:=`(
    method_label = factor(
      method_label,
      levels = method_labels[method_order]
    ),
    rank_label = factor(
      rank_label,
      levels = paste0(
        "Top-",
        rank_cutoffs
      )
    ),
    reference_area_million_km2 =
      reference_area_km2 /
      1e6
  )
]

zone_fill <- setNames(
  zone_metadata$zone_color,
  as.character(
    zone_metadata$zone
  )
)

figure_10d_topk <- ggplot(
  topk_bubble,
  aes(
    x = multiclass_f1,
    y = f1,
    size = reference_area_million_km2,
    fill = factor(zone)
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey48"
  ) +
  geom_point(
    shape = 21,
    color = "grey18",
    alpha = 0.72,
    stroke = 0.20
  ) +
  facet_grid(
    rank_label ~ method_label
  ) +
  scale_fill_manual(
    values = zone_fill,
    guide = "none"
  ) +
  scale_size_area(
    max_size = 6.5,
    name = expression(
      "Reference ecotype area (million km"^2*")"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +
  coord_fixed(
    ratio = 1
  ) +
  labs(
    x = "Multiclass RF zone-level F1",
    y = "Individual-zone workflow F1 under Top-k agreement",
    title = "Reference-period ecotype agreement under Top-1, Top-3, and Top-5 criteria",
    subtitle = paste(
      "Top-1 uses the assigned maps; Top-3 and Top-5 restore the observed ecotype",
      "only for reference-conditioned diagnostics. All panels use one common valid mask."
    )
  ) +
  theme_bw(
    base_size = 9.8
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 9.2
    ),
    plot.title = element_text(
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 8.8
    )
  )

ggsave(
  filename = output_bubble,
  plot = figure_10d_topk,
  width = 11.8,
  height = 13.0,
  units = "in",
  dpi = 300,
  bg = "white"
)


# 3. Confirm reused Top-k outputs ==============================================

reused_topk_files <- c(
  file.path(
    figure_dir,
    "Figure_var_8_assigned_species_area_topk.png"
  ),
  unlist(
    lapply(
      method_order,
      function(method) {
        unlist(
          lapply(
            c(
              "category",
              "zone"
            ),
            function(level) {
              file.path(
                chord_dir,
                paste0(
                  "reference_map_",
                  level,
                  "_chord_",
                  method,
                  "_top",
                  rank_cutoffs,
                  ".pdf"
                )
              )
            }
          )
        )
      }
    )
  )
)

missing_reused <- reused_topk_files[
  !file.exists(
    reused_topk_files
  )
]

if (length(
  missing_reused
) > 0L) {
  stop(
    "Previously generated Top-k outputs are missing: ",
    paste(
      missing_reused,
      collapse = "; "
    )
  )
}


# 4. Population-coloured species-niche maps ====================================

population_source <- fread(
  require_file(
    population_source_file
  )
)

required_population_columns <- c(
  "Species",
  "source_zone",
  "zone_name",
  "COLOR"
)

missing_population_columns <- setdiff(
  required_population_columns,
  names(population_source)
)

if (length(
  missing_population_columns
) > 0L) {
  stop(
    "Population-source table is missing columns: ",
    paste(
      missing_population_columns,
      collapse = ", "
    )
  )
}

population_source[
  ,
  `:=`(
    Species = as.character(Species),
    source_zone = as.integer(source_zone),
    zone_name = as.character(zone_name),
    COLOR = as.character(COLOR)
  )
]

if ("projected" %in% names(
  population_source
)) {
  projected_text <- tolower(
    as.character(
      population_source$projected
    )
  )
  
  population_source <- population_source[
    projected_text %in% c(
      "true",
      "t",
      "1",
      "yes"
    )
  ]
}

population_source <- unique(
  population_source[
    source_zone %in% model_zoneID,
    .(
      Species,
      source_zone,
      zone_name,
      COLOR
    )
  ]
)

if (!nrow(
  population_source
)) {
  stop(
    "No projection-eligible source populations were found."
  )
}

species_order <- c(
  "chamFor",
  "cyclLon",
  "lariGme",
  "lariOlg",
  "phylPub",
  "pinuSyl",
  "pinuYun",
  "querVar",
  "robiPse",
  "saliMat"
)

species_order <- species_order[
  species_order %in%
    unique(
      population_source$Species
    )
]

population_palette <- unique(
  population_source[
    order(source_zone),
    .(
      source_zone,
      zone_name,
      COLOR
    )
  ]
)

if (population_palette[
  ,
  any(
    uniqueN(COLOR) > 1L
  ),
  by = source_zone
][
  V1 == TRUE,
  .N
] > 0L) {
  stop(
    "A source ecotype has more than one colour in the lookup table."
  )
}

population_palette <- unique(
  population_palette,
  by = "source_zone"
)

if (anyNA(
  population_palette$COLOR
) || any(
  !nzchar(
    population_palette$COLOR
  )
)) {
  stop(
    "At least one source ecotype lacks a colour."
  )
}

make_plot_template <- function(reference_map) {
  aggregation_factor <- max(
    1L,
    as.integer(
      ceiling(
        sqrt(
          ncell(reference_map) /
            target_plot_cells
        )
      )
    )
  )
  
  if (aggregation_factor == 1L) {
    return(reference_map)
  }
  
  aggregate(
    reference_map,
    fact = aggregation_factor,
    fun = "modal",
    na.rm = TRUE
  )
}

keep_modeled_zones <- function(x) {
  subst(
    x,
    from = model_zoneID,
    to = model_zoneID,
    others = NA
  )
}

build_candidate_layers <- function(
    assigned_plot,
    ranked_zone_plot,
    ranked_suitability_plot) {
  assigned_candidate <- keep_modeled_zones(
    assigned_plot
  )
  
  rank_candidates <- vector(
    "list",
    5L
  )
  
  for (rank_index in 1:5) {
    candidate <- ifel(
      ranked_suitability_plot[[rank_index]] >=
        dual_threshold,
      ranked_zone_plot[[rank_index]],
      NA
    )
    
    candidate <- keep_modeled_zones(
      candidate
    )
    
    # The assigned Top-1 ecotype is already drawn. This removes its duplicate
    # from the ranked layers while retaining a different tie-selected ecotype.
    candidate <- ifel(
      !is.na(assigned_candidate) &
        !is.na(candidate) &
        candidate == assigned_candidate,
      NA,
      candidate
    )
    
    rank_candidates[[rank_index]] <- candidate
  }
  
  list(
    `1` = assigned_candidate,
    `3` = do.call(
      c,
      c(
        list(assigned_candidate),
        rank_candidates[1:3]
      )
    ),
    `5` = do.call(
      c,
      c(
        list(assigned_candidate),
        rank_candidates[1:5]
      )
    )
  )
}

plot_population_coloured_species_page <- function(
    background,
    candidate_layers,
    method,
    scenario,
    rank_cutoff,
    output_file) {
  png(
    output_file,
    width = 4200,
    height = 1320,
    res = 250,
    bg = "white"
  )
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  
  par(
    mfrow = c(
      2,
      5
    ),
    mar = c(
      0.15,
      0.15,
      1.72,
      0.15
    ),
    oma = c(
      0.15,
      0.15,
      2.25,
      0.15
    )
  )
  
  transparent_zone_colours <- grDevices::adjustcolor(
    population_palette$COLOR,
    alpha.f = population_alpha
  )
  
  for (species_code in species_order) {
    source_zones <- sort(
      unique(
        population_source[
          Species == species_code,
          source_zone
        ]
      )
    )
    
    plot(
      background,
      col = "grey94",
      legend = FALSE,
      axes = FALSE,
      box = FALSE,
      main = paste0(
        species_code,
        "\nZ",
        paste(
          source_zones,
          collapse = ", "
        )
      ),
      cex.main = 0.72,
      maxcell = target_plot_cells
    )
    
    for (layer_index in seq_len(
      nlyr(candidate_layers)
    )) {
      species_layer <- subst(
        candidate_layers[[layer_index]],
        from = source_zones,
        to = source_zones,
        others = NA
      )
      
      indexed_layer <- subst(
        species_layer,
        from = population_palette$source_zone,
        to = seq_len(
          nrow(population_palette)
        ),
        others = NA
      )
      
      plot(
        indexed_layer,
        add = TRUE,
        col = transparent_zone_colours,
        breaks = seq(
          0.5,
          nrow(population_palette) +
            0.5,
          by = 1
        ),
        legend = FALSE,
        axes = FALSE,
        box = FALSE,
        maxcell = target_plot_cells
      )
    }
  }
  
  mtext(
    paste0(
      "Species niches represented by source-ecotype populations | ",
      method_labels[[method]],
      " | ",
      sub(
        "SSP.*$",
        "",
        scenario
      ),
      " ",
      sub(
        "^.*(SSP[0-9]+)$",
        "\\1",
        scenario
      ),
      " | Top-",
      rank_cutoff
    ),
    outer = TRUE,
    line = 0.72,
    cex = 1.05,
    font = 2
  )
  
  mtext(
    paste0(
      "Colours identify source ecotypes; transparent overlap shows multiple contributing populations. ",
      "Ranked candidates require dual suitability >= ",
      dual_threshold,
      "."
    ),
    outer = TRUE,
    line = -0.25,
    cex = 0.70
  )
  
  invisible(TRUE)
}


# One reusable colour key for all population-coloured species maps.
png(
  output_population_colour_key,
  width = 3600,
  height = 1250,
  res = 250,
  bg = "white"
)

par(
  mar = c(
    0.2,
    0.2,
    1.6,
    0.2
  )
)

plot.new()

legend(
  "center",
  legend = paste0(
    "Zone ",
    population_palette$source_zone,
    ": ",
    substr(
      population_palette$zone_name,
      1,
      48
    )
  ),
  fill = grDevices::adjustcolor(
    population_palette$COLOR,
    alpha.f = 0.72
  ),
  border = NA,
  ncol = 3,
  cex = 0.66,
  bty = "n",
  x.intersp = 0.5,
  y.intersp = 0.82
)

title(
  main = "Source-ecotype colour key for population-coloured species-niche maps",
  cex.main = 1.0
)

dev.off()


reference_map <- rast(
  require_file(
    reference_file
  )
)[[1]]

plot_template <- make_plot_template(
  reference_map
)

reference_plot <- resample(
  reference_map,
  plot_template,
  method = "near"
)

background <- ifel(
  !is.na(reference_plot),
  1,
  NA
)

species_map_index <- list()

for (method in method_order) {
  for (scenario in future_order) {
    assigned_file <- require_file(
      assigned_map_file(
        method,
        scenario
      )
    )
    
    rank_zone_path <- require_file(
      ranked_zone_file(
        method,
        scenario
      )
    )
    
    rank_suitability_path <- require_file(
      ranked_suitability_file(
        method,
        scenario
      )
    )
    
    expected_outputs <- file.path(
      species_map_dir,
      paste0(
        "Figure_var_9_population_coloured_species_niches_",
        method,
        "_",
        scenario,
        "_top",
        rank_cutoffs,
        ".png"
      )
    )
    
    for (output_index in seq_along(
      rank_cutoffs
    )) {
      species_map_index[[
        length(species_map_index) + 1L
      ]] <- data.table(
        method = method,
        method_label = method_labels[[method]],
        scenario = scenario,
        rank_cutoff = rank_cutoffs[[output_index]],
        rank_label = paste0(
          "Top-",
          rank_cutoffs[[output_index]]
        ),
        figure_file = expected_outputs[[output_index]],
        criterion = paste0(
          "assigned Top-1 plus source ecotypes within first k ranks at dual suitability >= ",
          dual_threshold
        ),
        colour_definition = "source-ecotype colour with transparent overlap"
      )
    }
    
    if (
      reuse_existing_species_maps &&
      all(
        file.exists(
          expected_outputs
        )
      )
    ) {
      cat(
        "[REUSE] ",
        method,
        " | ",
        scenario,
        "\n",
        sep = ""
      )
      next
    }
    
    assigned <- rast(
      assigned_file
    )[[1]]
    
    ranked_zone <- rast(
      rank_zone_path
    )[[1:5]]
    
    ranked_suitability <- rast(
      rank_suitability_path
    )[[1:5]]
    
    check_geometry(
      assigned,
      reference_map,
      assigned_file
    )
    
    check_geometry(
      ranked_zone,
      reference_map,
      rank_zone_path
    )
    
    check_geometry(
      ranked_suitability,
      reference_map,
      rank_suitability_path
    )
    
    assigned_plot <- resample(
      assigned,
      plot_template,
      method = "near"
    )
    
    ranked_zone_plot <- resample(
      ranked_zone,
      plot_template,
      method = "near"
    )
    
    ranked_suitability_plot <- resample(
      ranked_suitability,
      plot_template,
      method = "near"
    )
    
    candidates_by_k <- build_candidate_layers(
      assigned_plot,
      ranked_zone_plot,
      ranked_suitability_plot
    )
    
    for (output_index in seq_along(
      rank_cutoffs
    )) {
      k <- rank_cutoffs[[output_index]]
      output_file <- expected_outputs[[output_index]]
      
      if (
        reuse_existing_species_maps &&
        file.exists(
          output_file
        )
      ) {
        next
      }
      
      plot_population_coloured_species_page(
        background = background,
        candidate_layers = candidates_by_k[[
          as.character(k)
        ]],
        method = method,
        scenario = scenario,
        rank_cutoff = k,
        output_file = output_file
      )
    }
    
    rm(
      assigned,
      ranked_zone,
      ranked_suitability,
      assigned_plot,
      ranked_zone_plot,
      ranked_suitability_plot,
      candidates_by_k
    )
    
    gc()
  }
}

species_map_index <- rbindlist(
  species_map_index,
  fill = TRUE
)

missing_species_maps <- species_map_index$figure_file[
  !file.exists(
    species_map_index$figure_file
  )
]

if (length(
  missing_species_maps
) > 0L) {
  stop(
    "Population-coloured species maps are missing after plotting: ",
    paste(
      missing_species_maps,
      collapse = "; "
    )
  )
}

fwrite(
  species_map_index,
  output_species_map_index
)

cat(
  "\nCOMPLETE\n",
  "Reused without recalculation:\n",
  "  Figure_var_8_assigned_species_area_topk.png\n",
  "  12 reference-map chord diagrams (zone/category x RF/MF x Top-1/3/5)\n",
  "Created:\n",
  "  ",
  output_bubble,
  "\n",
  "  ",
  output_bubble_table,
  "\n",
  "  ",
  output_bubble_summary,
  "\n",
  "  ",
  output_population_colour_key,
  "\n",
  "  ",
  nrow(species_map_index),
  " population-coloured species-map pages\n",
  "No RF fitting, projection, or area calculation was rerun.\n",
  sep = ""
)
