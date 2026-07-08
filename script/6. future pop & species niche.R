# Future population and species niches from available projected ecosystem maps
# Each retained species × reference-zone combination is one population.
# Future population niche = pixels assigned to that source zone.
# Future species niche = union of all population niches within the species.
#
# This version is resume-friendly:
#   1. Missing future assigned-zone maps are skipped rather than stopping.
#   2. Available method × scenario maps are processed immediately.
#   3. Existing population/species outputs can be reused on later reruns.

library(terra)
library(data.table)

rm(list = ls())
gc()

# 0. Paths and settings =========================================================

base_dir <- "H:/Jing/ecoChina2"

result_map_root <- file.path(base_dir, "result maps")
population_file <- file.path(
  base_dir,
  "species_zone_population_long.csv"
)
species_data_dir <- file.path(base_dir, "species data")
palette_file <- file.path(base_dir, "color_palette_China.csv")
reference_file <- file.path(base_dir, "raster/ecosys_ori.tif")

output_root <- file.path(base_dir, "future tree niche")
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

method_order <- c(
  "optimized_mf",
  "optimized_rf",
  "plain_mf",
  "plain_rf"
)

future_order <- c(
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

model_zoneID <- c(1:7, 9:50, 52:55)

threshold <- 0.2
tie_tol <- 1e-4
novel_value <- 99

# Reuse previously generated population/species rasters on later reruns.
reuse_existing_outputs <- TRUE

set.seed(49)

# Padding (degrees) around the focal species occurrence extent in the main panel.
focus_pad_deg <- 2

# Plot styling.
point_cex_main <- 0.65
point_cex_inset <- 0.30
point_alpha <- 0.30
legend_cex <- 0.62
legend_inset_x <- 0.985
legend_inset_y <- 0.985


# 1. Read population and palette tables ========================================

if (!file.exists(population_file)) {
  stop("Missing population table: ", population_file)
}

if (!file.exists(palette_file)) {
  stop("Missing color palette: ", palette_file)
}

if (!file.exists(reference_file)) {
  stop("Missing reference ecosystem raster: ", reference_file)
}

pop <- fread(population_file)

required_pop_cols <- c(
  "PopulationID",
  "Species",
  "Zone",
  "Population"
)

if (!all(required_pop_cols %in% names(pop))) {
  stop(
    "Population table must contain: ",
    paste(required_pop_cols, collapse = ", ")
  )
}

pop[, `:=`(
  PopulationID = as.character(PopulationID),
  Species = as.character(Species),
  Zone = as.character(Zone),
  source_zone = as.integer(sub("^Zone", "", Zone)),
  reference_abundance = as.numeric(Population)
)]

if (anyNA(pop$source_zone)) {
  stop("Some population zones could not be converted to numeric zone IDs.")
}

if (anyDuplicated(pop[, .(Species, source_zone)])) {
  stop("Duplicated species × source-zone populations were found.")
}

palette <- fread(palette_file)

if (!all(c("zoneID", "zone", "COLOR") %in% names(palette))) {
  stop("Palette must contain zoneID, zone and COLOR columns.")
}

palette[, zoneID := as.integer(zoneID)]

pop <- merge(
  pop,
  palette[, .(source_zone = zoneID, zone_name = zone, COLOR)],
  by = "source_zone",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(pop$COLOR)) {
  stop(
    "Missing palette colors for zones: ",
    paste(unique(pop[is.na(COLOR), source_zone]), collapse = ", ")
  )
}

# Future assigned maps do not contain Zones 8 or 51.
pop[, projected := source_zone %in% model_zoneID]

fwrite(
  pop,
  file.path(table_dir, "population_projection_lookup.csv")
)

fwrite(
  pop[projected == FALSE],
  file.path(table_dir, "population_projection_exclusions.csv")
)

projected_pop <- pop[projected == TRUE]

if (!nrow(projected_pop)) {
  stop("No populations occur in zones represented by the future ecosystem maps.")
}

species_names <- sort(unique(projected_pop$Species))


# 2. Find available future ecosystem maps ======================================

future_map_file <- function(method, scenario) {
  file.path(
    result_map_root,
    method,
    paste0(
      "assigned_zone_", scenario,
      "_threshold", threshold,
      "_tol", tie_tol,
      "_novel", novel_value,
      "_maskNA8_noNovelNormal.tif"
    )
  )
}

all_jobs <- CJ(
  method = method_order,
  scenario = future_order,
  sorted = FALSE
)

all_jobs[, map_file := mapply(
  future_map_file,
  method,
  scenario,
  USE.NAMES = FALSE
)]

available_jobs <- all_jobs[file.exists(map_file)]
missing_jobs <- all_jobs[!file.exists(map_file)]

fwrite(
  available_jobs,
  file.path(
    table_dir,
    "future_ecosystem_map_jobs_available.csv"
  )
)

fwrite(
  missing_jobs,
  file.path(
    table_dir,
    "future_ecosystem_map_jobs_missing.csv"
  )
)

if (nrow(missing_jobs) > 0) {
  cat(
    "\n[SKIP MISSING FUTURE MAPS]\n",
    paste(
      paste0(
        missing_jobs$method,
        " | ",
        missing_jobs$scenario,
        " | ",
        missing_jobs$map_file
      ),
      collapse = "\n"
    ),
    "\n",
    sep = ""
  )
}

if (nrow(available_jobs) == 0) {
  stop(
    "No future assigned-zone maps are currently available.\n",
    "Missing-map list: ",
    file.path(
      table_dir,
      "future_ecosystem_map_jobs_missing.csv"
    )
  )
}

jobs <- available_jobs

cat(
  "\n[AVAILABLE FUTURE MAPS]\n",
  "Available jobs: ", nrow(jobs), " of ", nrow(all_jobs), "\n",
  sep = ""
)

print(jobs[, .(method, scenario)])


# 3. Plot all species in the normal/reference period ============================

reference_map <- rast(reference_file)
names(reference_map) <- "zoneID"

reference_mask <- ifel(
  !is.na(reference_map),
  1,
  NA
)

normal_figure_dir <- file.path(
  figure_dir,
  "normal period"
)

dir.create(
  normal_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

normal_plot_index <- list()

for (species_name in species_names) {
  
  cat(
    "\n[NORMAL-PERIOD FIGURE]",
    species_name,
    "\n"
  )
  
  species_file <- file.path(
    species_data_dir,
    paste0(species_name, "_coord.csv")
  )
  
  if (!file.exists(species_file)) {
    warning(
      "Skipping normal-period figure; missing species file: ",
      species_file
    )
    next
  }
  
  species_points <- fread(species_file)
  
  if (!all(c("lon", "lat") %in% names(species_points))) {
    warning(
      "Skipping normal-period figure; lon/lat columns are missing: ",
      species_file
    )
    next
  }
  
  # Preserve the original presence rule: the second column equals "y".
  species_points <- species_points[
    species_points[[2]] == "y" &
      complete.cases(species_points[, .(lon, lat)])
  ]
  
  if (!nrow(species_points)) {
    warning(
      "Skipping normal-period figure; no presence records: ",
      species_name
    )
    next
  }
  
  points_v <- vect(
    species_points,
    geom = c("lon", "lat"),
    crs = "EPSG:4326"
  )
  
  if (!same.crs(points_v, reference_map)) {
    points_v <- project(
      points_v,
      crs(reference_map)
    )
  }
  
  point_cells <- cellFromXY(
    reference_map,
    crds(points_v)
  )
  
  point_zone <- extract(
    reference_map,
    points_v
  )$zoneID
  
  point_dt <- data.table(
    cell = point_cells,
    source_zone = as.integer(point_zone)
  )
  
  # Match the population definition:
  # one occupied raster cell counts once.
  point_dt <- unique(
    point_dt[
      !is.na(cell) &
        source_zone %in%
        projected_pop[
          Species == species_name,
          source_zone
        ]
    ],
    by = c("cell", "source_zone")
  )
  
  point_dt <- merge(
    point_dt,
    projected_pop[
      Species == species_name,
      .(
        source_zone,
        PopulationID,
        zone_name,
        COLOR
      )
    ],
    by = "source_zone",
    all.x = TRUE,
    sort = FALSE
  )
  
  if (!nrow(point_dt)) {
    warning(
      "Skipping normal-period figure; no retained population cells: ",
      species_name
    )
    next
  }
  
  point_xy <- xyFromCell(
    reference_map,
    point_dt$cell
  )
  
  point_dt[, `:=`(
    x = point_xy[, 1],
    y = point_xy[, 2],
    Species = species_name
  )]
  
  point_dt[, `:=`(
    COLOR_main = grDevices::adjustcolor(COLOR, alpha.f = point_alpha),
    COLOR_inset = grDevices::adjustcolor(COLOR, alpha.f = point_alpha)
  )]
  
  point_table_file <- file.path(
    table_dir,
    paste0(
      "normal_population_points_",
      species_name,
      ".csv"
    )
  )
  
  fwrite(
    point_dt,
    point_table_file
  )
  
  point_plot <- vect(
    point_dt,
    geom = c("x", "y"),
    crs = crs(reference_map)
  )
  
  # Main panel: occurrence extent plus two degrees on every side.
  map_ext <- ext(reference_map)
  
  focus_ext <- ext(
    max(xmin(map_ext), min(point_dt$x) - focus_pad_deg),
    min(xmax(map_ext), max(point_dt$x) + focus_pad_deg),
    max(ymin(map_ext), min(point_dt$y) - focus_pad_deg),
    min(ymax(map_ext), max(point_dt$y) + focus_pad_deg)
  )
  
  focus_mask <- crop(
    reference_mask,
    focus_ext,
    snap = "out"
  )
  
  legend_dt <- unique(
    point_dt[
      order(source_zone),
      .(
        PopulationID,
        source_zone,
        COLOR
      )
    ]
  )
  
  plot_normal_species <- function() {
    
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    
    # Main panel: focal occurrence region.
    par(
      fig = c(0, 1, 0, 1),
      mar = c(3.2, 3.2, 3.2, 1.2),
      new = FALSE
    )
    
    plot(
      focus_mask,
      col = "grey94",
      legend = FALSE,
      axes = TRUE,
      main = paste0(
        species_name,
        ": normal-period populations"
      )
    )
    
    plot(
      point_plot,
      add = TRUE,
      pch = 16,
      col = point_dt$COLOR_main,
      cex = point_cex_main
    )
    
    legend(
      x = "topright",
      inset = c(1 - legend_inset_x, 1 - legend_inset_y),
      legend = legend_dt$PopulationID,
      pch = 16,
      col = grDevices::adjustcolor(legend_dt$COLOR, alpha.f = point_alpha),
      cex = legend_cex,
      bty = "n",
      title = "Population",
      xpd = FALSE
    )
    
    # Inset panel: full China and the focal extent.
    par(
      fig = c(0.06, 0.23, 0.09, 0.26),
      mar = c(0.5, 0.5, 0.5, 0.5),
      new = TRUE
    )
    
    plot(
      reference_mask,
      col = "grey94",
      legend = FALSE,
      axes = FALSE,
      main = ""
    )
    
    plot(
      point_plot,
      add = TRUE,
      pch = 16,
      col = point_dt$COLOR_inset,
      cex = point_cex_inset
    )
    
    rect(
      xleft = xmin(focus_ext),
      ybottom = ymin(focus_ext),
      xright = xmax(focus_ext),
      ytop = ymax(focus_ext),
      border = "black",
      lwd = 1.5
    )
    
    box(
      col = "grey35",
      lwd = 1
    )
  }
  
  # Display each species sequentially in the RStudio plot history.
  plot_normal_species()
  
  normal_figure_file <- file.path(
    normal_figure_dir,
    paste0(
      species_name,
      "_normal_population_points.png"
    )
  )
  
  png(
    normal_figure_file,
    width = 2200,
    height = 1600,
    res = 250
  )
  
  plot_normal_species()
  
  dev.off()
  
  normal_plot_index[[length(normal_plot_index) + 1L]] <- data.table(
    Species = species_name,
    n_populations = nrow(legend_dt),
    n_occurrence_cells = nrow(point_dt),
    figure_file = normal_figure_file,
    point_table_file = point_table_file
  )
  
  rm(
    species_points,
    points_v,
    point_dt,
    point_plot,
    focus_mask,
    legend_dt
  )
  
  gc()
}

normal_plot_index <- rbindlist(
  normal_plot_index,
  fill = TRUE
)

fwrite(
  normal_plot_index,
  file.path(
    table_dir,
    "normal_population_figure_index.csv"
  )
)


# 4. Generate and plot future population and species niche maps ================

# Population raster:
#   cell value = source-zone ID of the suitable population.
#
# Species raster:
#   1  = suitable for at least one population;
#   NA = unsuitable for all populations.

future_population_figure_root <- file.path(
  figure_dir,
  "future period",
  "population niche"
)

future_species_figure_root <- file.path(
  figure_dir,
  "future period",
  "species niche"
)

dir.create(
  future_population_figure_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  future_species_figure_root,
  recursive = TRUE,
  showWarnings = FALSE
)

population_area_results <- list()
species_area_results <- list()
future_figure_index <- list()

result_index <- 0L

for (j in seq_len(nrow(jobs))) {
  
  method <- jobs$method[j]
  scenario <- jobs$scenario[j]
  map_file <- jobs$map_file[j]
  
  cat(
    "\n==============================\n",
    "FUTURE TREE NICHE: ",
    method,
    " | ",
    scenario,
    "\n",
    "==============================\n",
    sep = ""
  )
  
  future_map <- rast(map_file)
  names(future_map) <- "future_zone"
  
  if (!compareGeom(
    future_map,
    reference_map,
    stopOnError = FALSE
  )) {
    stop("Geometry mismatch: ", map_file)
  }
  
  # Each raster cell receives its own actual surface area in km2.
  cell_area <- cellSize(
    future_map,
    unit = "km"
  )
  
  population_dir <- file.path(
    output_root,
    method,
    scenario,
    "population niche"
  )
  
  species_dir_out <- file.path(
    output_root,
    method,
    scenario,
    "species niche"
  )
  
  dir.create(
    population_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    species_dir_out,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  for (species_name in species_names) {
    
    species_pop <- projected_pop[
      Species == species_name
    ][order(source_zone)]
    
    pop_zones <- species_pop$source_zone
    
    population_output <- file.path(
      population_dir,
      paste0(
        species_name,
        "_population_niche.tif"
      )
    )
    
    species_output <- file.path(
      species_dir_out,
      paste0(
        species_name,
        "_species_niche.tif"
      )
    )
    
    reuse_this_output <- (
      reuse_existing_outputs &&
        file.exists(population_output) &&
        file.exists(species_output)
    )
    
    if (reuse_this_output) {
      
      population_map <- rast(population_output)
      species_map <- rast(species_output)
      
      valid_existing <- (
        compareGeom(
          population_map,
          future_map,
          stopOnError = FALSE
        ) &&
          compareGeom(
            species_map,
            future_map,
            stopOnError = FALSE
          )
      )
      
      if (valid_existing) {
        names(population_map) <- "source_zone"
        names(species_map) <- "species_suitable"
        
        cat(
          "[REUSE]",
          species_name,
          "|",
          method,
          "|",
          scenario,
          "\n"
        )
        
      } else {
        reuse_this_output <- FALSE
        rm(population_map, species_map)
        gc()
      }
    }
    
    if (!reuse_this_output) {
      
      # Retain only future ecosystem zones that host a reference population
      # of the focal species. Each value identifies one population.
      population_map <- subst(
        future_map,
        from = pop_zones,
        to = pop_zones,
        others = NA
      )
      
      names(population_map) <- "source_zone"
      
      writeRaster(
        population_map,
        population_output,
        overwrite = TRUE,
        wopt = list(
          datatype = "INT2S",
          gdal = "COMPRESS=LZW"
        )
      )
      
      # Union of all population niches for the focal species.
      species_map <- ifel(
        !is.na(population_map),
        1,
        NA
      )
      
      names(species_map) <- "species_suitable"
      
      writeRaster(
        species_map,
        species_output,
        overwrite = TRUE,
        wopt = list(
          datatype = "INT1U",
          gdal = "COMPRESS=LZW"
        )
      )
    }
    
    # Population-level future suitable area.
    pop_area <- zonal(
      cell_area,
      population_map,
      fun = "sum",
      na.rm = TRUE
    )
    
    if (nrow(pop_area)) {
      
      names(pop_area) <- c(
        "source_zone",
        "future_area_km2"
      )
      
      pop_area <- as.data.table(pop_area)
      
      pop_area <- merge(
        species_pop[
          ,
          .(
            PopulationID,
            Species,
            source_zone,
            zone_name,
            COLOR,
            reference_abundance
          )
        ],
        pop_area,
        by = "source_zone",
        all.x = TRUE,
        sort = FALSE
      )
      
    } else {
      
      pop_area <- species_pop[
        ,
        .(
          PopulationID,
          Species,
          source_zone,
          zone_name,
          COLOR,
          reference_abundance,
          future_area_km2 = 0
        )
      ]
    }
    
    pop_area[
      is.na(future_area_km2),
      future_area_km2 := 0
    ]
    
    pop_area[, `:=`(
      method = method,
      scenario = scenario,
      population_raster = population_output
    )]
    
    result_index <- result_index + 1L
    population_area_results[[result_index]] <- pop_area
    
    # Species-level future suitable area.
    species_area_km2 <- global(
      cell_area * species_map,
      fun = "sum",
      na.rm = TRUE
    )[1, 1]
    
    species_area_results[[result_index]] <- data.table(
      Species = species_name,
      method = method,
      scenario = scenario,
      populations_projected = length(pop_zones),
      future_area_km2 = species_area_km2,
      species_raster = species_output
    )
    
    # Future figures ------------------------------------------------------------
    # Plot the full China map for both future population niche and future species niche.
    
    suitable_cells <- which(
      !is.na(values(species_map, mat = FALSE))
    )
    
    if (length(suitable_cells) > 0) {
      
      species_pop_legend <- species_pop[
        source_zone %in%
          unique(
            as.integer(
              na.omit(
                values(
                  population_map,
                  mat = FALSE
                )
              )
            )
          ),
        .(
          PopulationID,
          source_zone,
          COLOR
        )
      ][order(source_zone)]
      
      species_pop_legend[, COLOR_alpha := grDevices::adjustcolor(
        COLOR,
        alpha.f = point_alpha
      )]
      
      # Recode the categorical population map to consecutive indices so that
      # colors always match the actual source-zone populations.
      pop_index <- seq_len(nrow(species_pop_legend))
      
      population_full_index <- subst(
        population_map,
        from = species_pop_legend$source_zone,
        to = pop_index,
        others = NA
      )
      
      plot_future_population <- function() {
        
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)
        
        par(
          mar = c(3.2, 3.2, 3.2, 1.2)
        )
        
        plot(
          reference_mask,
          col = "grey94",
          legend = FALSE,
          axes = TRUE,
          main = paste0(
            species_name,
            ": future population niches\n",
            method,
            " | ",
            scenario
          )
        )
        
        plot(
          population_full_index,
          add = TRUE,
          col = species_pop_legend$COLOR_alpha,
          breaks = seq(
            0.5,
            nrow(species_pop_legend) + 0.5,
            by = 1
          ),
          legend = FALSE
        )
        
        legend(
          x = "topright",
          inset = c(1 - legend_inset_x, 1 - legend_inset_y),
          legend = species_pop_legend$PopulationID,
          pch = 15,
          col = species_pop_legend$COLOR_alpha,
          cex = legend_cex,
          bty = "n",
          title = "Population",
          xpd = FALSE
        )
      }
      
      plot_future_species <- function() {
        
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)
        
        par(
          mar = c(3.2, 3.2, 3.2, 1.2)
        )
        
        plot(
          reference_mask,
          col = "grey94",
          legend = FALSE,
          axes = TRUE,
          main = paste0(
            species_name,
            ": future species niche\n",
            method,
            " | ",
            scenario
          )
        )
        
        plot(
          species_map,
          add = TRUE,
          col = grDevices::adjustcolor("#2E8B57", alpha.f = point_alpha),
          legend = FALSE
        )
      }
      
      future_population_figure_dir <- file.path(
        future_population_figure_root,
        method,
        scenario
      )
      
      future_species_figure_dir <- file.path(
        future_species_figure_root,
        method,
        scenario
      )
      
      dir.create(
        future_population_figure_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      dir.create(
        future_species_figure_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      future_population_figure <- file.path(
        future_population_figure_dir,
        paste0(
          species_name,
          "_future_population_niche.png"
        )
      )
      
      future_species_figure <- file.path(
        future_species_figure_dir,
        paste0(
          species_name,
          "_future_species_niche.png"
        )
      )
      
      # Display the future population figure in the RStudio plot history.
      plot_future_population()
      
      png(
        future_population_figure,
        width = 2200,
        height = 1600,
        res = 250
      )
      
      plot_future_population()
      
      dev.off()
      
      # Display the future species figure after the population figure.
      plot_future_species()
      
      png(
        future_species_figure,
        width = 2200,
        height = 1600,
        res = 250
      )
      
      plot_future_species()
      
      dev.off()
      
      future_figure_index[[length(future_figure_index) + 1L]] <- data.table(
        Species = species_name,
        method = method,
        scenario = scenario,
        n_populations = nrow(species_pop_legend),
        population_figure = future_population_figure,
        species_figure = future_species_figure
      )
      
      rm(
        suitable_cells,
        species_pop_legend,
        population_full_index
      )
      
      gc()
      
    } else {
      
      warning(
        "No future suitable cells for ",
        species_name,
        " | ",
        method,
        " | ",
        scenario
      )
    }
    
    cat(
      "[SAVED]",
      species_name,
      "| populations:",
      length(pop_zones),
      "| area km2:",
      round(species_area_km2, 2),
      "\n"
    )
    
    rm(
      population_map,
      species_map,
      pop_area
    )
    
    gc()
  }
  
  rm(
    future_map,
    cell_area
  )
  
  gc()
}


# 5. Save summary tables for currently available jobs ===========================

population_area <- rbindlist(
  population_area_results,
  fill = TRUE
)

species_area <- rbindlist(
  species_area_results,
  fill = TRUE
)

setcolorder(
  population_area,
  c(
    "Species",
    "PopulationID",
    "source_zone",
    "zone_name",
    "reference_abundance",
    "method",
    "scenario",
    "future_area_km2",
    "COLOR",
    "population_raster"
  )
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

fwrite(
  population_area,
  file.path(
    table_dir,
    "future_population_niche_area.csv"
  )
)

fwrite(
  species_area,
  file.path(
    table_dir,
    "future_species_niche_area.csv"
  )
)

future_figure_index <- rbindlist(
  future_figure_index,
  fill = TRUE
)

fwrite(
  future_figure_index,
  file.path(
    table_dir,
    "future_niche_figure_index.csv"
  )
)

cat(
  "\nCOMPLETE\n",
  "Available method-scenario jobs processed: ",
  nrow(jobs),
  " of ",
  nrow(all_jobs),
  "\n",
  "Methods represented: ",
  paste(unique(jobs$method), collapse = ", "),
  "\n",
  "Scenarios represented: ",
  paste(unique(jobs$scenario), collapse = ", "),
  "\n",
  "Species: ",
  length(species_names),
  "\n",
  "Population map definition: categorical source-zone raster\n",
  "Species map definition: union of all population niches\n",
  "Normal-period figures: ",
  normal_figure_dir,
  "\n",
  "Future-period figures: ",
  file.path(figure_dir, "future period"),
  "\n",
  "Available-job table: ",
  file.path(
    table_dir,
    "future_ecosystem_map_jobs_available.csv"
  ),
  "\n",
  "Missing-job table: ",
  file.path(
    table_dir,
    "future_ecosystem_map_jobs_missing.csv"
  ),
  "\n",
  "Outputs: ",
  output_root,
  "\n",
  sep = ""
)
