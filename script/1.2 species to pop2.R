# Species to Populations
# Each species × reference-period ecosystem zone is treated as one population.

library(terra)
library(CEMT)
library(reshape2)
library(stringr)
library(grid)

rm(list = ls())
gc()

base_dir <- "H:/Jing/ecoChina2"
species_dir <- file.path(base_dir, "species data")  # Change only if needed.
setwd(base_dir)


# 1. Species information =======================================================

species_info <- data.frame(
  SpeciesCode = c(
    "querVar", "pinuYun", "lariGme", "phylPub", "lariOlg",
    "robiPse", "chamFor", "cyclLon", "pinuSyl", "saliMat"
  ),
  LatinName = c(
    "Quercus variabilis",
    "Pinus yunnanensis",
    "Larix gmelinii",
    "Phyllostachys pubescens",
    "Larix olgensis var. changpaiensis",
    "Robinia pseudoacacia",
    "Chamaecyparis formosensis",
    "Cyclobalanopsis longinux",
    "Pinus sylvestris var. mongolica",
    "Salix matsudana"
  ),
  CommonName = c(
    "Chinese Cork Oak",
    "Yunnan Pine",
    "Dahurian Larch",
    "Moso Bamboo",
    "Hinggan Larch",
    "Black Locust",
    "Taiwan Cypress",
    "Longinux Oak",
    "Mongolian Pine",
    "Chinese Willow"
  ),
  stringsAsFactors = FALSE
)

species_names <- c(
  "lariGme", "lariOlg", "pinuSyl", "querVar", "robiPse",
  "saliMat", "phylPub", "pinuYun", "chamFor", "cyclLon"
)

# Preserve the original population definition.
# This includes Zone8 when a species occurs there.
zoneID <- 1:55

# Minimum occupied raster cells required to define a population.
min_population_cells <- 10L


# 2. Load the reference-period ecosystem raster ================================

ecotype_file <- file.path(
  base_dir,
  "raster/ecosys_ori.tif"
)

if (!file.exists(ecotype_file)) {
  stop("Missing ecosystem raster: ", ecotype_file)
}

ecotype_raster <- rast(ecotype_file)
names(ecotype_raster) <- "zoneID"

ecotype_raster


# 3. Helper: count occupied raster cells by zone ===============================

count_zone_cells <- function(species_data) {
  
  species_points <- vect(
    species_data,
    geom = c("lon", "lat"),
    crs = "EPSG:4326"
  )
  
  if (!same.crs(species_points, ecotype_raster)) {
    species_points <- project(
      species_points,
      crs(ecotype_raster)
    )
  }
  
  cells <- cellFromXY(
    ecotype_raster,
    crds(species_points)
  )
  
  zones <- extract(
    ecotype_raster,
    species_points
  )$zoneID
  
  # Count each occupied raster cell only once for each species.
  occupied <- unique(
    data.frame(
      cell = cells,
      zone = zones
    )
  )
  
  occupied <- occupied[
    !is.na(occupied$cell) &
      !is.na(occupied$zone) &
      occupied$zone %in% zoneID,
    ,
    drop = FALSE
  ]
  
  table(occupied$zone)
}


# 4. Test one species ===========================================================

species_data <- read.csv(
  file.path(species_dir, "lariGme_coord.csv")
)

head(species_data)

# Keep presence records only.
species_data <- species_data[
  !is.na(species_data[, 2]) &
    species_data[, 2] == "y",
  ,
  drop = FALSE
]

if (nrow(species_data) == 0) {
  
  cat("No presence data for species: lariGme\n")
  
} else {
  
  tabulated_values <- count_zone_cells(species_data)
  
  cat("Presence-cell counts for each zone:\n")
  print(tabulated_values)
}


# 5. Overlay all species with reference-period zones ===========================

results <- matrix(
  NA_integer_,
  nrow = length(zoneID),
  ncol = length(species_names)
)

rownames(results) <- paste0("Zone", zoneID)
colnames(results) <- species_names

for (i in seq_along(species_names)) {
  
  species_file <- file.path(
    species_dir,
    paste0(species_names[i], "_coord.csv")
  )
  
  if (!file.exists(species_file)) {
    stop("Missing species file: ", species_file)
  }
  
  species_data <- read.csv(species_file)
  
  species_data <- species_data[
    !is.na(species_data[, 2]) &
      species_data[, 2] == "y",
    ,
    drop = FALSE
  ]
  
  if (nrow(species_data) == 0) {
    cat("No presence data for species:", species_names[i], "\n")
    next
  }
  
  tabulated_values <- count_zone_cells(species_data)
  
  for (zone in names(tabulated_values)) {
    
    results[
      match(as.numeric(zone), zoneID),
      i
    ] <- as.integer(tabulated_values[zone])
  }
}

# Keep only zones occupied by at least one species.
results <- results[
  apply(
    results,
    1,
    function(row) any(!is.na(row))
  ),
  ,
  drop = FALSE
]

write.csv(
  results,
  file.path(base_dir, "specie-zone-pop.csv"),
  row.names = TRUE
)


# 6. Retain valid populations and calculate Shannon diversity ==================

# Population abundance is the number of occupied raster cells.
population_abundance <- results

# Retain populations represented by at least 10 occupied cells.
population_abundance[
  !is.na(population_abundance) &
    population_abundance < min_population_cells
] <- NA

# Remove zones containing no retained population.
population_abundance <- population_abundance[
  apply(
    population_abundance,
    1,
    function(row) any(!is.na(row))
  ),
  ,
  drop = FALSE
]

shannon_index <- function(x) {
  
  x <- x[
    !is.na(x) &
      is.finite(x) &
      x > 0
  ]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  p <- x / sum(x)
  
  -sum(p * log(p))
}

species_summary <- data.frame(
  SpeciesCode = species_names,
  NumberPopulations = colSums(
    !is.na(population_abundance)
  ),
  H_species = apply(
    population_abundance,
    2,
    shannon_index
  ),
  stringsAsFactors = FALSE
)

species_summary <- merge(
  species_info,
  species_summary,
  by = "SpeciesCode",
  all.y = TRUE,
  sort = FALSE
)

species_summary <- species_summary[
  order(
    -species_summary$NumberPopulations,
    species_summary$SpeciesCode
  ),
  ,
  drop = FALSE
]

species_summary$H_species <- round(
  species_summary$H_species,
  4
)

write.csv(
  species_summary,
  file.path(base_dir, "species_population_summary.csv"),
  row.names = FALSE
)

cat("\nSpecies population summary:\n")
print(species_summary)


# 7. Create population tables ==================================================

df <- as.data.frame(population_abundance)
df$Zone <- rownames(df)

df_long <- na.omit(
  melt(
    df,
    id.vars = "Zone",
    variable.name = "Species",
    value.name = "Population"
  )
)

df_long$PopulationID <- paste(
  df_long$Species,
  df_long$Zone,
  sep = "-"
)

write.csv(
  df_long[
    ,
    c(
      "PopulationID",
      "Species",
      "Zone",
      "Population"
    )
  ],
  file.path(base_dir, "species_zone_population_long.csv"),
  row.names = FALSE
)

wide_df <- dcast(
  df_long,
  PopulationID ~ Zone,
  value.var = "Population"
)

zone_cols <- setdiff(
  names(wide_df),
  "PopulationID"
)

wide_df[zone_cols] <- lapply(
  wide_df[zone_cols],
  function(x) {
    x[is.na(x)] <- 0
    x
  }
)

sorted_columns <- zone_cols[
  order(
    as.numeric(
      str_extract(zone_cols, "\\d+")
    )
  )
]

wide_df <- wide_df[
  ,
  c(
    "PopulationID",
    sorted_columns
  ),
  drop = FALSE
]

wide_df$Species <- sub(
  "-Zone.*",
  "",
  wide_df$PopulationID
)

wide_df$ZoneNum <- as.numeric(
  sub(
    ".*-Zone",
    "",
    wide_df$PopulationID
  )
)

wide_df <- wide_df[
  order(
    wide_df$Species,
    wide_df$ZoneNum
  ),
  ,
  drop = FALSE
]

wide_df$Species <- NULL
wide_df$ZoneNum <- NULL

write.csv(
  wide_df,
  file.path(base_dir, "sorted_population_zone.csv"),
  row.names = FALSE
)


# 8. Generate the species × zone abundance heatmap =============================

heatmap_abundance <- population_abundance[
  ,
  species_names,
  drop = FALSE
]

heatmap_abundance <- heatmap_abundance[
  apply(
    heatmap_abundance,
    1,
    function(row) any(!is.na(row))
  ),
  ,
  drop = FALSE
]

heatmap_abundance <- as.matrix(heatmap_abundance)

if (
  nrow(heatmap_abundance) == 0 ||
  ncol(heatmap_abundance) == 0 ||
  all(is.na(heatmap_abundance))
) {
  stop(
    "The heatmap is empty. Check whether min_population_cells = ",
    min_population_cells,
    " removed all populations."
  )
}

# Natural logarithm: ln(abundance + 1).
heatmap_matrix <- log1p(
  heatmap_abundance
)

label_matrix <- matrix(
  "",
  nrow = nrow(heatmap_matrix),
  ncol = ncol(heatmap_matrix),
  dimnames = dimnames(heatmap_matrix)
)

label_matrix[!is.na(heatmap_matrix)] <- format(
  round(
    heatmap_matrix[!is.na(heatmap_matrix)],
    2
  ),
  trim = TRUE,
  nsmall = 2
)

max_value <- max(
  heatmap_matrix,
  na.rm = TRUE
)

if (!is.finite(max_value) || max_value <= 0) {
  stop("The heatmap contains no positive finite values.")
}

# Round the legend maximum upward.
# For example, 10.99 becomes 11.
legend_max <- ceiling(max_value)

# Detailed color gradient.
heatmap_colors <- grDevices::colorRampPalette(
  c(
    "#FDE725",
    "#5DC863",
    "#21918C",
    "#3B528B",
    "#440154"
  )
)(2001)

value_to_color <- function(x) {
  
  color_index <- 1L + round(
    pmax(
      0,
      pmin(legend_max, x)
    ) /
      legend_max *
      (length(heatmap_colors) - 1L)
  )
  
  heatmap_colors[color_index]
}


# Draw heatmap and custom legend ===============================================

draw_species_heatmap <- function() {
  
  nr <- nrow(heatmap_matrix)
  nc <- ncol(heatmap_matrix)
  
  grid::grid.newpage()
  
  # Heatmap and legend share the same layout row and therefore the same height.
  plot_layout <- grid::grid.layout(
    nrow = 5,
    ncol = 6,
    widths = grid::unit.c(
      grid::unit(0.55, "cm"),
      grid::unit(1.55, "cm"),
      grid::unit(1, "null"),
      grid::unit(0.55, "cm"),
      grid::unit(1.55, "cm"),
      grid::unit(0.20, "cm")
    ),
    heights = grid::unit.c(
      grid::unit(0.25, "cm"),
      grid::unit(1, "null"),
      grid::unit(1.35, "cm"),
      grid::unit(0.50, "cm"),
      grid::unit(0.20, "cm")
    )
  )
  
  grid::pushViewport(
    grid::viewport(
      layout = plot_layout,
      gp = grid::gpar(fill = "white")
    )
  )
  
  
  # 8.1 Heatmap body -----------------------------------------------------------
  
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 2,
      layout.pos.col = 3,
      clip = "off"
    )
  )
  
  for (i in seq_len(nr)) {
    
    for (j in seq_len(nc)) {
      
      value <- heatmap_matrix[i, j]
      
      x_position <- (j - 0.5) / nc
      y_position <- 1 - (i - 0.5) / nr
      
      if (is.na(value)) {
        fill_color <- "white"
      } else {
        fill_color <- value_to_color(value)
      }
      
      grid::grid.rect(
        x = grid::unit(x_position, "npc"),
        y = grid::unit(y_position, "npc"),
        width = grid::unit(1 / nc, "npc"),
        height = grid::unit(1 / nr, "npc"),
        gp = grid::gpar(
          fill = fill_color,
          col = "white",
          lwd = 0.7
        )
      )
      
      if (!is.na(value)) {
        
        rgb_value <- grDevices::col2rgb(fill_color) / 255
        
        luminance <- (
          0.2126 * rgb_value[1] +
            0.7152 * rgb_value[2] +
            0.0722 * rgb_value[3]
        )
        
        text_color <- ifelse(
          luminance < 0.52,
          "white",
          "black"
        )
        
        grid::grid.text(
          label_matrix[i, j],
          x = grid::unit(x_position, "npc"),
          y = grid::unit(y_position, "npc"),
          gp = grid::gpar(
            fontsize = 7.3,
            col = text_color
          )
        )
      }
    }
  }
  
  grid::popViewport()
  
  
  # 8.2 Row names --------------------------------------------------------------
  
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 2,
      layout.pos.col = 2,
      clip = "off"
    )
  )
  
  grid::grid.text(
    rownames(heatmap_matrix),
    x = grid::unit(0.98, "npc"),
    y = grid::unit(
      1 - (seq_len(nr) - 0.5) / nr,
      "npc"
    ),
    just = "right",
    gp = grid::gpar(
      fontsize = 8.5
    )
  )
  
  grid::popViewport()
  
  
  # 8.3 Y-axis title -----------------------------------------------------------
  
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 2,
      layout.pos.col = 1
    )
  )
  
  grid::grid.text(
    "Ecotypes",
    rot = 90,
    gp = grid::gpar(
      fontsize = 10
    )
  )
  
  grid::popViewport()
  
  
  # 8.4 Column names -----------------------------------------------------------
  
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 3,
      layout.pos.col = 3,
      clip = "off"
    )
  )
  
  grid::grid.text(
    colnames(heatmap_matrix),
    x = grid::unit(
      (seq_len(nc) - 0.5) / nc,
      "npc"
    ),
    y = grid::unit(0.95, "npc"),
    rot = 45,
    just = c("right", "centre"),
    gp = grid::gpar(
      fontsize = 8.5
    )
  )
  
  grid::popViewport()
  
  
  # 8.5 X-axis title -----------------------------------------------------------
  
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 4,
      layout.pos.col = 3
    )
  )
  
  grid::grid.text(
    "Species",
    gp = grid::gpar(
      fontsize = 10
    )
  )
  
  grid::popViewport()
  
  
  # 8.6 Legend title -----------------------------------------------------------
  
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 2,
      layout.pos.col = 4
    )
  )
  
  grid::grid.text(
    "ln abundance",
    rot = 90,
    gp = grid::gpar(
      fontsize = 9
    )
  )
  
  grid::popViewport()
  
  
  # 8.7 Continuous legend with fine ticks -------------------------------------
  
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 2,
      layout.pos.col = 5,
      xscale = c(0, 1),
      yscale = c(0, legend_max),
      clip = "off"
    )
  )
  
  # Thin sections produce a smooth continuous gradient.
  number_legend_sections <- 1000L
  
  legend_edges <- seq(
    0,
    legend_max,
    length.out = number_legend_sections + 1L
  )
  
  legend_centres <- (
    legend_edges[-1] +
      legend_edges[-length(legend_edges)]
  ) / 2
  
  legend_section_height <- legend_max /
    number_legend_sections
  
  grid::grid.rect(
    x = grid::unit(0.18, "npc"),
    y = grid::unit(legend_centres, "native"),
    width = grid::unit(0.24, "npc"),
    height = grid::unit(
      legend_section_height * 1.02,
      "native"
    ),
    gp = grid::gpar(
      fill = value_to_color(legend_centres),
      col = NA
    )
  )
  
  # Legend border.
  grid::grid.rect(
    x = grid::unit(0.18, "npc"),
    y = grid::unit(legend_max / 2, "native"),
    width = grid::unit(0.24, "npc"),
    height = grid::unit(legend_max, "native"),
    gp = grid::gpar(
      fill = NA,
      col = "grey45",
      lwd = 0.7
    )
  )
  
  # Twenty subdivisions per unit.
  legend_step <- 0.05
  
  minor_ticks <- seq(
    0,
    legend_max,
    by = legend_step
  )
  
  major_ticks <- seq(
    0,
    legend_max,
    by = 1
  )
  
  is_major <- abs(
    minor_ticks - round(minor_ticks)
  ) < 1e-8
  
  fine_ticks <- minor_ticks[
    !is_major
  ]
  
  # Short minor ticks.
  grid::grid.segments(
    x0 = grid::unit(0.30, "npc"),
    x1 = grid::unit(0.35, "npc"),
    y0 = grid::unit(fine_ticks, "native"),
    y1 = grid::unit(fine_ticks, "native"),
    gp = grid::gpar(
      col = "grey55",
      lwd = 0.45
    )
  )
  
  # Longer integer ticks.
  grid::grid.segments(
    x0 = grid::unit(0.30, "npc"),
    x1 = grid::unit(0.41, "npc"),
    y0 = grid::unit(major_ticks, "native"),
    y1 = grid::unit(major_ticks, "native"),
    gp = grid::gpar(
      col = "grey20",
      lwd = 0.7
    )
  )
  
  # Integer labels from 0 to the rounded upper limit.
  grid::grid.text(
    as.character(major_ticks),
    x = grid::unit(0.46, "npc"),
    y = grid::unit(major_ticks, "native"),
    just = "left",
    gp = grid::gpar(
      fontsize = 7.5
    )
  )
  
  grid::popViewport()
  
  grid::popViewport()
}


# Display directly in the RStudio Plots panel.
draw_species_heatmap()


# Save the same heatmap.
heatmap_file <- file.path(
  base_dir,
  "species_zone_population_heatmap.png"
)

png(
  filename = heatmap_file,
  width = 2600,
  height = 2200,
  res = 250
)

draw_species_heatmap()

dev.off()

# 9. Summary ===================================================================

cat(
  "\nCOMPLETE\n",
  "Species processed: ",
  length(species_names),
  "\n",
  "Populations retained: ",
  sum(!is.na(population_abundance)),
  "\n",
  "Zones shown in heatmap: ",
  nrow(heatmap_abundance),
  "\n",
  "Heatmap saved to: ",
  heatmap_file,
  "\n",
  "Outputs written to: ",
  base_dir,
  "\n",
  sep = ""
)