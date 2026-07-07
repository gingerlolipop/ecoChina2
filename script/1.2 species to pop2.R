# Species to Populations
# Each species × reference-period ecosystem zone is treated as one population.

library(terra)
library(CEMT)
library(reshape2)
library(stringr)

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

# Preserve the original population definition used in the previous results.
# This includes Zone8 when a species occurs there.
zoneID <- 1:55

# Minimum number of occupied raster cells required to define a population.
min_population_cells <- 10L


# 2. Load the reference-period ecosystem raster ================================

ecotype_file <- file.path(base_dir, "raster/ecosys_ori.tif")

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

# Keep presence records only, following the original file structure.
species_data <- species_data[species_data[, 2] == "y", ]

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
  species_data <- species_data[species_data[, 2] == "y", ]
  
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

# Retain only species × zone populations represented by at least 10 cells.
population_abundance[
  !is.na(population_abundance) &
    population_abundance < min_population_cells
] <- NA

# Remove zones that contain no retained population.
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

# Sort as in the previous table: number of populations from high to low.
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

# Long table: one row for every retained species × zone population.
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

# Natural-log transformation.
heatmap_matrix <- log1p(heatmap_abundance)

# Blank cells represent species that do not form a retained population in a zone.
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
  nsmall = 1
)

max_value <- max(
  heatmap_matrix,
  na.rm = TRUE
)

# Use many color steps to create a refined continuous gradient.
heatmap_colors <- grDevices::colorRampPalette(
  c(
    "#FDE725",
    "#5DC863",
    "#21918C",
    "#3B528B",
    "#440154"
  )
)(501)

col_fun <- circlize::colorRamp2(
  seq(
    0,
    max_value,
    length.out = length(heatmap_colors)
  ),
  heatmap_colors
)

# Fine legend ticks: 20 subdivisions per unit.
legend_at <- unique(
  c(
    seq(
      0,
      floor(max_value / 0.05) * 0.05,
      by = 0.05
    ),
    max_value
  )
)

# Show labels only at whole-number values.
legend_labels <- ifelse(
  abs(legend_at - round(legend_at)) < 1e-8,
  sprintf("%.0f", round(legend_at)),
  ""
)

# Values above this threshold use white text.
text_threshold <- max_value * 0.42

ht <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "ln abundance",
  col = col_fun,
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = grid::gpar(
    col = "white",
    lwd = 1
  ),
  row_names_side = "left",
  row_names_gp = grid::gpar(
    fontsize = 10
  ),
  column_names_gp = grid::gpar(
    fontsize = 10
  ),
  column_names_rot = 45,
  row_title = "Ecotypes",
  row_title_gp = grid::gpar(
    fontsize = 11
  ),
  column_title = "Species",
  column_title_side = "bottom",
  column_title_gp = grid::gpar(
    fontsize = 11
  ),
  heatmap_legend_param = list(
    title = "ln abundance",
    color_bar = "continuous",
    at = legend_at,
    labels = legend_labels,
    title_gp = grid::gpar(
      fontsize = 9
    ),
    labels_gp = grid::gpar(
      fontsize = 7.5
    ),
    legend_height = grid::unit(
      9,
      "cm"
    ),
    legend_width = grid::unit(
      3.5,
      "mm"
    ),
    border = "grey55"
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    value <- heatmap_matrix[i, j]
    
    if (!is.na(value)) {
      
      text_color <- ifelse(
        value >= text_threshold,
        "white",
        "black"
      )
      
      grid::grid.text(
        label_matrix[i, j],
        x = x,
        y = y,
        gp = grid::gpar(
          fontsize = 7.5,
          col = text_color
        )
      )
    }
  }
)

# Show in the RStudio Plots panel.
ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "right"
)

# Save the same heatmap.
heatmap_file <- file.path(
  base_dir,
  "species_zone_population_heatmap.png"
)

png(
  heatmap_file,
  width = 2600,
  height = 2200,
  res = 250
)

ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "right"
)

dev.off()


# 9. Summary ==================================================================

cat(
  "\nCOMPLETE\n",
  "Species processed: ", length(species_names), "\n",
  "Populations retained: ",
  sum(!is.na(population_abundance)), "\n",
  "Zones shown in heatmap: ",
  nrow(heatmap_abundance), "\n",
  "Outputs written to: ", base_dir, "\n",
  sep = ""
)