# China vegetation-zone color palette using terra
#
# This script creates:
#   1. color_palette_China.json
#   2. color_palette_China.csv
#   3. map and legend examples in "result maps/"
#
# Important:
#   - Colors are matched by the real raster value zoneID.
#   - zoneID 8 is kept in the palette as a safeguard, but never shown in legends.
#   - zoneID 99 is novel ecosystem / no suitable current zone.
#   - Zone names are hard-coded here, so the name legend does not depend on
#     unstable JSON/CSV column names.

library(terra)
library(jsonlite)

setwd("H:/Jing/ecoChina2")

# ------------------------------------------------------------
# 1. Hard-coded palette table
# ------------------------------------------------------------

color_df <- data.frame(
  zoneID = c(
    1, 2, 3, 4, 5, 6,
    7, 8, 9, 10, 11, 12,
    13, 14, 15, 16, 17, 18,
    19, 20, 21, 22, 23, 24,
    25, 26, 27, 28, 29, 30,
    31, 32, 33, 34, 35, 36,
    37, 38, 39, 40, 41, 42,
    43, 44, 45, 46, 47, 48,
    49, 50, 51, 52, 53, 54,
    55, 56, 99
  ),
  zone = c(
    "One year one ripe fields with cold-tolerant crops of short growth duration",
    "Needleleaf forests in cold-temperate zone and on mountains in temperate zone",
    "Broadleaf deciduous forests in temperate zone",
    "Grass, Carex and forb swamp meadows",
    "Grass and forb meadows",
    "Cold-temperate and temperate marshes",
    "Deciduous scrubs in temperate zone",
    "Rare / original raster odd value",
    "Needleleaf and deciduous broadleaf mixed forests in temperate zone",
    "Needleleaf forests on mountains in subtropical and tropical zones",
    "Temperate grass and forb meadow steppes",
    "One year one ripe grain fields and cold-tolerant economic crop fields",
    "Temperate tufted grass steppes",
    "Microphyllous deciduous woodlands in temperate zone",
    "Grass and forb halophytic meadows",
    "Kobresia spp. and forb high-cold meadows",
    "Alpine tundra",
    "Subalpine broadleaf deciduous scrubs",
    "No vegetation",
    "Alpine sparse vegetation",
    "One year two ripes or three ripes grain fields; evergreen orchards and subtropical economic tree plantations",
    "Alpine cushion vegetation",
    "Shrub deserts",
    "Semi-shrub and dwarf semi-shrub deserts",
    "Temperate tufted low grass and nano-semi-shrub desert steppes",
    "Bamboo forests and scrubs in subtropical and tropical zones",
    "Two years three ripes or one year two ripes grain fields and deciduous orchards",
    "Broadleaf deciduous forests in subtropical zone",
    "Dwarf semi-arboreous deserts",
    "Needleleaf forests in subtropical zone",
    "Grass and Carex spp. high-cold steppes",
    "Succulent halophytic dwarf semi-shrub deserts",
    "One year one ripe grain fields and cold-tolerant economic crop fields; deciduous orchards",
    "Needleleaf, evergreen and deciduous broadleaf mixed forests on mountains in subtropical zone",
    "Steppe shrub deserts",
    "Broadleaf evergreen and deciduous mixed forests in subtropical zone",
    "Needleleaf forests in temperate zone",
    "Annual herb deserts",
    "Sclerophyllous broadleaf evergreen forests in subtropical zone",
    "Subtropical and tropical grasslands",
    "Temperate grasslands",
    "Cushion nano-semi-shrub high-cold deserts",
    "High-cold marshes",
    "Broadleaf evergreen and deciduous scrubs in subtropical and tropical zones",
    "Subalpine needleleaf evergreen scrubs",
    "Subalpine sclerophyllous broadleaf evergreen scrubs",
    "One year two ripes grain fields; evergreen and deciduous orchards; economic tree plantations",
    "Subtropical and tropical marshes",
    "Broadleaf evergreen forests in subtropical zone",
    "Evergreen xeromorphic succulent thorny scrubs in subtropical and tropical zones",
    "Needleleaf forests in tropical zone",
    "Tropical rain forests",
    "One year three ripes grain fields; tropical evergreen orchards and economic tree plantations",
    "Tropical mangroves",
    "Tropical monsoon forests",
    "Broadleaf evergreen succulent scrub and dwarf forest on coral islands in tropical zone",
    "Novel zone / no suitable current zone"
  ),
  category = c(
    "cropland", "cold_temperate_needleleaf_forest", "temperate_deciduous_forest", "swamp_meadow", "meadow", "marsh",
    "temperate_scrub", "rare_or_error", "mixed_forest", "mountain_needleleaf_forest", "meadow_steppe", "cropland",
    "tufted_grass_steppe", "temperate_woodland", "halophytic_meadow", "high_cold_meadow", "alpine_tundra", "subalpine_scrub",
    "no_vegetation", "alpine_sparse", "cropland_orchard_plantation", "alpine_cushion", "shrub_desert", "semi_shrub_desert",
    "desert_steppe", "bamboo", "cropland_orchard", "subtropical_deciduous_forest", "semi_arboreous_desert", "subtropical_needleleaf_forest",
    "high_cold_steppe", "halophytic_desert", "cropland_orchard", "subtropical_mixed_forest", "steppe_shrub_desert", "subtropical_mixed_broadleaf_forest",
    "temperate_needleleaf_forest", "annual_herb_desert", "sclerophyll_broadleaf_forest", "subtropical_tropical_grassland", "temperate_grassland", "high_cold_desert",
    "high_cold_marsh", "subtropical_tropical_scrub", "subalpine_needleleaf_scrub", "subalpine_sclerophyll_scrub", "cropland_orchard_plantation", "subtropical_tropical_marsh",
    "subtropical_evergreen_forest", "xeromorphic_succulent_thorny_scrub", "tropical_needleleaf_forest", "tropical_rain_forest", "tropical_cropland_orchard_plantation", "mangrove",
    "tropical_monsoon_forest", "coral_island_scrub_dwarf_forest", "novel"
  ),
  COLOR = c(
    "#FFD23F", "#005A32", "#33A02C", "#00B4D8", "#B8DE29", "#0077B6",
    "#C77DFF", "#BDBDBD", "#2CA25F", "#006D77", "#D9ED92", "#F8961E",
    "#DDA15E", "#8AC926", "#FF70A6", "#90BE6D", "#A2D2FF", "#9D4EDD",
    "#E0E0E0", "#B8C0FF", "#F3722C", "#CDB4DB", "#B5651D", "#9C6644",
    "#BC6C25", "#00A878", "#FFB703", "#52B788", "#7F5539", "#1B9AAA",
    "#48CAE4", "#E63946", "#FDB833", "#2A9D8F", "#D62828", "#40916C",
    "#0B3D2E", "#F4A261", "#118AB2", "#AACC00", "#70E000", "#8D99AE",
    "#4CC9F0", "#F72585", "#7209B7", "#B5179E", "#FF6D00", "#00BFA6",
    "#006400", "#FF006E", "#003F5C", "#00C853", "#FF9F1C", "#006D5B",
    "#00A676", "#FF477E", "#1A1A1A"
  ),
  count = c(
    15353, 284343, 553420, 92312, 191879, 92003,
    221602, NA_integer_, 30540, 147373, 147426, 375765,
    620911, 33001, 242556, 917167, 2032, 113046,
    800196, 411347, 610148, 43164, 330737, 846694,
    237969, 40966, 866811, 46263, 164699, 580358,
    848105, 62847, 391144, 5163, 32917, 22763,
    193911, 30035, 22033, 341246, 91481, 157267,
    5708, 576554, 16418, 136995, 359409, 2400,
    158047, 2515, 40, 24103, 67859, 11883,
    6461, 3, NA_integer_
  ),
  stringsAsFactors = FALSE
)

# Save files used by later mapping scripts.
writeLines(
  toJSON(color_df, dataframe = "rows", pretty = TRUE, auto_unbox = TRUE, na = "null"),
  "color_palette_China.json"
)
write.csv(color_df, "color_palette_China.csv", row.names = FALSE, na = "")

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

read_zone_palette <- function(json_file = "color_palette_China.json") {
  pal <- fromJSON(json_file)
  pal$zoneID <- as.integer(pal$zoneID)
  pal$zone <- as.character(pal$zone)
  pal$category <- as.character(pal$category)
  pal$COLOR <- as.character(pal$COLOR)
  
  # This prevents the old "1. NA" name-legend problem if an older JSON file
  # is accidentally read. The hard-coded names here are the authority.
  zone_lut <- color_df[, c("zoneID", "zone", "category", "count")]
  names(zone_lut) <- c("zoneID", "zone_hard", "category_hard", "count_hard")
  
  pal <- merge(pal[, c("zoneID", "COLOR")], zone_lut,
               by = "zoneID", all.x = TRUE, sort = FALSE)
  names(pal)[names(pal) == "zone_hard"] <- "zone"
  names(pal)[names(pal) == "category_hard"] <- "category"
  names(pal)[names(pal) == "count_hard"] <- "count"
  
  pal[order(pal$zoneID), c("zoneID", "zone", "category", "COLOR", "count")]
}

get_zoneID <- function(r = NULL, zoneID = NULL) {
  if (!is.null(zoneID)) {
    return(sort(unique(as.integer(zoneID))))
  }
  
  if (!is.null(r)) {
    z <- values(r, mat = FALSE)
    z <- z[!is.na(z)]
    return(sort(unique(as.integer(z))))
  }
  
  NULL
}

filter_palette <- function(pal, r = NULL, zoneID = NULL) {
  keep <- get_zoneID(r = r, zoneID = zoneID)
  
  if (is.null(keep)) {
    return(pal[order(pal$zoneID), ])
  }
  
  missing <- setdiff(keep, pal$zoneID)
  if (length(missing) > 0) {
    stop(
      "Some raster values do not have colors in the palette JSON: ",
      paste(missing, collapse = ", ")
    )
  }
  
  out <- pal[pal$zoneID %in% keep, ]
  out[order(out$zoneID), ]
}

open_device <- function(file, width, height, res = 300) {
  ext <- tolower(tools::file_ext(file))
  
  if (ext %in% c("tif", "tiff")) {
    tiff(file, width = width, height = height, units = "in",
         res = res, compression = "lzw")
  } else if (ext == "png") {
    png(file, width = width, height = height, units = "in", res = res)
  } else if (ext %in% c("jpg", "jpeg")) {
    jpeg(file, width = width, height = height, units = "in", res = res)
  } else {
    stop("Unsupported output format: ", ext)
  }
}

# ------------------------------------------------------------
# 3. Plot map
# ------------------------------------------------------------

plot_zone_raster <- function(r, pal, out_file, main = "Vegetation zones") {
  zoneID <- get_zoneID(r = r)
  pal_now <- filter_palette(pal, zoneID = zoneID)
  
  # terra::plot() uses sequential color positions, so recode only for plotting.
  r_plot <- subst(r, from = zoneID, to = seq_along(zoneID), others = NA)
  
  open_device(out_file, width = 10, height = 8, res = 300)
  plot(r_plot, col = pal_now$COLOR, main = main, axes = TRUE, legend = FALSE)
  dev.off()
}

# ------------------------------------------------------------
# 4. Plot legends
# ------------------------------------------------------------

plot_zone_legend <- function(pal, out_file, type = c("ID", "name"),
                             r = NULL, zoneID = NULL, drop_zoneID = 8) {
  type <- match.arg(type)
  
  legend_df <- filter_palette(pal, r = r, zoneID = zoneID)
  legend_df <- legend_df[!legend_df$zoneID %in% drop_zoneID, ]
  legend_df <- legend_df[order(legend_df$zoneID), ]
  
  if (nrow(legend_df) == 0) {
    stop("No zones left to plot after filtering the legend.")
  }
  
  if (any(is.na(legend_df$zone))) {
    stop("Some zone names are NA: ",
         paste(legend_df$zoneID[is.na(legend_df$zone)], collapse = ", "))
  }
  
  if (type == "ID") {
    n <- nrow(legend_df)
    n_col <- 2
    n_row <- ceiling(n / n_col)
    
    col_id <- rep(seq_len(n_col), length.out = n)
    row_id <- rep(n_row:1, each = n_col)[seq_len(n)]
    
    # The ID legend is compact. The color bar is intentionally narrow.
    col_start <- c(0.10, 0.62)
    bar_w <- 0.018
    text_gap <- 0.025
    
    x0 <- col_start[col_id]
    x1 <- x0 + bar_w
    tx <- x1 + text_gap
    y <- row_id
    
    open_device(out_file, width = 7.0, height = max(2.2, 0.25 * n_row + 0.8), res = 300)
    par(bg = "white", fg = "black", col.axis = "black", col.main = "black",
        mar = c(0.2, 0.2, 1.0, 0.2), xpd = NA)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0.5, n_row + 0.5),
                xaxs = "i", yaxs = "i")
    title(main = "Legend: Vegetation Type ID", adj = 0, cex.main = 1.0)
    
    rect(xleft = x0, ybottom = y - 0.32,
         xright = x1, ytop = y + 0.32,
         col = legend_df$COLOR, border = NA)
    text(tx, y, labels = legend_df$zoneID,
         adj = c(0, 0.5), cex = 0.85, col = "black")
    dev.off()
    
  } else {
    n <- nrow(legend_df)
    n_col <- 2
    n_row <- ceiling(n / n_col)
    
    col_id <- rep(seq_len(n_col), length.out = n)
    row_id <- rep(n_row:1, each = n_col)[seq_len(n)]
    
    # Two columns keep the name legend readable without forcing a huge image.
    col_start <- c(0.02, 0.51)
    bar_w <- 0.012
    text_gap <- 0.012
    
    x0 <- col_start[col_id]
    x1 <- x0 + bar_w
    tx <- x1 + text_gap
    y <- row_id
    
    label <- paste0(legend_df$zoneID, ". ", legend_df$zone)
    
    open_device(out_file, width = 20, height = max(6, 0.34 * n_row + 0.8), res = 300)
    par(bg = "white", fg = "black", col.axis = "black", col.main = "black",
        mar = c(0.2, 0.2, 1.0, 0.2), xpd = NA)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0.5, n_row + 0.5),
                xaxs = "i", yaxs = "i")
    title(main = "Legend: Vegetation Type Names", adj = 0, cex.main = 1.0)
    
    rect(xleft = x0, ybottom = y - 0.30,
         xright = x1, ytop = y + 0.30,
         col = legend_df$COLOR, border = NA)
    text(tx, y, labels = label,
         adj = c(0, 0.5), cex = 0.60, col = "black")
    dev.off()
  }
}

# ------------------------------------------------------------
# 5. Example use
# ------------------------------------------------------------

pal <- read_zone_palette()

dir.create("result maps", showWarnings = FALSE, recursive = TRUE)

r <- rast("raster/veg_3")
plot_zone_raster(
  r, pal,
  "result maps/veg_3_color.tif",
  "Vegetation Types of China"
)

plot_zone_legend(pal, "result maps/legend_ID.png", type = "ID", r = r)
plot_zone_legend(pal, "result maps/legend_names.png", type = "name", r = r)

# Future predicted map example:
# pred <- rast("path/to/predicted_zone_map.tif")
# plot_zone_raster(pred, pal,
#                  "result maps/pred_zone_color.tif",
#                  "Predicted Vegetation Zones")
# plot_zone_legend(pal, "result maps/pred_legend_ID.png", type = "ID", r = pred)
# plot_zone_legend(pal, "result maps/pred_legend_names.png", type = "name", r = pred)

# ------------------------------------------------------------
# Plot original map and predicted normal-period map
# ------------------------------------------------------------

dir.create("result maps", showWarnings = FALSE, recursive = TRUE)

# Original vegetation raster
r_ori <- rast("raster/veg_3")

plot_zone_raster(
  r_ori, pal,
  "result maps/veg_3_original_color.tif",
  "Original Vegetation Types of China"
)

plot_zone_legend(
  pal,
  "result maps/legend_original_ID.png",
  type = "ID",
  r = r_ori
)

plot_zone_legend(
  pal,
  "result maps/legend_original_names.png",
  type = "name",
  r = r_ori
)


# Predicted normal-period vegetation raster
r_pred <- rast("result maps/assigned_zone_normal_threshold0.1_novel99_originalTie.tif")

plot_zone_raster(
  r_pred, pal,
  "result maps/pred_normal_color.tif",
  "Predicted Vegetation Zones: Normal Period"
)

