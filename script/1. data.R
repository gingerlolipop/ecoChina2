# Prepare the reference ecotype raster and reference-period climate table.
# Existing legacy outputs are reused by default because raster extraction is slow.

library(data.table)
library(terra)

rm(list = ls())
gc()

find_project_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "script", "1. data.R"))) return(path)
    parent <- dirname(path)
    if (parent == path) stop("Run inside the repository or set ECOCHINA2_DIR.")
    path <- parent
  }
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "true" else "false")
  tolower(trimws(value)) %in% c("1", "true", "yes", "y")
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[1] else NA_character_
}

env_root <- Sys.getenv("ECOCHINA2_DIR", unset = "")
base_dir <- if (nzchar(env_root)) {
  normalizePath(env_root, winslash = "/", mustWork = TRUE)
} else {
  find_project_root()
}

raw_dir <- file.path(base_dir, "data raw")
raster_dir <- file.path(base_dir, "raster")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)

force_data <- env_flag("ECOCHINA2_FORCE_DATA", FALSE)

reference_file <- file.path(raster_dir, "ecosys_ori.tif")
zone_lookup_file <- file.path(raw_dir, "1. zoneID_zone_count.csv")
dem_file <- file.path(raw_dir, "1. zoneID_dem.csv")
climate_file <- file.path(
  raw_dir,
  "1. zoneID_Clm_800m_Normal_1961_1990SY.csv"
)

# 1. Reference raster ---------------------------------------------------------

if (!force_data && file.exists(reference_file)) {
  cat("[USE EXISTING]", reference_file, "\n")
  r <- rast(reference_file)
} else {
  veg_file <- first_existing(c(
    file.path(raster_dir, "veg_3"),
    file.path(raster_dir, "veg_3.tif"),
    file.path(base_dir, "data", "raw", "veg_3.tif")
  ))
  if (is.na(veg_file)) stop("Missing reference vegetation raster: raster/veg_3")

  r <- rast(veg_file)
  if (!nzchar(crs(r))) crs(r) <- "EPSG:4326"
  writeRaster(r, reference_file, filetype = "GTiff", overwrite = TRUE)
  cat("[SAVED]", reference_file, "\n")
}

if (!force_data && file.exists(zone_lookup_file)) {
  cat("[USE EXISTING]", zone_lookup_file, "\n")
} else {
  zone_lookup <- tryCatch(as.data.frame(cats(r)[[1]]), error = function(e) NULL)
  if (is.null(zone_lookup) || !nrow(zone_lookup)) {
    zone_lookup <- as.data.frame(freq(r))
  }
  fwrite(zone_lookup, zone_lookup_file)
  cat("[SAVED]", zone_lookup_file, "\n")
}

# 2. Coordinates and DEM ------------------------------------------------------

if (!force_data && file.exists(dem_file)) {
  cat("[USE EXISTING]", dem_file, "\n")
} else {
  coord_file <- first_existing(c(
    file.path(raw_dir, "1. coord.csv"),
    file.path(base_dir, "data", "raw", "coord.csv")
  ))
  if (is.na(coord_file)) stop("Missing DEM coordinates: data raw/1. coord.csv")

  id <- as.data.table(as.data.frame(
    r,
    cells = TRUE,
    xy = TRUE,
    na.rm = TRUE
  ))
  setnames(id, names(id)[4], "zoneID")
  id <- id[complete.cases(id[, .(cell, x, y, zoneID)])]
  setcolorder(id, c("cell", "zoneID", "y", "x"))

  coord <- fread(coord_file)
  dem_col <- intersect(c("china_90m", "elevation", "dem"), names(coord))[1]
  if (is.na(dem_col)) stop("No elevation column found in: ", coord_file)
  if (nrow(coord) != nrow(id)) {
    stop("DEM row count does not match the non-NA reference raster cells.")
  }

  id[, china_90m := coord[[dem_col]]]
  fwrite(id, dem_file)
  cat("[SAVED]", dem_file, "\n")
  rm(id, coord)
  gc()
}

# 3. Reference-period climate -------------------------------------------------

if (!force_data && file.exists(climate_file)) {
  cat("[USE EXISTING]", climate_file, "\n")
} else {
  # Accept the short-lived clean-layout table as a fallback, but always restore
  # the established filename used by the analysis.
  fallback_table <- file.path(base_dir, "data", "climate_reference.csv")
  if (!force_data && file.exists(fallback_table)) {
    fwrite(fread(fallback_table), climate_file)
    cat("[RESTORED LEGACY NAME]", climate_file, "\n")
  } else {
    env_climate <- Sys.getenv("ECOCHINA2_CLIMATE_DIR", unset = "")
    climate_root <- if (nzchar(env_climate)) {
      env_climate
    } else {
      file.path(base_dir, "data", "rasters", "climate")
    }
    climate_dir <- if (basename(climate_root) == "Normal_1961_1990") {
      climate_root
    } else {
      file.path(climate_root, "Normal_1961_1990")
    }
    if (!dir.exists(climate_dir)) {
      stop("Missing reference climate raster directory. Set ECOCHINA2_CLIMATE_DIR.")
    }

    dem <- fread(dem_file)
    climate_files <- list.files(climate_dir, pattern = "\\.tif$", full.names = TRUE)
    if (!length(climate_files)) stop("No climate .tif files found in: ", climate_dir)

    climate_stack <- rast(climate_files)
    names(climate_stack) <- tools::file_path_sans_ext(basename(climate_files))
    climate_values <- as.data.table(terra::extract(
      climate_stack,
      dem[, .(x, y)]
    ))
    climate_values[, ID := NULL]
    fwrite(cbind(dem, climate_values), climate_file)
    cat("[SAVED]", climate_file, "\n")
  }
}

cat(
  "\nCOMPLETE\n",
  "Reference raster: ", reference_file, "\n",
  "Reference climate table: ", climate_file, "\n",
  "Set ECOCHINA2_FORCE_DATA=true only to rebuild existing outputs.\n",
  sep = ""
)
