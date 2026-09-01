# Define populations as species x reference-period ecotype combinations.
# Existing population tables are reused before the reference raster is opened.

library(data.table)
library(terra)

rm(list = ls())
gc()

find_project_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "script", "1.2 species to pop2.R"))) return(path)
    parent <- dirname(path)
    if (parent == path) stop("Run inside the repository or set ECOCHINA2_DIR.")
    path <- parent
  }
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "true" else "false")
  tolower(trimws(value)) %in% c("1", "true", "yes", "y")
}

env_root <- Sys.getenv("ECOCHINA2_DIR", unset = "")
base_dir <- if (nzchar(env_root)) {
  normalizePath(env_root, winslash = "/", mustWork = TRUE)
} else {
  find_project_root()
}

env_species <- Sys.getenv("ECOCHINA2_SPECIES_DIR", unset = "")
species_dir <- if (nzchar(env_species)) {
  env_species
} else {
  file.path(base_dir, "species data")
}

reference_file <- file.path(base_dir, "raster", "ecosys_ori.tif")
force_populations <- env_flag("ECOCHINA2_FORCE_POPULATIONS", FALSE)
min_population_cells <- 10L
zoneID <- 1:55

output_files <- c(
  file.path(base_dir, "specie-zone-pop.csv"),
  file.path(base_dir, "species_population_summary.csv"),
  file.path(base_dir, "species_zone_population_long.csv"),
  file.path(base_dir, "sorted_population_zone.csv")
)

if (!force_populations && all(file.exists(output_files))) {
  cat(
    "[USE EXISTING] Population tables are complete.\n",
    paste(output_files, collapse = "\n"), "\n",
    "Set ECOCHINA2_FORCE_POPULATIONS=true to rebuild them.\n",
    sep = ""
  )
} else {
  if (!file.exists(reference_file)) stop("Run script/1. data.R first.")

  species_info <- data.table(
    SpeciesCode = c(
      "querVar", "pinuYun", "lariGme", "phylPub", "lariOlg",
      "robiPse", "chamFor", "cyclLon", "pinuSyl", "saliMat"
    ),
    LatinName = c(
      "Quercus variabilis", "Pinus yunnanensis", "Larix gmelinii",
      "Phyllostachys pubescens", "Larix olgensis var. changpaiensis",
      "Robinia pseudoacacia", "Chamaecyparis formosensis",
      "Cyclobalanopsis longinux", "Pinus sylvestris var. mongolica",
      "Salix matsudana"
    ),
    CommonName = c(
      "Chinese Cork Oak", "Yunnan Pine", "Dahurian Larch", "Moso Bamboo",
      "Hinggan Larch", "Black Locust", "Taiwan Cypress", "Longinux Oak",
      "Mongolian Pine", "Chinese Willow"
    )
  )
  species_names <- c(
    "lariGme", "lariOlg", "pinuSyl", "querVar", "robiPse",
    "saliMat", "phylPub", "pinuYun", "chamFor", "cyclLon"
  )

  reference <- rast(reference_file)
  names(reference) <- "zoneID"

  presence_rows <- function(x) {
    presence_col <- intersect(
      c("presence", "Presence", "present", "Present", "status", "Status"),
      names(x)
    )[1]
    if (is.na(presence_col)) presence_col <- names(x)[2]
    value <- tolower(trimws(as.character(x[[presence_col]])))
    x[value %in% c("y", "yes", "1", "true", "presence", "present")]
  }

  count_zone_cells <- function(species_data) {
    lon_col <- intersect(c("lon", "longitude", "x"), names(species_data))[1]
    lat_col <- intersect(c("lat", "latitude", "y"), names(species_data))[1]
    if (is.na(lon_col) || is.na(lat_col)) {
      stop("Species coordinates must include lon/lat columns.")
    }

    points <- vect(
      species_data,
      geom = c(lon_col, lat_col),
      crs = "EPSG:4326"
    )
    if (!same.crs(points, reference)) points <- project(points, crs(reference))

    occupied <- unique(data.table(
      cell = cellFromXY(reference, crds(points)),
      zone = as.integer(terra::extract(reference, points)$zoneID)
    ))
    occupied <- occupied[!is.na(cell) & zone %in% zoneID]
    occupied[, .N, by = zone]
  }

  results <- matrix(
    NA_integer_,
    nrow = length(zoneID),
    ncol = length(species_names),
    dimnames = list(paste0("Zone", zoneID), species_names)
  )

  for (species in species_names) {
    species_file <- file.path(species_dir, paste0(species, "_coord.csv"))
    if (!file.exists(species_file)) stop("Missing species file: ", species_file)

    species_data <- presence_rows(fread(species_file))
    if (!nrow(species_data)) {
      cat("[NO PRESENCE]", species, "\n")
      next
    }

    counts <- count_zone_cells(species_data)
    results[match(paste0("Zone", counts$zone), rownames(results)), species] <- counts$N
    cat("[OVERLAY]", species, "\n")
  }

  results <- results[
    apply(results, 1, function(row) any(!is.na(row))),
    ,
    drop = FALSE
  ]
  write.csv(results, output_files[1], row.names = TRUE)

  population_abundance <- results
  population_abundance[
    !is.na(population_abundance) &
      population_abundance < min_population_cells
  ] <- NA
  population_abundance <- population_abundance[
    apply(population_abundance, 1, function(row) any(!is.na(row))),
    ,
    drop = FALSE
  ]

  shannon_index <- function(x) {
    x <- x[!is.na(x) & is.finite(x) & x > 0]
    if (!length(x)) return(NA_real_)
    p <- x / sum(x)
    -sum(p * log(p))
  }

  species_summary <- data.table(
    SpeciesCode = species_names,
    NumberPopulations = colSums(!is.na(population_abundance)),
    H_species = round(apply(population_abundance, 2, shannon_index), 4)
  )
  species_summary <- merge(
    species_info,
    species_summary,
    by = "SpeciesCode",
    all.y = TRUE,
    sort = FALSE
  )
  setorder(species_summary, -NumberPopulations, SpeciesCode)
  fwrite(species_summary, output_files[2])

  long_rows <- vector("list", length(species_names))
  for (j in seq_along(species_names)) {
    keep <- which(!is.na(population_abundance[, species_names[j]]))
    long_rows[[j]] <- data.table(
      PopulationID = paste(species_names[j], rownames(population_abundance)[keep], sep = "-"),
      Species = species_names[j],
      Zone = rownames(population_abundance)[keep],
      Population = as.numeric(population_abundance[keep, species_names[j]])
    )
  }
  population_long <- rbindlist(long_rows, use.names = TRUE)
  fwrite(population_long, output_files[3])

  population_wide <- dcast(
    population_long,
    PopulationID ~ Zone,
    value.var = "Population",
    fill = 0
  )
  zone_columns <- grep("^Zone[0-9]+$", names(population_wide), value = TRUE)
  zone_columns <- zone_columns[order(as.integer(sub("Zone", "", zone_columns)))]
  population_wide[, SpeciesOrder := sub("-Zone.*", "", PopulationID)]
  population_wide[, ZoneOrder := as.integer(sub(".*-Zone", "", PopulationID))]
  setorder(population_wide, SpeciesOrder, ZoneOrder)
  population_wide[, c("SpeciesOrder", "ZoneOrder") := NULL]
  setcolorder(population_wide, c("PopulationID", zone_columns))
  fwrite(population_wide, output_files[4])

  cat(
    "\nCOMPLETE\n",
    "Populations retained: ", nrow(population_long), "\n",
    "Outputs retain their established repository-root filenames.\n",
    sep = ""
  )
}
