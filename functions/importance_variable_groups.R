# Shared feature-importance annotation; no model fitting or raster operations.
# Temperature statistic and climate period are separate dimensions.
importance_grouping_version <- "climate_statistic_period_v2"
importance_climate_periods <- c("Annual", "Winter", "Spring", "Summer", "Autumn")
importance_climate_groups <- c(
  "Mean temperature", "Mean maximum temperature", "Mean minimum temperature",
  "Extreme maximum temperature", "Extreme minimum temperature",
  "Warmest-month mean temperature", "Coldest-month mean temperature",
  "Annual temperature range", "Growing season and degree-days",
  "Radiation and evaporative demand", "Precipitation", "Climatic moisture balance",
  "Humidity", "Snowfall"
)
importance_soil_groups <- c(
  "Nutrient retention and base status", "Bulk density", "Coarse fragments",
  "Texture", "Organic carbon", "Soil reaction", "Carbonates", "Gypsum",
  "Sodicity", "Salinity", "Other"
)
importance_semantic_colours <- c(
  "Mean temperature" = "#E69F00",
  "Mean maximum temperature" = "#D55E00",
  "Mean minimum temperature" = "#F6C667",
  "Extreme maximum temperature" = "#9E2A2B",
  "Extreme minimum temperature" = "#A66B56",
  "Warmest-month mean temperature" = "#F46D43",
  "Coldest-month mean temperature" = "#D9A17B",
  "Annual temperature range" = "#B8860B",
  "Precipitation" = "#2C7FB8",
  "Climatic moisture balance" = "#1D91A0",
  "Growing season and degree-days" = "#A6761D",
  "Radiation and evaporative demand" = "#9B8A3C",
  "Humidity" = "#67A9CF", "Snowfall" = "#A6DCEB",
  "Nutrient retention and base status" = "#2E8B57",
  "Bulk density" = "#6B4C3B", "Coarse fragments" = "#7A7F87",
  "Texture" = "#C17C2C", "Organic carbon" = "#4E342E",
  "Soil reaction" = "#7B3294", "Carbonates" = "#C9A227",
  "Gypsum" = "#CC79A7", "Sodicity" = "#E76F51",
  "Salinity" = "#008C95", "Other" = "#BDBDBD"
)

annotate_importance_variables <- function(niche, variable, strict = TRUE) {
  niche <- as.character(niche)
  variable <- as.character(variable)
  if (length(niche) != length(variable) || anyNA(niche) || anyNA(variable)) {
    stop("niche and variable must have equal lengths and no missing names.")
  }
  if (any(!niche %in% c("climate", "soil"))) stop("Unknown niche component.")
  # Preserve DD18 versus DD_18. Only case and hyphens are standardized.
  standard <- toupper(gsub("-", "_", trimws(variable), fixed = TRUE))
  seasonal <- grepl("_(WT|SP|SM|AT)$", standard)
  base <- sub("_(WT|SP|SM|AT)$", "", standard)
  climate <- niche == "climate"
  soil <- niche == "soil"
  group <- rep("Other", length(variable))
  matches <- integer(length(variable))
  climate_rules <- list(
    "Mean temperature" = "^(MAT|TAVE)$",
    "Mean maximum temperature" = "^TMAX$",
    "Mean minimum temperature" = "^TMIN$",
    "Extreme maximum temperature" = "^EXT$",
    "Extreme minimum temperature" = "^EMT$",
    "Warmest-month mean temperature" = "^MWMT$",
    "Coldest-month mean temperature" = "^MCMT$",
    "Annual temperature range" = "^TD$",
    "Precipitation" = "^(MAP|MSP|PPT)$",
    "Climatic moisture balance" = "^(AHM|SHM|CMD|CMI)$",
    "Growing season and degree-days" = "^(BFFP|EFFP|FFP|NFFD|DD_?[0-9]+)$",
    "Radiation and evaporative demand" = "^(EREF|RSDS)$",
    "Humidity" = "^RH$", "Snowfall" = "^PAS$"
  )
  for (label in names(climate_rules)) {
    hit <- climate & grepl(climate_rules[[label]], base)
    group[hit] <- label
    matches[hit] <- matches[hit] + 1L
  }
  if (any(climate & matches > 1L)) stop("A climate variable matched multiple groups.")
  if (strict && any(climate & matches != 1L)) {
    stop("Unclassified climate variable(s): ",
         paste(unique(variable[climate & matches != 1L]), collapse = ", "),
         ". Add an explicit rule before interpreting grouped results.")
  }

  # Keep the original soil rules and their precedence unchanged.
  group[soil & grepl("CEC|TEB|BASE|(^|_)BS($|_)", standard)] <-
    "Nutrient retention and base status"
  group[soil & grepl("BLD|BULK|DENS", standard)] <- "Bulk density"
  group[soil & grepl("CRF|GRAVEL|COARSE", standard)] <- "Coarse fragments"
  group[soil & group == "Other" & grepl("SND|SLT|CLY|SAND|SILT|CLAY|TEXT", standard)] <- "Texture"
  group[soil & grepl("(^|_)OC($|_)|ORC|ORGANIC", standard)] <- "Organic carbon"
  group[soil & grepl("PH", standard)] <- "Soil reaction"
  group[soil & grepl("TCEQ|CACO3|CARBONATE", standard)] <- "Carbonates"
  group[soil & grepl("GY|CASO4", standard)] <- "Gypsum"
  group[soil & grepl("ESP|SOD", standard)] <- "Sodicity"
  group[soil & grepl("(^|_)ECE($|_)|SALIN|ELECTRICAL", standard)] <- "Salinity"

  period <- rep(NA_character_, length(variable))
  period[climate] <- "Annual"
  suffix <- sub("^.*_", "", standard)
  period_labels <- c(WT = "Winter", SP = "Spring", SM = "Summer", AT = "Autumn")
  for (suffix_value in names(period_labels)) {
    period[climate & seasonal & suffix == suffix_value] <- period_labels[[suffix_value]]
  }
  detail <- rep(NA_character_, length(variable))
  detail[climate] <- "Annual-scale variable or index; not necessarily an annual mean"
  detail[climate & base == "MAT" & !seasonal] <- "Mean annual temperature"
  detail[climate & base == "MWMT"] <- "Mean temperature of the warmest month"
  detail[climate & base == "MCMT"] <- "Mean temperature of the coldest month"
  detail[climate & base == "TD"] <- "Warmest-month minus coldest-month mean temperature"
  detail[climate & base == "MSP"] <- "Growing-season precipitation index"
  detail[climate & seasonal] <- paste(period[climate & seasonal], "seasonal variable")
  data.table::data.table(
    niche = niche, variable = variable, variable_group = group,
    climate_period = period, climate_statistic = ifelse(climate, group, NA_character_),
    climate_period_detail = detail
  )
}

reannotate_importance_table <- function(x) {
  x <- data.table::copy(data.table::as.data.table(x))
  if (!all(c("niche", "variable") %in% names(x))) {
    stop("Importance table must contain niche and variable.")
  }
  old <- intersect(c("variable_group", "climate_period", "climate_statistic",
                     "climate_period_detail"), names(x))
  if (length(old)) x[, (old) := NULL]
  annotation <- unique(annotate_importance_variables(x$niche, x$variable))
  answer <- merge(x, annotation, by = c("niche", "variable"), all.x = TRUE, sort = FALSE)
  if (nrow(answer) != nrow(x)) stop("Annotation unexpectedly changed the row count.")
  answer
}

write_climate_group_period_tables <- function(aligned, table_dir) {
  keep <- function(x) {
    tolower(gsub("[ -]+", "_", trimws(as.character(x)))) != "no_vegetation"
  }
  climate <- data.table::copy(aligned[niche == "climate"])
  n_climate_zones <- data.table::uniqueN(climate$zoneID)
  coverage <- climate[, .(n_zones = data.table::uniqueN(zoneID)), by = variable]
  if (!nrow(climate) || any(coverage$n_zones != n_climate_zones)) {
    stop("The aligned climate profile is incomplete: retain zero rows for unselected variables.")
  }
  zone_table <- climate[, .(
    n_variables = .N, n_binary_selected = sum(binary_selected),
    n_multiclass_selected = sum(multiclass_selected),
    binary_share = sum(binary_share), multiclass_share = sum(multiclass_share),
    consensus_share = sum(consensus_share)
  ), by = .(zoneID, zone, category, category2, variable_group, climate_period)]
  # Complete observed group-period cells across zones; unselected predictors
  # already have zero shares in the aligned profile. Structural empty cells
  # (e.g. seasonal MWMT) are deliberately absent, not fabricated predictors.
  totals <- zone_table[, .(group_period_sum = sum(consensus_share)), by = zoneID]
  expected <- climate[, .(profile_sum = sum(consensus_share)), by = zoneID]
  check <- merge(totals, expected, by = "zoneID")
  if (any(abs(check$group_period_sum - check$profile_sum) > 1e-10)) {
    stop("Group-period aggregation changed a climatic importance total.")
  }
  se <- function(x) if (length(x) < 2L) NA_real_ else stats::sd(x) / sqrt(length(x))
  category_table <- zone_table[, .(
    n_zones = .N, binary_mean = mean(binary_share), binary_se = se(binary_share),
    multiclass_mean = mean(multiclass_share), multiclass_se = se(multiclass_share),
    consensus_mean = mean(consensus_share), consensus_se = se(consensus_share)
  ), by = .(category2, variable_group, climate_period)]
  data.table::setorder(zone_table, category2, zoneID, variable_group, climate_period)
  data.table::setorder(category_table, category2, variable_group, climate_period)
  tables <- list(
    FI_12_zone_climate_group_period = zone_table,
    FI_13_category_climate_group_period = category_table,
    FI_reader_10_zone_climate_group_period = zone_table[keep(category2)],
    FI_reader_11_category_climate_group_period = category_table[keep(category2)]
  )
  for (name in names(tables)) {
    data.table::fwrite(tables[[name]], file.path(table_dir, paste0(name, ".csv")))
  }
  invisible(tables)
}
