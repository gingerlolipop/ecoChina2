# Small tests; no fitted models or raster files are needed.
library(data.table)
configured <- Sys.getenv("ECOCHINA2_DIR", unset = "")
root <- if (nzchar(configured)) configured else getwd()
source(file.path(root, "functions", "importance_variable_groups.R"), local = TRUE)

names_to_check <- c("MAT", "Tave_sm", "TMAX_WT", "Tmin-at", "EXT", "EMT",
                    "MWMT", "MCMT", "TD", "DD18", "DD_18", "DD_0_wt", "PPT_sp")
expected <- c("Mean temperature", "Mean temperature", "Mean maximum temperature",
              "Mean minimum temperature", "Extreme maximum temperature",
              "Extreme minimum temperature", "Warmest-month mean temperature",
              "Coldest-month mean temperature", "Annual temperature range",
              "Growing season and degree-days", "Growing season and degree-days",
              "Growing season and degree-days", "Precipitation")
a <- annotate_importance_variables(rep("climate", length(names_to_check)), names_to_check)
stopifnot(identical(a$variable_group, expected),
          a$climate_period[2] == "Summer", a$climate_period[3] == "Winter",
          a$climate_period[4] == "Autumn", a$climate_period[13] == "Spring",
          identical(a$variable, names_to_check),
          a$climate_period[7] == "Annual")
soil <- annotate_importance_variables(rep("soil", 5),
                                      c("T_CEC_CLAY", "T_CLY", "T_BS", "T_OC", "T_PH_H2O"))
stopifnot(identical(soil$variable_group,
  c("Nutrient retention and base status", "Texture", "Nutrient retention and base status",
    "Organic carbon", "Soil reaction")))
unknown_stopped <- tryCatch({
  annotate_importance_variables("climate", "temperature_unknown")
  FALSE
}, error = function(e) TRUE)
stopifnot(unknown_stopped)

fixture <- data.table(
  zoneID = rep(c(1L, 2L), each = 3), zone = rep(c("A", "B"), each = 3),
  category = "Forest", category2 = "forest", niche = "climate",
  variable = rep(c("Tave_sm", "Tmin_wt", "MWMT"), 2),
  binary_share = c(.6, .3, .1, 0, .5, .5),
  multiclass_share = c(.4, .4, .2, .2, .4, .4),
  binary_selected = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
  multiclass_selected = TRUE
)
fixture[, consensus_share := (binary_share + multiclass_share) / 2]
x <- reannotate_importance_table(fixture)
stopifnot(nrow(x) == nrow(fixture))
temporary <- tempfile("climate_group_test_")
dir.create(temporary)
tryCatch({
  tables <- write_climate_group_period_tables(x, temporary)
  z <- tables$FI_12_zone_climate_group_period
  stopifnot(all(abs(z[, sum(consensus_share), by = zoneID]$V1 - 1) < 1e-12))
  cat_table <- tables$FI_13_category_climate_group_period
  stopifnot(all(cat_table$n_zones == 2L),
            abs(sum(cat_table$consensus_mean) - 1) < 1e-12,
            abs(cat_table[variable_group == "Mean temperature", consensus_mean] - .3) < 1e-12)
}, finally = unlink(temporary, recursive = TRUE))
cat("All variable-group and aggregation tests passed.\n")
