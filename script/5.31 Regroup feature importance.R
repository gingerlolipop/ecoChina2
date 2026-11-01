# 5.31 Regroup existing feature-importance CSVs; no RF objects or rasters.
# Run this after installing the revised grouping helper and script 11.
library(data.table)

find_importance_project <- function(path = getwd()) {
  configured <- Sys.getenv("ECOCHINA2_DIR", unset = "")
  if (nzchar(configured)) return(normalizePath(configured, winslash = "/", mustWork = TRUE))
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "functions", "importance_variable_groups.R"))) return(path)
    parent <- dirname(path)
    if (parent == path) stop("Run from the project or set ECOCHINA2_DIR.")
    path <- parent
  }
}
project_root <- find_importance_project()
source(file.path(project_root, "functions", "importance_variable_groups.R"), local = TRUE)
table_dir <- file.path(project_root, "assessment", "feature_importance", "tables")
aligned_file <- file.path(table_dir, "FI_05_aligned_zone_profiles.csv")
if (!file.exists(aligned_file)) {
  stop("Missing ", aligned_file,
       ". Run the revised script/5.3 feature importance.R once; it reads existing fitted models.")
}
aligned_before <- fread(aligned_file)
required <- c("zoneID", "zone", "category", "category2", "niche", "variable",
              "variable_key", "binary_share", "multiclass_share", "consensus_share",
              "binary_selected", "multiclass_selected")
if (!all(required %in% names(aligned_before))) {
  stop("The aligned file is not the binary-Multi-Forest/multiclass format from script 5.3.")
}
if (anyDuplicated(aligned_before[, .(zoneID, niche, variable)])) {
  stop("Duplicate zone/niche/variable rows in aligned input.")
}
share_columns <- c("binary_share", "multiclass_share", "consensus_share")
if (any(!is.finite(as.matrix(aligned_before[, ..share_columns])))) {
  stop("Aligned importance shares contain non-finite values.")
}
# Save the previous annotations once. All raw importance values are preserved.
backup_dir <- file.path(table_dir, "before_climate_statistic_period_v2")
dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
for (name in list.files(table_dir, pattern = "^FI.*[.]csv$")) {
  backup_file <- file.path(backup_dir, name)
  if (!file.exists(backup_file)) file.copy(file.path(table_dir, name), backup_file)
}
unlink(file.path(table_dir, "FI_grouping_settings.csv"))
aligned <- reannotate_importance_table(aligned_before)
keys <- c("zoneID", "niche", "variable")
before <- copy(aligned_before)
after <- copy(aligned)
setorderv(before, keys)
setorderv(after, keys)
for (name in share_columns) {
  if (!identical(before[[name]], after[[name]])) stop("Regrouping changed ", name)
}
zone_metadata <- unique(aligned[, .(zoneID, zone, category, category2)])
standard_error <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x) / sqrt(length(x))
}
is_no_vegetation <- function(x) {
  tolower(gsub("[ -]+", "_", trimws(as.character(x)))) == "no_vegetation"
}

# 7. Sum variables into groups, then average zones within category =============

# Every variable is present in `aligned` for every zone, with zero importance
# when absent from a fitted model. Consequently the following sums do not omit
# unselected variables or small categories.
zone_group <- aligned[, .(
  n_variables = .N,
  n_binary_selected = sum(binary_selected),
  n_multiclass_selected = sum(multiclass_selected),
  binary_share = sum(binary_share),
  multiclass_share = sum(multiclass_share),
  consensus_share = sum(consensus_share)
), by = .(zoneID, zone, category, category2, niche, variable_group)]

category_group <- zone_group[, .(
  n_zones = .N,
  binary_mean = mean(binary_share),
  binary_se = standard_error(binary_share),
  multiclass_mean = mean(multiclass_share),
  multiclass_se = standard_error(multiclass_share),
  consensus_mean = mean(consensus_share),
  consensus_se = standard_error(consensus_share)
), by = .(category2, niche, variable_group)]

category_variable <- aligned[, .(
  n_zones = .N,
  binary_mean = mean(binary_share),
  binary_se = standard_error(binary_share),
  multiclass_mean = mean(multiclass_share),
  multiclass_se = standard_error(multiclass_share),
  consensus_mean = mean(consensus_share),
  consensus_se = standard_error(consensus_share),
  binary_selection_frequency = mean(binary_selected),
  multiclass_selection_frequency = mean(multiclass_selected)
), by = .(category2, niche, variable, variable_group, climate_period, climate_statistic, climate_period_detail)]
category_variable[, consensus_rank :=
                    frank(-consensus_mean, ties.method = "min"), by = .(category2, niche)]

zone_period <- aligned[
  niche == "climate" & !is.na(climate_period),
  .(
    binary_share = sum(binary_share),
    multiclass_share = sum(multiclass_share),
    consensus_share = sum(consensus_share)
  ),
  by = .(zoneID, zone, category, category2, climate_period)
]
category_period <- zone_period[, .(
  n_zones = .N,
  binary_mean = mean(binary_share),
  binary_se = standard_error(binary_share),
  multiclass_mean = mean(multiclass_share),
  multiclass_se = standard_error(multiclass_share),
  consensus_mean = mean(consensus_share),
  consensus_se = standard_error(consensus_share)
), by = .(category2, climate_period)]

setorder(zone_group, niche, category2, zoneID, -consensus_share)
setorder(category_group, niche, category2, -consensus_mean)
setorder(category_variable, niche, category2, consensus_rank, variable)
setorder(zone_period, category2, zoneID, climate_period)
setorder(category_period, category2, climate_period)

fwrite(zone_group, file.path(table_dir, "FI_07_zone_grouped_importance.csv"))
fwrite(category_group, file.path(table_dir, "FI_08_category_grouped_importance.csv"))
fwrite(category_variable, file.path(table_dir, "FI_09_category_variable_importance.csv"))
fwrite(zone_period, file.path(table_dir, "FI_10_zone_climate_period.csv"))
fwrite(category_period, file.path(table_dir, "FI_11_category_climate_period.csv"))


write_climate_group_period_tables(aligned, table_dir)


reader_zones <- zone_metadata[!is_no_vegetation(category2), zoneID]
outputs <- list(
  FI_05_aligned_zone_profiles = aligned,
  FI_reader_03_aligned_zone_profiles = aligned[zoneID %in% reader_zones],
  FI_reader_05_zone_grouped_importance = zone_group[zoneID %in% reader_zones],
  FI_reader_06_category_grouped_importance = category_group[!is_no_vegetation(category2)],
  FI_reader_07_category_variable_importance = category_variable[!is_no_vegetation(category2)],
  FI_reader_08_zone_climate_period = zone_period[zoneID %in% reader_zones],
  FI_reader_09_category_climate_period = category_period[!is_no_vegetation(category2)]
)
for (name in names(outputs)) fwrite(outputs[[name]], file.path(table_dir, paste0(name, ".csv")))

# Refresh annotations on raw and reader variable tables when they already exist.
# Agreement scores/ranks do not depend on the physical grouping and are untouched.
raw_files <- c("FI_02_binary_mf_importance.csv", "FI_03_multiclass_global_importance.csv",
               "FI_04_multiclass_class_importance.csv", "FI_reader_01_binary_mf_importance.csv",
               "FI_reader_02_multiclass_class_importance.csv")
for (name in raw_files) {
  path <- file.path(table_dir, name)
  if (file.exists(path)) fwrite(reannotate_importance_table(fread(path)), path)
}
lookup <- unique(annotate_importance_variables(aligned$niche, aligned$variable))
setorder(lookup, niche, variable_group, climate_period, variable)
fwrite(lookup, file.path(table_dir, "FI_00_variable_group_lookup.csv"))
fwrite(data.table(grouping_version = importance_grouping_version),
       file.path(table_dir, "FI_grouping_settings.csv"))
audit <- aligned[, lapply(.SD, sum), by = .(zoneID, niche), .SDcols = share_columns]
fwrite(audit, file.path(table_dir, "FI_regrouping_share_audit.csv"))
cat("\nREGROUPING COMPLETE: importance shares unchanged.\n",
    "Tables: ", table_dir, "\n",
    "Next run script 11 with ECOCHINA2_FIGURE_SECTIONS=importance.\n", sep = "")
