# 5.3 Feature importance: binary Multi-Forest and multiclass RF
# =============================================================================
# Run after script 7. This script reads fitted models, writes complete audit and
# reader-facing tables, and creates the core feature-importance figures.
#
# Primary comparison
# ------------------
# Negative permutation importance is retained in the audit tables but set to
# zero before normalization. Binary importance is normalized within each
# zone-specific climate or soil model. Multiclass class-specific importance is
# first normalized within class, and is also re-normalized within niche for the
# like-for-like climate and soil comparisons. Variables absent from one model
# are explicitly assigned zero before agreement and consensus are calculated.
# =============================================================================

library(data.table)
library(ggplot2)
library(randomForest)

rm(list = ls())
gc()


# 0. Paths and settings =========================================================

find_project_root <- function(path = getwd()) {
  configured <- Sys.getenv("ECOCHINA2_DIR", unset = "")
  if (nzchar(configured)) {
    return(normalizePath(configured, winslash = "/", mustWork = TRUE))
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "script", "5.3 feature importance.R"))) return(path)
    parent <- dirname(path)
    if (parent == path) stop("Run inside the repository or set ECOCHINA2_DIR.")
    path <- parent
  }
}

project_root <- find_project_root()

zones <- c(1:7, 9:50, 52:55)

first_existing <- function(paths, label) {
  paths <- paths[nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop("Missing ", label, ":\n", paste(paths, collapse = "\n"))
  }
  hit[[1]]
}

# Keep the fitted objects produced by the original workflow as the primary
# inputs. The clean-branch RDS paths are accepted only as a compatibility
# fallback, so existing 53-zone fits are never duplicated merely because the
# scripts were reorganized.
binary_model_file <- function(niche, zone) {
  if (niche == "climate") {
    first_existing(
      c(
        file.path(project_root, "rf", paste0("clm_mfVar_zone", zone, ".Rdata")),
        file.path(project_root, "models", "climate", paste0("climate_mf_zone", zone, ".rds"))
      ),
      paste("climate Multi-Forest model for Zone", zone)
    )
  } else {
    first_existing(
      c(
        file.path(project_root, "rf_soil", paste0("soil_mf_zone", zone, ".Rdata")),
        file.path(project_root, "models", "soil", paste0("soil_mf_zone", zone, ".rds"))
      ),
      paste("soil Multi-Forest model for Zone", zone)
    )
  }
}

multiclass_model_file <- first_existing(
  c(
    file.path(project_root, "rf_multiclass", "multiclass_climate_soil_rf.Rdata"),
    file.path(project_root, "multiclass", "models", "multiclass_rf.rds")
  ),
  "multiclass RF model"
)
table_dir <- file.path(
  project_root, "assessment", "feature_importance", "tables"
)
figure_dir <- file.path(
  project_root, "assessment", "feature_importance", "figures"
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

palette_candidates <- c(
  file.path(project_root, "color_palette_China.csv"),
  file.path(project_root, "data", "zone_palette.csv")
)
palette_file <- palette_candidates[file.exists(palette_candidates)][1]
if (is.na(palette_file)) {
  stop("Zone palette not found. Run script/color_palette.R first.")
}


# 1. Helpers ====================================================================

read_rf <- function(file) {
  if (tolower(tools::file_ext(file)) == "rds") {
    objects <- list(readRDS(file))
  } else {
    environment <- new.env(parent = emptyenv())
    object_names <- load(file, envir = environment)
    objects <- mget(object_names, envir = environment, inherits = FALSE)
  }
  direct <- Filter(function(x) inherits(x, "randomForest"), objects)
  if (length(direct) == 1L) return(direct[[1]])
  candidates <- Filter(
    function(x) inherits(x, "randomForest"),
    unlist(lapply(objects, function(x) if (is.list(x)) x else list()), recursive = FALSE)
  )
  if (length(candidates) != 1L) {
    stop("Could not identify one randomForest object in: ", file)
  }
  candidates[[1]]
}

importance_matrix <- function(model) {
  answer <- randomForest::importance(model, scale = FALSE)
  if (is.null(dim(answer))) {
    answer <- matrix(
      answer,
      ncol = 1L,
      dimnames = list(names(answer), "MeanDecreaseAccuracy")
    )
  }
  answer
}

column_or_na <- function(x, column) {
  if (column %in% colnames(x)) as.numeric(x[, column]) else rep(NA_real_, nrow(x))
}

normalize_positive <- function(x) {
  positive <- pmax(x, 0)
  total <- sum(positive, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) rep(0, length(x)) else positive / total
}

safe_cor <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = method))
}

cosine_similarity <- function(x, y) {
  denominator <- sqrt(sum(x ^ 2)) * sqrt(sum(y ^ 2))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  sum(x * y) / denominator
}

standard_error <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x) / sqrt(length(x))
}

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

is_no_vegetation <- function(x) {
  tolower(gsub("[ -]+", "_", trimws(as.character(x)))) == "no_vegetation"
}

annotate_variables <- function(niche, variable) {
  niche <- as.character(niche)
  variable <- as.character(variable)
  standard <- toupper(gsub("-", "_", variable, fixed = TRUE))
  group <- rep("Other", length(variable))

  climate <- niche == "climate"
  soil <- niche == "soil"

  group[climate & grepl("^(MAT|MWMT|MCMT|TD$|EMT|EXT|TMAX|TMIN|TAVE)", standard)] <-
    "Temperature"
  group[climate & grepl("^(MAP|MSP|PPT)", standard)] <- "Precipitation"
  group[climate & grepl("^(AHM|SHM|CMD|CMI)", standard)] <-
    "Climatic moisture balance"
  group[climate & grepl("^(BFFP|EFFP|FFP|NFFD|DD)", standard)] <-
    "Growing season and degree-days"
  group[climate & grepl("^(EREF|RSDS)", standard)] <-
    "Radiation and evaporative demand"
  group[climate & grepl("^RH", standard)] <- "Humidity"
  group[climate & grepl("^PAS", standard)] <- "Snowfall"

  group[soil & grepl("CEC|TEB|BASE|(^|_)BS($|_)", standard)] <-
    "Nutrient retention and base status"
  group[soil & grepl("BLD|BULK|DENS", standard)] <- "Bulk density"
  group[soil & grepl("CRF|GRAVEL|COARSE", standard)] <- "Coarse fragments"
  group[soil & group == "Other" & grepl("SND|SLT|CLY|SAND|SILT|CLAY|TEXT", standard)] <-
    "Texture"
  group[soil & grepl("(^|_)OC($|_)|ORC|ORGANIC", standard)] <- "Organic carbon"
  group[soil & grepl("PH", standard)] <- "Soil reaction"
  group[soil & grepl("TCEQ|CACO3|CARBONATE", standard)] <- "Carbonates"
  group[soil & grepl("GY|CASO4", standard)] <- "Gypsum"
  group[soil & grepl("ESP|SOD", standard)] <- "Sodicity"
  group[soil & grepl("(^|_)ECE($|_)|SALIN|ELECTRICAL", standard)] <- "Salinity"

  suffix <- tolower(sub("^.*_(wt|sp|sm|at)$", "\\1", variable))
  seasonal <- grepl("_(wt|sp|sm|at)$", tolower(variable))
  period <- rep(NA_character_, length(variable))
  period[climate] <- "Annual"
  period[climate & seasonal & suffix == "wt"] <- "Winter"
  period[climate & seasonal & suffix == "sp"] <- "Spring"
  period[climate & seasonal & suffix == "sm"] <- "Summer"
  period[climate & seasonal & suffix == "at"] <- "Autumn"

  data.table(
    niche = niche,
    variable = variable,
    variable_group = group,
    climate_period = period
  )
}


# 2. Palette and expected model inventory ======================================

palette <- fread(palette_file)
required_palette <- c("zoneID", "zone", "category", "category2")
if (!all(required_palette %in% names(palette))) {
  stop("Palette must contain: ", paste(required_palette, collapse = ", "))
}
palette[, zoneID := as.integer(zoneID)]
zone_metadata <- unique(palette[zoneID %in% zones, ..required_palette])
if (nrow(zone_metadata) != length(zones)) {
  stop("Palette does not contain every modeled zone.")
}

binary_specs <- rbindlist(list(
  data.table(
    zoneID = zones,
    niche = "climate",
    model_file = vapply(zones, function(z) binary_model_file("climate", z), character(1))
  ),
  data.table(
    zoneID = zones,
    niche = "soil",
    model_file = vapply(zones, function(z) binary_model_file("soil", z), character(1))
  )
))


# 3. Binary Multi-Forest importance ============================================

extract_binary <- function(specification) {
  file <- specification$model_file
  if (!file.exists(file)) {
    return(list(
      importance = NULL,
      audit = data.table(
        workflow = "binary_mf", zoneID = specification$zoneID,
        niche = specification$niche, scope = "zone_model", model_file = file,
        file_exists = FALSE, load_ok = FALSE, n_predictors = NA_integer_,
        n_trees = NA_integer_, n_negative_mda = NA_integer_,
        positive_mda_total = NA_real_, normalized_sum = NA_real_,
        error_message = "model file not found"
      )
    ))
  }

  tryCatch({
    model <- read_rf(file)
    imp <- importance_matrix(model)
    if (!("MeanDecreaseAccuracy" %in% colnames(imp))) {
      stop("MeanDecreaseAccuracy is absent.")
    }
    variable <- rownames(imp)
    mda <- as.numeric(imp[, "MeanDecreaseAccuracy"])
    share <- normalize_positive(mda)
    result <- data.table(
      workflow = "binary_mf",
      scope = "zone_model",
      zoneID = as.integer(specification$zoneID),
      niche = specification$niche,
      variable = variable,
      variable_key = paste(specification$niche, variable, sep = "::"),
      mda_unscaled = mda,
      presence_mda_unscaled = column_or_na(imp, "1"),
      mean_decrease_gini = column_or_na(imp, "MeanDecreaseGini"),
      positive_mda = pmax(mda, 0),
      importance_share_model = share,
      importance_share_within_niche = share,
      importance_rank = frank(-mda, ties.method = "min"),
      selected = TRUE,
      model_file = file
    )
    list(
      importance = result,
      audit = data.table(
        workflow = "binary_mf", zoneID = specification$zoneID,
        niche = specification$niche, scope = "zone_model", model_file = file,
        file_exists = TRUE, load_ok = TRUE, n_predictors = nrow(imp),
        n_trees = as.integer(model$ntree), n_negative_mda = sum(mda < 0),
        positive_mda_total = sum(pmax(mda, 0)),
        normalized_sum = sum(share), error_message = NA_character_
      )
    )
  }, error = function(e) {
    list(
      importance = NULL,
      audit = data.table(
        workflow = "binary_mf", zoneID = specification$zoneID,
        niche = specification$niche, scope = "zone_model", model_file = file,
        file_exists = TRUE, load_ok = FALSE, n_predictors = NA_integer_,
        n_trees = NA_integer_, n_negative_mda = NA_integer_,
        positive_mda_total = NA_real_, normalized_sum = NA_real_,
        error_message = conditionMessage(e)
      )
    )
  })
}

binary_result <- lapply(seq_len(nrow(binary_specs)), function(i) {
  extract_binary(binary_specs[i])
})
binary_audit <- rbindlist(lapply(binary_result, `[[`, "audit"), fill = TRUE)
binary_importance <- rbindlist(
  lapply(binary_result, `[[`, "importance"), fill = TRUE
)


# 4. Multiclass global and class-specific importance ===========================

multiclass_audit <- data.table(
  workflow = "multiclass_rf", zoneID = NA_integer_, niche = "joint",
  scope = "global_and_class", model_file = multiclass_model_file,
  file_exists = file.exists(multiclass_model_file), load_ok = FALSE,
  n_predictors = NA_integer_, n_trees = NA_integer_,
  n_negative_mda = NA_integer_, positive_mda_total = NA_real_,
  normalized_sum = NA_real_,
  error_message = if (file.exists(multiclass_model_file)) {
    NA_character_
  } else {
    "model file not found"
  }
)

multiclass_importance <- NULL
if (file.exists(multiclass_model_file)) {
  multiclass_extraction <- tryCatch({
    model <- read_rf(multiclass_model_file)
    imp <- importance_matrix(model)
    if (!("MeanDecreaseAccuracy" %in% colnames(imp))) {
      stop("Global MeanDecreaseAccuracy is absent.")
    }
    missing_classes <- setdiff(as.character(zones), colnames(imp))
    if (length(missing_classes)) {
      stop("Missing class-specific importance for Zones: ",
           paste(missing_classes, collapse = ", "))
    }

    model_variable <- rownames(imp)
    niche <- ifelse(grepl("^soil_", model_variable), "soil", "climate")
    variable <- sub("^soil_", "", model_variable)

    global_mda <- as.numeric(imp[, "MeanDecreaseAccuracy"])
    global <- data.table(
      workflow = "multiclass_rf", scope = "global", zoneID = NA_integer_,
      niche = niche, model_variable = model_variable, variable = variable,
      variable_key = paste(niche, variable, sep = "::"),
      mda_unscaled = global_mda,
      positive_mda = pmax(global_mda, 0),
      importance_share_model = normalize_positive(global_mda),
      mean_decrease_gini = column_or_na(imp, "MeanDecreaseGini"),
      model_file = multiclass_model_file
    )
    global[, importance_share_within_niche :=
      normalize_positive(mda_unscaled), by = niche]
    global[, importance_rank := frank(-mda_unscaled, ties.method = "min")]

    class_specific <- rbindlist(lapply(zones, function(z) {
      class_mda <- as.numeric(imp[, as.character(z)])
      data.table(
        workflow = "multiclass_rf", scope = "class_specific", zoneID = z,
        niche = niche, model_variable = model_variable, variable = variable,
        variable_key = paste(niche, variable, sep = "::"),
        mda_unscaled = class_mda, positive_mda = pmax(class_mda, 0),
        importance_share_model = normalize_positive(class_mda),
        mean_decrease_gini = column_or_na(imp, "MeanDecreaseGini"),
        model_file = multiclass_model_file
      )
    }))
    class_specific[, importance_share_within_niche :=
      normalize_positive(mda_unscaled), by = .(zoneID, niche)]
    class_specific[, importance_rank :=
      frank(-mda_unscaled, ties.method = "min"), by = zoneID]

    list(
      importance = rbindlist(list(global, class_specific), fill = TRUE),
      audit = data.table(
        workflow = "multiclass_rf", zoneID = NA_integer_, niche = "joint",
        scope = "global_and_class", model_file = multiclass_model_file,
        file_exists = TRUE, load_ok = TRUE, n_predictors = nrow(imp),
        n_trees = as.integer(model$ntree),
        n_negative_mda = sum(global_mda < 0),
        positive_mda_total = sum(pmax(global_mda, 0)),
        normalized_sum = sum(normalize_positive(global_mda)),
        error_message = NA_character_
      )
    )
  }, error = function(e) {
    list(
      importance = NULL,
      audit = data.table(
        workflow = "multiclass_rf", zoneID = NA_integer_, niche = "joint",
        scope = "global_and_class", model_file = multiclass_model_file,
        file_exists = TRUE, load_ok = FALSE, n_predictors = NA_integer_,
        n_trees = NA_integer_, n_negative_mda = NA_integer_,
        positive_mda_total = NA_real_, normalized_sum = NA_real_,
        error_message = conditionMessage(e)
      )
    )
  })
  multiclass_importance <- multiclass_extraction$importance
  multiclass_audit <- multiclass_extraction$audit
}

model_audit <- rbindlist(list(binary_audit, multiclass_audit), fill = TRUE)
fwrite(model_audit, file.path(table_dir, "FI_01_model_audit.csv"))

if (any(model_audit$load_ok != TRUE) || is.null(multiclass_importance) ||
    !nrow(binary_importance)) {
  stop("Feature-importance extraction failed. See FI_01_model_audit.csv.")
}


# 5. Annotate and save the full extracted importance ===========================

annotation <- unique(annotate_variables(
  c(binary_importance$niche, multiclass_importance$niche),
  c(binary_importance$variable, multiclass_importance$variable)
))
binary_importance <- merge(
  binary_importance, annotation,
  by = c("niche", "variable"), all.x = TRUE, sort = FALSE
)
multiclass_importance <- merge(
  multiclass_importance, annotation,
  by = c("niche", "variable"), all.x = TRUE, sort = FALSE
)
binary_importance <- merge(
  binary_importance, zone_metadata,
  by = "zoneID", all.x = TRUE, sort = FALSE
)
multiclass_importance <- merge(
  multiclass_importance, zone_metadata,
  by = "zoneID", all.x = TRUE, sort = FALSE
)

setorder(binary_importance, niche, zoneID, importance_rank, variable)
setorder(multiclass_importance, scope, zoneID, importance_rank, variable)
fwrite(binary_importance, file.path(table_dir, "FI_02_binary_mf_importance.csv"))
fwrite(
  multiclass_importance[scope == "global"],
  file.path(table_dir, "FI_03_multiclass_global_importance.csv")
)
fwrite(
  multiclass_importance[scope == "class_specific"],
  file.path(table_dir, "FI_04_multiclass_class_importance.csv")
)


# 6. Align absent variables, agreement, and consensus ==========================

multiclass_class <- multiclass_importance[scope == "class_specific"]
variable_lookup <- unique(rbindlist(list(
  binary_importance[, .(niche, variable, variable_key, variable_group, climate_period)],
  multiclass_class[, .(niche, variable, variable_key, variable_group, climate_period)]
)))

alignment_template <- rbindlist(lapply(unique(variable_lookup$niche), function(n) {
  data.table(
    zoneID = rep(zones, each = variable_lookup[niche == n, .N]),
    variable_key = rep(variable_lookup[niche == n, variable_key], times = length(zones))
  )
}))
alignment_template <- merge(
  alignment_template, variable_lookup,
  by = "variable_key", all.x = TRUE, sort = FALSE
)

binary_for_merge <- binary_importance[, .(
  zoneID, variable_key,
  binary_mda = mda_unscaled,
  binary_share = importance_share_within_niche,
  binary_selected = TRUE
)]
multiclass_for_merge <- multiclass_class[, .(
  zoneID, variable_key,
  multiclass_mda = mda_unscaled,
  multiclass_share_model = importance_share_model,
  multiclass_share = importance_share_within_niche,
  multiclass_selected = TRUE
)]

aligned <- merge(
  alignment_template, binary_for_merge,
  by = c("zoneID", "variable_key"), all.x = TRUE, sort = FALSE
)
aligned <- merge(
  aligned, multiclass_for_merge,
  by = c("zoneID", "variable_key"), all.x = TRUE, sort = FALSE
)
for (column in c("binary_share", "multiclass_share", "multiclass_share_model")) {
  set(aligned, which(is.na(aligned[[column]])), column, 0)
}
for (column in c("binary_selected", "multiclass_selected")) {
  set(aligned, which(is.na(aligned[[column]])), column, FALSE)
}
# Consensus gives the two workflows equal weight within each niche; raw MDA is
# never averaged across model types because its scale is model-dependent.
aligned[, `:=`(
  consensus_share = (binary_share + multiclass_share) / 2,
  binary_rank = frank(-binary_share, ties.method = "min"),
  multiclass_rank = frank(-multiclass_share, ties.method = "min"),
  consensus_rank = frank(-consensus_share, ties.method = "min")
), by = .(zoneID, niche)]
aligned <- merge(aligned, zone_metadata, by = "zoneID", all.x = TRUE, sort = FALSE)

agreement <- aligned[, {
  binary_top5 <- variable_key[binary_share > 0 & binary_rank <= 5L]
  multiclass_top5 <- variable_key[multiclass_share > 0 & multiclass_rank <= 5L]
  union_top5 <- union(binary_top5, multiclass_top5)
  .(
    n_variables_aligned = .N,
    n_binary_selected = sum(binary_selected),
    n_multiclass_selected = sum(multiclass_selected),
    pearson = safe_cor(binary_share, multiclass_share, "pearson"),
    spearman = safe_cor(binary_share, multiclass_share, "spearman"),
    cosine = cosine_similarity(binary_share, multiclass_share),
    top5_jaccard = if (!length(union_top5)) NA_real_ else
      length(intersect(binary_top5, multiclass_top5)) / length(union_top5)
  )
}, by = .(zoneID, zone, category, category2, niche)]

setorder(aligned, niche, zoneID, consensus_rank, variable)
setorder(agreement, niche, category2, zoneID)
category_agreement <- agreement[, .(
  n_zones = .N,
  pearson_mean = mean_or_na(pearson),
  pearson_se = standard_error(pearson),
  spearman_mean = mean_or_na(spearman),
  spearman_se = standard_error(spearman),
  cosine_mean = mean_or_na(cosine),
  cosine_se = standard_error(cosine),
  top5_jaccard_mean = mean_or_na(top5_jaccard),
  top5_jaccard_se = standard_error(top5_jaccard)
), by = .(category2, niche)]
fwrite(aligned, file.path(table_dir, "FI_05_aligned_zone_profiles.csv"))
fwrite(agreement, file.path(table_dir, "FI_06_zone_agreement.csv"))
fwrite(
  category_agreement,
  file.path(table_dir, "FI_06b_category_agreement.csv")
)


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
), by = .(category2, niche, variable, variable_group, climate_period)]
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


# 8. Reader-facing tables: exclude the no-vegetation ecotype ===================

reader_zones <- zone_metadata[!is_no_vegetation(category2), zoneID]
reader_tables <- list(
  FI_reader_01_binary_mf_importance = binary_importance[zoneID %in% reader_zones],
  FI_reader_02_multiclass_class_importance = multiclass_class[zoneID %in% reader_zones],
  FI_reader_03_aligned_zone_profiles = aligned[zoneID %in% reader_zones],
  FI_reader_04_zone_agreement = agreement[zoneID %in% reader_zones],
  FI_reader_04b_category_agreement = category_agreement[
    !is_no_vegetation(category2)
  ],
  FI_reader_05_zone_grouped_importance = zone_group[zoneID %in% reader_zones],
  FI_reader_06_category_grouped_importance = category_group[
    !is_no_vegetation(category2)
  ],
  FI_reader_07_category_variable_importance = category_variable[
    !is_no_vegetation(category2)
  ],
  FI_reader_08_zone_climate_period = zone_period[zoneID %in% reader_zones],
  FI_reader_09_category_climate_period = category_period[
    !is_no_vegetation(category2)
  ]
)

for (name in names(reader_tables)) {
  fwrite(reader_tables[[name]], file.path(table_dir, paste0(name, ".csv")))
}


# 9. Core manuscript-ready figures ============================================

category_label <- function(x) {
  tools::toTitleCase(gsub("_", " ", as.character(x), fixed = TRUE))
}

percent_axis <- function(x) paste0(round(100 * x), "%")

agreement_plot_data <- copy(zone_group[
  zoneID %in% reader_zones & (binary_share > 0 | multiclass_share > 0)
])
agreement_plot_data[, category_label := category_label(category2)]

agreement_plot <- ggplot(
  agreement_plot_data,
  aes(
    x = binary_share,
    y = multiclass_share,
    color = category_label,
    size = consensus_share
  )
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
  geom_point(alpha = 0.60) +
  facet_wrap(~niche, nrow = 1) +
  scale_x_continuous(labels = percent_axis) +
  scale_y_continuous(labels = percent_axis) +
  scale_size_continuous(labels = percent_axis, range = c(1, 4)) +
  coord_equal() +
  labs(
    x = "Binary Multi-Forest grouped importance",
    y = "Multiclass RF grouped importance",
    color = "Vegetation category",
    size = "Consensus",
    title = "Agreement in ecotype-level permutation importance"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

agreement_figure <- file.path(
  figure_dir, "Figure_FI_1_binary_vs_multiclass_agreement.png"
)
ggsave(
  agreement_figure, agreement_plot,
  width = 12, height = 7, units = "in", dpi = 300, bg = "white"
)

group_figures <- character()
for (niche_value in c("climate", "soil")) {
  plot_data <- copy(category_group[
    niche == niche_value & !is_no_vegetation(category2)
  ])
  plot_data[, category_label := category_label(category2)]
  group_order <- plot_data[
    , .(total = sum(consensus_mean)), by = variable_group
  ][order(-total), variable_group]
  plot_data[, variable_group := factor(variable_group, levels = group_order)]
  plot_data[, category_label := factor(
    category_label, levels = rev(sort(unique(category_label)))
  )]

  group_plot <- ggplot(
    plot_data,
    aes(x = consensus_mean, y = category_label, fill = variable_group)
  ) +
    geom_col(width = 0.75) +
    scale_x_continuous(labels = percent_axis, expand = c(0, 0)) +
    labs(
      x = "Consensus permutation-importance share",
      y = NULL,
      fill = "Variable group",
      title = paste(
        if (niche_value == "climate") "Climatic" else "Topsoil",
        "importance by vegetation category"
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))

  group_file <- file.path(
    figure_dir,
    paste0("Figure_FI_2_category_variable_groups_", niche_value, ".png")
  )
  ggsave(
    group_file, group_plot,
    width = 11, height = 6.5, units = "in", dpi = 300, bg = "white"
  )
  group_figures <- c(group_figures, group_file)
}

period_plot_data <- copy(category_period[!is_no_vegetation(category2)])
period_plot_data[, category_label := category_label(category2)]
period_plot_data[, category_label := factor(
  category_label, levels = rev(sort(unique(category_label)))
)]
period_plot_data[, climate_period := factor(
  climate_period,
  levels = c("Annual", "Winter", "Spring", "Summer", "Autumn")
)]

period_plot <- ggplot(
  period_plot_data,
  aes(x = consensus_mean, y = category_label, fill = climate_period)
) +
  geom_col(width = 0.75) +
  scale_x_continuous(labels = percent_axis, expand = c(0, 0)) +
  labs(
    x = "Consensus climatic-importance share",
    y = NULL,
    fill = "Climate period",
    title = "Annual and seasonal climatic importance"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

period_figure <- file.path(
  figure_dir, "Figure_FI_3_category_climate_period.png"
)
ggsave(
  period_figure, period_plot,
  width = 10, height = 5.5, units = "in", dpi = 300, bg = "white"
)

fwrite(
  data.table(
    figure = c(agreement_figure, group_figures, period_figure),
    content = c(
      "Binary-versus-multiclass grouped-importance agreement",
      "Category-level climatic variable-group consensus",
      "Category-level topsoil variable-group consensus",
      "Category-level annual and seasonal climatic consensus"
    )
  ),
  file.path(table_dir, "FI_figure_inventory.csv")
)

cat(
  "\nSECTION 5.3 COMPLETE\n",
  "Full audit and reader-facing tables: ", table_dir, "\n",
  "Core figures: ", figure_dir, "\n",
  "Reader-facing summaries exclude the no-vegetation ecotype.\n",
  sep = ""
)
