# One-command entry point: validate the taxonomy, regroup cached importance,
# and redraw importance figures. Fitted models are read only if the cache is absent.
local({
  refresh_importance <- function() {
    find_root <- function(path = getwd()) {
      configured <- Sys.getenv("ECOCHINA2_DIR", unset = "")
      if (nzchar(configured)) return(normalizePath(configured, winslash = "/", mustWork = TRUE))
      path <- normalizePath(path, winslash = "/", mustWork = TRUE)
      repeat {
        if (file.exists(file.path(path, "functions", "importance_variable_groups.R"))) return(path)
        parent <- dirname(path)
        if (parent == path) stop("Run inside ecoChina2 or set ECOCHINA2_DIR.")
        path <- parent
      }
    }
    root <- find_root()
    env_names <- c("ECOCHINA2_DIR", "ECOCHINA2_FIGURE_SECTIONS", "ECOCHINA2_FIGURE_IDS")
    previous <- Sys.getenv(env_names, unset = NA_character_)
    on.exit({
      Sys.unsetenv(env_names[is.na(previous)])
      present <- previous[!is.na(previous)]
      if (length(present)) do.call(Sys.setenv, as.list(present))
    }, add = TRUE)
    Sys.setenv(ECOCHINA2_DIR = root)
    run_file <- function(name) {
      environment <- new.env(parent = globalenv())
      source(file.path(root, "script", name), local = environment)
      invisible(environment)
    }
    run_file("test_importance_variable_groups.R")
    aligned <- file.path(root, "assessment", "feature_importance", "tables",
                         "FI_05_aligned_zone_profiles.csv")
    if (file.exists(aligned)) {
      run_file("5.31 Regroup feature importance.R")
    } else {
      cat("No current aligned importance cache; reading fitted RF objects with script 5.3.\n")
      run_file("5.3 feature importance.R")
    }
    Sys.setenv(ECOCHINA2_FIGURE_SECTIONS = "importance", ECOCHINA2_FIGURE_IDS = "all")
    viz <- run_file("11. visualization.R")
    audit <- data.table::rbindlist(viz$figure_log, fill = TRUE)
    essential <- c("importance_category_groups", "importance_climate_group_period")
    completed <- audit[status == "complete", figure]
    if (!all(essential %in% completed)) {
      stop("Some required grouped figures were not regenerated; inspect figures/figure_log.csv.")
    }
    cat("\nUPDATED TABLES: assessment/feature_importance/tables/\n",
        "UPDATED FIGURES: figures/importance/\n",
        "Please push both the tables and figures for manuscript revision.\n", sep = "")
  }
  refresh_importance()
})
