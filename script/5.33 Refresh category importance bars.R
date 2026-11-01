# Redraw only the requested horizontal category bars, using cached importance.
# No fitted RF objects, raster layers, or full script-11 startup are required.
local({
  refresh_category_bars <- function() {
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
    previous <- Sys.getenv("ECOCHINA2_DIR", unset = NA_character_)
    on.exit({
      if (is.na(previous)) Sys.unsetenv("ECOCHINA2_DIR") else Sys.setenv(ECOCHINA2_DIR = previous)
    }, add = TRUE)
    Sys.setenv(ECOCHINA2_DIR = root)
    library(data.table)
    library(ggplot2)
    run_script <- function(name) {
      source(file.path(root, "script", name), local = new.env(parent = globalenv()))
    }
    aligned <- file.path(root, "assessment", "feature_importance", "tables",
                         "FI_05_aligned_zone_profiles.csv")
    if (!file.exists(aligned)) {
      stop("No aligned importance cache. Run script/5.32 Refresh grouped importance figures.R once first.")
    }
    run_script("test_importance_variable_groups.R")
    run_script("5.31 Regroup feature importance.R")
    source(file.path(root, "functions", "importance_variable_groups.R"), local = TRUE)
    source(file.path(root, "functions", "plot_importance_category_groups.R"), local = TRUE)
    grouped <- fread(file.path(root, "assessment", "feature_importance", "tables",
                               "FI_reader_06_category_grouped_importance.csv"))
    plots <- importance_category_group_plots(grouped)
    write_importance_category_group_plots(plots, root)
    cat("\nUPDATED: detailed horizontal climate bars and topsoil bars.\n",
        "Both figures/importance/ and the manuscript's existing figure directory are current.\n",
        "No models were refitted and no suitability rasters were recalculated.\n", sep = "")
  }
  refresh_category_bars()
})
