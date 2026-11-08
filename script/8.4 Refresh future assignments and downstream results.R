# Run from the ecoChina2 root after copying this update into script/functions.
# Reuses suitability and raw ranks. Never trains models or predicts suitability.
# Completed stages and per-scenario raster caches survive an interrupted run.
local({
  refresh_future_assignments <- function() {
    root <- normalizePath(Sys.getenv("ECOCHINA2_DIR", getwd()), winslash = "/", mustWork = TRUE)
    while (!file.exists(file.path(root, "script", "8. future top k.R"))) {
      parent <- dirname(root)
      if (parent == root) stop("Run inside ecoChina2 or set ECOCHINA2_DIR.")
      root <- parent
    }
    helper <- file.path(root, "functions", "future_topk_rule.R")
    source(helper, local = environment())
    report <- file.path(root, "assessment", "future_assignment_refresh")
    logs <- file.path(report, "logs")
    checkpoints <- file.path(report, "checkpoints", eco_topk_rule_version)
    for (d in c(report, logs, checkpoints)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    progress <- file.path(report, "latest_progress.txt")
    status_file <- file.path(report, "stage_status.csv")
    old_wd <- getwd(); setwd(root); on.exit(setwd(old_wd), add = TRUE)
    settings <- c("ECOCHINA2_DIR", "ECOCHINA2_PROGRESS_FILE", "ECOCHINA2_FORCE_RANKING",
      "ECOCHINA2_FORCE_MAPS", "ECOCHINA2_REUSE_MAPS", "ECOCHINA2_REUSE_POPULATION",
      "ECOCHINA2_FORCE_MULTICLASS", "ECOCHINA2_FIGURE_SECTIONS", "ECOCHINA2_FIGURE_IDS",
      "ECOCHINA2_TOPK_MAP_K", "ECOCHINA2_BRIEF_METHOD", "ECOCHINA2_BRIEF_OUT")
    old_env <- Sys.getenv(settings, unset = NA_character_)
    on.exit({
      for (j in seq_along(settings)) {
        if (is.na(old_env[[j]])) Sys.unsetenv(settings[[j]]) else
          do.call(Sys.setenv, setNames(list(old_env[[j]]), settings[[j]]))
      }
    }, add = TRUE)
    Sys.setenv(ECOCHINA2_DIR = root, ECOCHINA2_PROGRESS_FILE = progress,
      ECOCHINA2_FORCE_RANKING = "false", ECOCHINA2_FORCE_MAPS = "false",
      ECOCHINA2_REUSE_MAPS = "true", ECOCHINA2_REUSE_POPULATION = "true",
      ECOCHINA2_FORCE_MULTICLASS = "false", ECOCHINA2_FIGURE_SECTIONS = "maps,niche,topk",
      ECOCHINA2_FIGURE_IDS = "all", ECOCHINA2_TOPK_MAP_K = "all",
      ECOCHINA2_BRIEF_METHOD = "mf_var",
      ECOCHINA2_BRIEF_OUT = file.path(root, "assessment", "collaborator_update_fast_v1"))
    eco_topk_milestone("REFRESH START", eco_topk_rule_version)
    cat("Live progress: ", progress, "\nDetailed logs: ", logs, "\n", sep = "")

    zones <- c(1:7, 9:50, 52:55)
    scenarios <- c("normal", "2011-2040SSP245", "2041-2070SSP245", "2071-2100SSP245",
      "2011-2040SSP585", "2041-2070SSP585", "2071-2100SSP585")
    legacy <- dir.exists(file.path(root, "dual suit", "rf_var"))
    methods <- if (legacy) c("rf_var", "mf_var") else "mf_var"
    dual_files <- unlist(lapply(methods, function(m) unlist(lapply(scenarios, function(s)
      file.path(root, "dual suit", m, s, paste0("dual_suitability_zone", zones, ".tif"))))))
    rank_files <- unlist(lapply(methods, function(m) unlist(lapply(scenarios, function(s)
      file.path(root, "dual suit ranking var", m, s,
        c("ranked_zone.tif", "ranked_suitability.tif", "ranked_summary.tif"))))))
    missing <- c(dual_files, rank_files)[!file.exists(c(dual_files, rank_files))]
    if (length(missing)) stop("Required existing suitability/rank checkpoints are missing.\n",
      "This refresh will not regenerate them:\n", paste(head(missing, 20L), collapse = "\n"))
    for (p in c("terra", "data.table", "ggplot2", "patchwork", "scales", "circlize", "jsonlite")) {
      if (!requireNamespace(p, quietly = TRUE)) stop("Missing R package: ", p)
    }
    stage_scripts <- c(
      if (legacy) "4.1 Assign zones for rf_var and mf_var.R" else "4. assign zone.R",
      if (legacy) "6.1 Future results for rf_var and mf_var.R" else "6. future pop & species niche.R",
      "7.1 Refresh assigned-map multiclass comparison.R", "8. future top k.R",
      if (legacy) c("8.2future top k.R", "11.2 Visualization for corrected climate models.R",
        "8.3 Matched Top-k reference-to-future comparison.R") else character(),
      "11. visualization.R",
      if (legacy) c("11.4 Top-k downstream figures.R", "11.5 Matched Top-k and novel-ecosystem figures.R",
        "11.7 Manuscript Top-k diagnostics and population-coloured species maps.R") else character())
    if (file.exists(file.path(root, "assessment", "missing_climate_repair_v1", "repair_completed.rds"))) {
      stage_scripts <- c(stage_scripts, "quick_collaborator_update.R")
    } else {
      eco_topk_milestone("REPORT NOT PRESENT", "No previous collaborator repair report; canonical results will be refreshed.")
    }
    files <- file.path(root, "script", stage_scripts)
    if (!all(file.exists(files))) stop("Missing script(s):\n", paste(files[!file.exists(files)], collapse = "\n"))
    writeLines(stage_scripts, file.path(report, "run_order.txt"))
    # File signatures keep checkpoint decisions cheap: no national raster values
    # are read merely to decide whether a completed stage needs to run again.
    signature <- function(paths) {
      paths <- sort(unique(paths))
      if (!length(paths) || !all(file.exists(paths))) return(NA_character_)
      z <- file.info(paths)
      paste(paste(normalizePath(paths, winslash = "/"), z$size,
        sprintf("%.9f", as.numeric(z$mtime)), sep = "|"), collapse = "\n")
    }
    map_files <- unlist(lapply(methods, function(m) file.path(root, "result maps", m,
      paste0("assigned_zone_", scenarios, "_threshold0.4_tol1e-04_novel99_maskNA8_noNovelNormal.tif"))))
    original_files <- file.path(root, c("raster/ecosys_ori.tif", "color_palette_China.csv",
      "species_zone_population_long.csv", "species_population_summary.csv"))
    original_files <- original_files[file.exists(original_files)]
    external_dirs <- file.path(root, c("result maps/multiclass_rf", "multiclass/maps",
      "future tree niche dual suitability var/tables"))
    repair_record <- file.path(root, "assessment", "missing_climate_repair_v1", "repair_completed.rds")
    external_inputs <- function() {
      paths <- unique(unlist(lapply(external_dirs, list.files,
        pattern = "\\.(tif|csv)$", recursive = TRUE, full.names = TRUE)))
      c(paths, repair_record[file.exists(repair_record)])
    }
    helper_files <- list.files(file.path(root, "functions"), pattern = "\\.R$", full.names = TRUE)
    input_signature <- function(i) {
      # Stage 1 creates the maps, so only its existing source inputs enter its key.
      inputs <- c(files[seq_len(i)], helper_files, dual_files, rank_files, original_files)
      if (i > 1L) inputs <- c(inputs, map_files)
      # Canonical stage 6 may itself create continuous-area tables. Capture
      # these only after it completes, avoiding self-invalidating checkpoints.
      if (i > 2L) inputs <- c(inputs, external_inputs())
      # Include completed upstream analysis products. This detects a changed
      # source CSV even when the national raster inputs themselves are stable.
      if (i > 1L) for (j in seq_len(i - 1L)) {
        if (!grepl("^[4678]", stage_scripts[[j]])) next
        prior_file <- file.path(checkpoints, paste0("stage", j, ".rds"))
        if (file.exists(prior_file)) {
          prior <- tryCatch(readRDS(prior_file), error = function(e) NULL)
          if (!is.null(prior)) inputs <- c(inputs, prior$output_files)
        }
      }
      signature(inputs)
    }
    output_dirs <- function(script) {
      if (grepl("^4", script)) return(file.path(root, "result maps", methods))
      if (grepl("^6", script)) return(file.path(root, c("future tree niche var", "future tree niche dual suitability var")))
      if (grepl("^7", script)) return(file.path(root, "assessment", "multiclass_rf"))
      if (identical(script, "8. future top k.R")) return(file.path(root, "assessment", "topk"))
      if (grepl("^8.2", script)) return(file.path(root, "assessment_var", "future_topk_analysis"))
      if (grepl("^8.3", script)) return(file.path(root, "assessment_var", "future_topk_matched"))
      if (identical(script, "11. visualization.R")) return(file.path(root, "figures"))
      if (grepl("^quick", script)) return(file.path(root, "assessment", "collaborator_update_fast_v1"))
      file.path(root, "visualization var threshold0.4")
    }
    snapshot <- function(dirs) {
      paths <- unique(unlist(lapply(dirs, list.files, recursive = TRUE, full.names = TRUE)))
      paths <- paths[file.exists(paths) & !file.info(paths)$isdir]
      info <- file.info(paths)
      setNames(paste(info$size, sprintf("%.9f", as.numeric(info$mtime))), paths)
    }
    status <- data.frame(stage = seq_along(stage_scripts), script = stage_scripts,
      status = "pending", started = "", finished = "", log = "", stringsAsFactors = FALSE)
    save_status <- function() write.csv(status, status_file, row.names = FALSE)
    save_status()
    for (i in seq_along(stage_scripts)) {
      script <- stage_scripts[[i]]
      checkpoint <- file.path(checkpoints, paste0("stage", i, ".rds"))
      current_key <- input_signature(i)
      saved <- if (file.exists(checkpoint)) tryCatch(readRDS(checkpoint), error = function(e) NULL) else NULL
      current <- !is.null(saved) && identical(saved$script, script) &&
        identical(saved$inputs, current_key) && length(saved$output_files) > 0L &&
        identical(saved$output_signature, signature(saved$output_files))
      if (isTRUE(current)) {
        status$status[[i]] <- "reused"; status$finished[[i]] <- as.character(Sys.time()); save_status()
        eco_topk_milestone(sprintf("MILESTONE %02d/%02d REUSED", i, length(stage_scripts)), script)
        next
      }
      log_file <- file.path(logs, sprintf("%02d_%s.log", i, gsub("[^A-Za-z0-9]+", "_", script)))
      status$status[[i]] <- "running"; status$started[[i]] <- as.character(Sys.time())
      status$log[[i]] <- log_file; save_status()
      eco_topk_milestone(sprintf("MILESTONE %02d/%02d START", i, length(stage_scripts)), script)
      before <- snapshot(output_dirs(script))
      run_stage <- function() {
        depth <- sink.number(type = "output")
        sink(log_file, append = TRUE, split = TRUE)
        on.exit(while (sink.number(type = "output") > depth) sink(), add = TRUE)
        cat("\nStarted: ", as.character(Sys.time()), "\nScript: ", script, "\n", sep = "")
        # Legacy scripts clear their evaluation environment. Keep the driver
        # outside it; the global search path supplies their attached packages.
        stage_env <- new.env(parent = globalenv())
        withCallingHandlers(source(files[[i]], local = stage_env, encoding = "UTF-8"),
          warning = function(w) cat("[WARNING] ", conditionMessage(w), "\n", sep = ""),
          message = function(m) cat("[MESSAGE] ", conditionMessage(m), sep = ""))
        rm(stage_env); gc()
        cat("Finished: ", as.character(Sys.time()), "\n", sep = "")
      }
      failure <- tryCatch({run_stage(); NULL}, error = function(e) e)
      if (!is.null(failure)) {
        status$status[[i]] <- "failed"; status$finished[[i]] <- as.character(Sys.time()); save_status()
        eco_topk_milestone("REFRESH STOPPED", paste(script, conditionMessage(failure)))
        cat("\n[ERROR] ", conditionMessage(failure), "\n", file = log_file, append = TRUE, sep = "")
        stop("Stopped at ", script, ".\n", conditionMessage(failure),
          "\nSee ", log_file, ". Completed checkpoints are retained; rerun this driver to resume.", call. = FALSE)
      }
      after <- snapshot(output_dirs(script))
      changed <- names(after)[is.na(before[names(after)]) | before[names(after)] != after]
      outputs <- if (length(changed)) changed else names(after)
      if (!length(outputs)) stop("No outputs recorded for stage: ", script)
      saveRDS(list(script = script, inputs = input_signature(i), output_files = outputs,
        output_signature = signature(outputs), completed = Sys.time()), checkpoint)
      status$status[[i]] <- "complete"; status$finished[[i]] <- as.character(Sys.time()); save_status()
      eco_topk_milestone(sprintf("MILESTONE %02d/%02d COMPLETE", i, length(stage_scripts)), script)
    }
    writeLines(c(paste("Completed:", Sys.time()), paste("Rule:", eco_topk_rule_version),
      "Future observed-zone restoration requires rank <= k and own dual suitability >= 0.4.",
      "Reference-period rules and raw dual-suitability ranks are unchanged.",
      stage_scripts), file.path(report, "refresh_completed.txt"))
    eco_topk_milestone("REFRESH COMPLETE", report)
  }
  refresh_future_assignments()
})
