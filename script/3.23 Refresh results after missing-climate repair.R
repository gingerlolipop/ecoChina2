# Refresh map-derived analyses and plots after 3.22; no fitting/importance work.
# Whole stages checkpoint; 6.21 additionally checks each cached raster's inputs.
local({
  refresh <- function() {
    root <- normalizePath(Sys.getenv("ECOCHINA2_DIR", getwd()), winslash = "/", mustWork = TRUE)
    while (!file.exists(file.path(root, "script", "1. data.R"))) {
      q <- dirname(root); if (q == root) stop("Run inside ecoChina2."); root <- q
    }
    source(file.path(root, "functions", "dual_missing_policy.R"), local = environment())
    report <- file.path(root, "assessment", "missing_climate_repair_v1")
    completed <- file.path(report, "repair_completed.rds")
    if (!file.exists(completed)) stop("Complete script 3.22 first.")
    state <- readRDS(completed)
    if (!identical(state$policy, eco_missing_policy) ||
        !identical(state$signature, eco_input_signature(state$input_files)) ||
        length(state$zones) != 53L || length(state$methods) != 2L || length(state$scenarios) != 7L ||
        !isTRUE(state$sensitivity_ready)) stop("Repair is incomplete/stale; run 3.22 first.")
    stages <- c("5.11 Refresh map assessment after missing-climate repair.R",
      "6.1 Future results for rf_var and mf_var.R",
      "6.21 Future population and species niches from dual suitability.R",
      "8.1 Rank dual suitability by pixel.R", "8. future top k.R",
      "8.2future top k.R",
      "11.2 Visualization for corrected climate models.R",
      "8.3 Matched Top-k reference-to-future comparison.R",
      "11.3 Reference-period Top-k confusion flows on a common scale.R",
      "11.4 Top-k downstream figures.R", "11.5 Matched Top-k and novel-ecosystem figures.R",
      "11.7 Manuscript Top-k diagnostics and population-coloured species maps.R")
    files <- file.path(root, "script", stages)
    if (!all(file.exists(files))) stop("Missing refresh scripts: ", paste(files[!file.exists(files)], collapse = "\n"))
    output_roots <- c("assessment_var", "future tree niche var",
      "future tree niche dual suitability var", "dual suit ranking var", "assessment/topk",
      "assessment_var/future_topk_analysis", "visualization var threshold0.4",
      "assessment_var/future_topk_matched", rep("visualization var threshold0.4", 4L))
    key_inputs <- c(state$input_files, file.path(root, "functions", "dual_missing_policy.R"),
      file.path(root, "assessment", c("soil_gate_sensitivity.csv", "soil_gate_sensitivity_common_domain.csv")))
    old_wd <- getwd(); setwd(root); on.exit(setwd(old_wd), add = TRUE)
    old_env <- Sys.getenv("ECOCHINA2_DIR", unset = NA_character_)
    Sys.setenv(ECOCHINA2_DIR = root)
    on.exit(if (is.na(old_env)) Sys.unsetenv("ECOCHINA2_DIR") else Sys.setenv(ECOCHINA2_DIR = old_env), add = TRUE)
    dir.create(file.path(report, "refresh_checkpoints"), showWarnings = FALSE)
    force <- tolower(Sys.getenv("ECOCHINA2_FORCE_RESULT_REFRESH", "false")) %in% c("1", "true", "yes")
    for (i in seq_along(stages)) {
      checkpoint <- file.path(report, "refresh_checkpoints", paste0("stage", i, ".rds"))
      signature <- eco_input_signature(c(key_inputs, files[seq_len(i)]))
      current <- FALSE
      if (!force && file.exists(checkpoint)) {
        z <- readRDS(checkpoint)
        current <- identical(z$signature, signature) && all(file.exists(z$output_files)) &&
          identical(z$output_signature, eco_input_signature(z$output_files))
      }
      cat(sprintf("\n[REFRESH %d/%d] %s\n", i, length(stages), stages[i])); flush.console()
      if (isTRUE(current)) { cat("[CURRENT STAGE]\n"); next }
      # Isolated evaluation prevents legacy rm(list=ls()) from erasing the driver.
      stage_env <- new.env(parent = globalenv())
      tryCatch(source(files[i], local = stage_env, encoding = "UTF-8"), error = function(e) {
        stop("Stopped at ", stages[i], ": ", conditionMessage(e),
             "\nCompleted stages are retained; run this driver again after resolving the error.")
      })
      rm(stage_env); gc()
      ff <- list.files(file.path(root, output_roots[i]), recursive = TRUE, full.names = TRUE)
      ff <- ff[file.exists(ff) & !file.info(ff)$isdir]
      if (!length(ff)) stop("Stage produced no recorded outputs: ", stages[i])
      saveRDS(list(signature = signature, output_files = ff,
        output_signature = eco_input_signature(ff)), checkpoint)
    }
    writeLines(c("All refresh stages completed (refresh chain v2).", paste("Time:", Sys.time()),
      "Training, held-out tests, feature importance and multiclass fitting were not rerun.",
      stages), file.path(report, "refresh_completed.txt"))
    cat("\n[RESULTS REFRESH COMPLETE]\n")
  }
  refresh()
})
