# Focused collaborator update from completed 3.22 outputs. Separate from 3.23.
# Default: main Multi-Forest, normal + late SSP245/585; full-resolution statistics.
# One read pass per period over UNIQUE population-source dual surfaces.
# No model fitting, national prediction, population sorting, or species TIFF writes.
# Existing project outputs are read only; report/cache/figures stay in a new folder.

eco_brief_plots <- function(out, normal_agreement, ecosystem_area, transition,
                            species_niche, population_niche, display,
                            plot_template, palette, species_info) {
  figure_dir <- file.path(out, "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  saved <- character()
  write_png <- function(name, width, height, draw) {
    file <- file.path(figure_dir, name)
    grDevices::png(file, width = width, height = height, units = "in", res = 180,
                   pointsize = 12, bg = "white")
    tryCatch(draw(), finally = grDevices::dev.off())
    if (!file.exists(file) || file.info(file)$size <= 0) stop("Empty figure: ", file)
    saved <<- c(saved, file)
    invisible(file)
  }
  dt <- data.table::as.data.table
  pal <- data.table::copy(dt(palette))
  pal[, zoneID := as.integer(zoneID)]
  pal <- unique(pal[, .(zoneID, zone = as.character(zone), COLOR = as.character(COLOR))])
  if (anyDuplicated(pal$zoneID) || anyNA(pal$COLOR)) stop("Invalid ecotype plotting palette.")
  pal <- pal[order(zoneID)]
  zone_codes <- pal$zoneID
  zone_colors <- pal$COLOR
  no_data_color <- "#DDDDDD"
  novel_color <- if (99L %in% zone_codes) pal[zoneID == 99L, COLOR][1] else "#1A1A1A"
  change_colors <- c(retained = "#2C7FB8", different_current = "#E69F00", novel = novel_color)
  change_labels <- c("Same ecotype niche", "Different current ecotype", "Novel ecosystem niche")
  scenario_label <- function(s) {
    if (identical(s, "normal")) return("Reference period")
    sub("([0-9]{4}-[0-9]{4})(SSP[0-9]+)", "\\1 | \\2", s)
  }
  scenario_title <- function(s) {
    if (identical(s, "normal")) return("Reference period")
    sub("([0-9]{4}-[0-9]{4})(SSP[0-9]+)", "\\1, \\2", s)
  }
  scenarios <- names(display)
  if (is.null(scenarios) || !"normal" %in% scenarios) stop("Display data need a named normal period.")
  futures <- setdiff(scenarios, "normal")
  all_species <- colnames(display[["normal"]]$species_suitability)
  if (is.null(all_species)) all_species <- unique(as.character(species_niche$Species))
  info <- dt(species_info)
  scientific <- if (all(c("SpeciesCode", "LatinName") %in% names(info))) {
    as.character(info$LatinName[match(all_species, info$SpeciesCode)])
  } else rep(NA_character_, length(all_species))
  scientific[is.na(scientific) | !nzchar(scientific)] <- all_species[is.na(scientific) | !nzchar(scientific)]
  names(scientific) <- all_species
  species_title <- function(s) paste(strwrap(scientific[[s]], width = 34), collapse = "\n")
  map_plot <- function(values, title, colors, breaks) {
    if (length(values) != terra::ncell(plot_template)) stop("Plot vector does not match display geometry.")
    graphics::par(mar = c(0.4, 0.4, 2.6, 0.4))
    if (!any(is.finite(values))) {
      graphics::plot.new()
      graphics::plot.window(xlim = c(terra::xmin(plot_template), terra::xmax(plot_template)),
                            ylim = c(terra::ymin(plot_template), terra::ymax(plot_template)))
      graphics::rect(terra::xmin(plot_template), terra::ymin(plot_template),
                     terra::xmax(plot_template), terra::ymax(plot_template), col = no_data_color, border = NA)
      graphics::title(main = title, cex.main = 1)
      return(invisible(NULL))
    }
    r <- terra::rast(plot_template)
    terra::values(r) <- as.numeric(values)
    terra::plot(r, col = colors, breaks = breaks, colNA = no_data_color,
                axes = FALSE, legend = FALSE, main = title,
                cex.main = 1, maxcell = terra::ncell(r))
    invisible(NULL)
  }
  source_indices <- function(z, suitable = NULL) {
    ans <- match(as.integer(z), zone_codes)
    if (!is.null(suitable)) {
      ans[is.finite(suitable) & suitable < 0.4] <- 0L
      ans[!is.finite(suitable)] <- NA_integer_
    }
    ans
  }
  color_key <- function(codes, name, title) {
    key <- pal[zoneID %in% codes]
    if (!nrow(key)) return(invisible(NULL))
    nc <- if (nrow(key) > 20) 3L else if (nrow(key) > 8) 2L else 1L
    nr <- ceiling(nrow(key) / nc)
    write_png(name, 17, max(4.5, nr * 0.48 + 1.2), function() {
      graphics::par(mar = c(1.7, 0.4, 2.5, 0.4))
      graphics::plot.new(); graphics::plot.window(xlim = c(0, nc), ylim = c(0, nr))
      for (i in seq_len(nrow(key))) {
        col <- (i - 1L) %/% nr; row <- (i - 1L) %% nr
        y <- nr - row - 0.5
        graphics::rect(col + 0.02, y - 0.17, col + 0.13, y + 0.17,
                       col = key$COLOR[i], border = "#777777")
        label <- paste0("Zone ", key$zoneID[i], ": ", gsub("_", " ", key$zone[i]))
        graphics::text(col + 0.17, y, labels = paste(strwrap(label, width = 49), collapse = "\n"),
                       adj = c(0, 0.5), cex = 0.72)
      }
      graphics::title(main = title, cex.main = 1.15)
      graphics::mtext("Colours identify reference ecotypes; species x ecotype combinations are population proxies.",
                      side = 1, line = 0.8, cex = 0.75)
    })
  }

  agreement <- data.table::copy(dt(normal_agreement))[order(k)]
  if (any(agreement$pixel_agreement > 1.000001 | agreement$area_agreement > 1.000001, na.rm = TRUE)) {
    stop("Agreement plotting expects proportions between zero and one.")
  }
  write_png("01_reference_topk_agreement.png", 8.2, 5.3, function() {
    graphics::par(mar = c(4.5, 4.8, 3, 1))
    graphics::plot(agreement$k, agreement$pixel_agreement * 100, type = "b", pch = 16,
                   col = "#2166AC", lwd = 2, ylim = c(0, 100), xaxt = "n",
                   xlab = "Reference-assisted Top-k rule", ylab = "Reference-map agreement (%)",
                   main = "Multi-Forest reference-map agreement")
    graphics::axis(1, at = agreement$k, labels = paste0("Top-", agreement$k))
    graphics::lines(agreement$k, agreement$area_agreement * 100, type = "b", pch = 17,
                    col = "#B35806", lwd = 2)
    graphics::legend("bottomright", c("Cells weighted equally", "Cells weighted by area"),
                     col = c("#2166AC", "#B35806"), lty = 1, pch = c(16, 17), bty = "n")
  })
  write_png("02_reference_and_future_top3_ecotype_maps.png", 16, 6.5, function() {
    graphics::par(mfrow = c(1, length(scenarios)), oma = c(3.2, 0, 2, 0))
    for (s in scenarios) map_plot(source_indices(display[[s]]$ecotype_top3), scenario_title(s),
                                  zone_colors, seq(0.5, length(zone_codes) + 0.5, by = 1))
    graphics::mtext("Multi-Forest, reference-assisted Top-3 ecotype niches", side = 3, outer = TRUE, cex = 1.2)
    graphics::mtext("Grey: unavailable or outside reference footprint. Ecotype colour key supplied separately.",
                    side = 1, outer = TRUE, line = 1, cex = 0.8)
  })
  color_key(zone_codes, "08_ecotype_colour_key.png", "Ecotype colour key")

  if (length(futures)) {
    write_png("03_future_top3_matched_change_maps.png", 13, 6.8, function() {
      graphics::layout(rbind(seq_along(futures), rep(length(futures) + 1L, length(futures))), heights = c(1, 0.18))
      for (s in futures) map_plot(display[[s]]$change_top3, scenario_title(s),
                                   unname(change_colors), seq(0.5, 3.5, by = 1))
      graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
      graphics::legend("center", c(change_labels, "Unavailable / outside footprint"),
                       fill = c(change_colors, no_data_color), ncol = 2, bty = "n", cex = 0.92)
    })
    tr <- dt(transition)[k == 3L]
    shares <- matrix(0, nrow = 3, ncol = length(futures), dimnames = list(names(change_colors), futures))
    for (s in futures) for (cl in names(change_colors)) {
      a <- tr[scenario == s & change_class == cl, area_share]
      if (length(a)) shares[cl, s] <- a[1] * 100
    }
    write_png("04_future_top3_change_area_shares.png", 11.5, 5.7, function() {
      graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(1, 0.16))
      graphics::par(mar = c(4.5, 6.5, 3, 1))
      mids <- graphics::barplot(shares, horiz = TRUE, col = change_colors, border = NA, xlim = c(0, 100),
        names.arg = gsub(", ", "\n", vapply(futures, scenario_title, character(1)), fixed = TRUE), las = 1,
        xlab = "Share of common reference-future comparison area (%)",
        main = "Multi-Forest Top-3 matched niche change")
      for (j in seq_along(futures)) {
        centers <- cumsum(shares[, j]) - shares[, j] / 2
        keep <- shares[, j] >= 4
        graphics::text(centers[keep], mids[j], sprintf("%.1f%%", shares[keep, j]), col = "white", cex = 0.9)
      }
      graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
      graphics::legend("center", legend = change_labels, fill = change_colors,
                       horiz = TRUE, bty = "n", cex = 0.8)
    })
  }

  sn <- data.table::copy(dt(species_niche))
  baseline <- sn[scenario == "normal", .(Species, reference_area_km2 = suitable_area_km2)]
  projected <- merge(sn[scenario %in% futures], baseline, by = "Species", all.x = TRUE, sort = FALSE)
  projected[, change_percent := ifelse(is.finite(reference_area_km2) & reference_area_km2 > 0,
    100 * (suitable_area_km2 / reference_area_km2 - 1), NA_real_)]
  change_matrix <- matrix(NA_real_, length(futures), length(all_species), dimnames = list(futures, all_species))
  for (s in futures) for (sp in all_species) {
    a <- projected[scenario == s & Species == sp, change_percent]
    if (length(a)) change_matrix[s, sp] <- a[1]
  }
  if (length(futures) && any(is.finite(change_matrix))) {
    write_png("05_species_suitable_niche_area_change.png", 12, max(6, length(all_species) * 0.62), function() {
      graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(1, 0.12))
      graphics::par(mar = c(4.5, 15, 3.5, 1))
      keep <- colSums(is.finite(change_matrix)) > 0
      x <- change_matrix[, keep, drop = FALSE]
      scenario_colors <- c("#4575B4", "#D73027", "#66A61E")[seq_along(futures)]
      xr <- range(c(0, x), na.rm = TRUE); span <- max(10, diff(xr)); xr <- xr + c(-0.04, 0.08) * span
      graphics::barplot(x, horiz = TRUE, beside = TRUE, border = NA, col = scenario_colors,
        names.arg = unname(scientific[colnames(x)]), las = 1, cex.names = 0.85,
        xlim = xr, xlab = "Change in suitable niche area relative to reference (%)",
        main = "Species envelopes: maximum source-ecotype dual suitability >= 0.40")
      graphics::abline(v = 0, col = "#555555")
      graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
      graphics::legend("center", legend = vapply(futures, scenario_title, character(1)),
        fill = scenario_colors, bty = "n", horiz = TRUE, cex = 0.8)
    })
  }

  ns <- length(all_species)
  if (ns) {
    nc <- min(5L, ceiling(sqrt(ns * 1.7))); nr <- ceiling(ns / nc)
    suit_colors <- rev(grDevices::hcl.colors(64, "YlGnBu"))
    for (s in scenarios) {
      zz <- display[[s]]$species_suitability
      pp <- display[[s]]$population_best
      if (nrow(zz) != terra::ncell(plot_template) || !all(all_species %in% colnames(zz))) {
        stop("Invalid species display matrix for ", s)
      }
      write_png(paste0("06_species_suitable_niches_", s, ".png"), 18, max(5.5, nr * 3.5 + 1.1), function() {
        graphics::par(mfrow = c(nr, nc), oma = c(4, 0, 3, 0))
        for (sp in all_species) {
          value <- zz[, sp]
          level <- rep(NA_integer_, length(value))
          level[is.finite(value) & value < 0.4] <- 0L
          good <- is.finite(value) & value >= 0.4
          level[good] <- 1L + pmin(63L, floor((pmin(1, value[good]) - 0.4) / 0.6 * 64))
          map_plot(level, species_title(sp), c("white", suit_colors), seq(-0.5, 64.5, by = 1))
        }
        if (nr * nc > ns) for (i in seq_len(nr * nc - ns)) graphics::plot.new()
        graphics::mtext(paste0("Continuous species niche envelopes | ", scenario_title(s)),
          side = 3, outer = TRUE, cex = 1.15, line = 0.8)
        graphics::par(fig = c(0.36, 0.64, 0.045, 0.075), mar = c(0, 0, 0, 0), new = TRUE)
        graphics::plot.new(); graphics::plot.window(xlim = c(0.4, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
        at <- seq(0.4, 1, length.out = 65)
        graphics::rect(at[-65], 0, at[-1], 1, col = suit_colors, border = NA)
        graphics::axis(1, at = c(0.4,0.6,0.8,1), cex.axis = 0.8, tcl = -0.2, mgp = c(1,0.3,0))
        graphics::mtext("Dual suitability; white: below 0.40; grey: unavailable / outside footprint",
          side = 1, outer = TRUE, line = 0.6, cex = 0.8)
      })
      if (!is.null(pp) && all(all_species %in% colnames(pp))) {
        write_png(paste0("07_most_suitable_source_population_", s, ".png"), 18,
          max(5.5, nr * 3.5 + 1.1), function() {
            graphics::par(mfrow = c(nr, nc), oma = c(3.2, 0, 3, 0))
            for (sp in all_species) map_plot(source_indices(pp[, sp], zz[, sp]), species_title(sp),
              c("white", zone_colors), seq(-0.5, length(zone_codes) + 0.5, by = 1))
            if (nr * nc > ns) for (i in seq_len(nr * nc - ns)) graphics::plot.new()
            graphics::mtext(paste0("Most suitable source-ecotype population | ", scenario_title(s)),
              side = 3, outer = TRUE, cex = 1.15, line = 0.8)
            graphics::mtext("Maximum dual suitability >= 0.40. White: below threshold; grey: unavailable. Source colour key supplied separately.",
              side = 1, outer = TRUE, line = 0.6, cex = 0.8)
          })
      }
    }
    color_key(unique(as.integer(population_niche$source_zone)), "09_population_source_colour_key.png",
              "Source-ecotype population colour key")
  }
  invisible(saved)
}


local({
  run_brief <- function() {
    for (p in c("terra", "data.table")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
    library(terra); library(data.table)
    root <- normalizePath(Sys.getenv("ECOCHINA2_DIR", getwd()), winslash = "/", mustWork = TRUE)
    while (!file.exists(file.path(root, "script", "1. data.R"))) {
      parent <- dirname(root); if (parent == root) stop("Run inside ecoChina2 or set ECOCHINA2_DIR."); root <- parent
    }
    signature <- function(files) {
      if (!length(files) || !all(file.exists(files))) return(NA_character_)
      z <- file.info(files)
      paste(paste(normalizePath(files, winslash = "/", mustWork = TRUE), z$size,
        sprintf("%.6f", as.numeric(z$mtime)), sep = "|"), collapse = "\n")
    }
    policy <- "complete_climate_soil_candidates_v1"
    completion <- file.path(root, "assessment/missing_climate_repair_v1/repair_completed.rds")
    if (!file.exists(completion)) stop("Complete 3.22 first.")
    state <- readRDS(completion)
    if (!identical(state$policy, policy) || !isTRUE(state$sensitivity_ready) ||
        !identical(state$signature, signature(state$input_files))) stop("3.22 outputs are incomplete or stale.")
    method <- Sys.getenv("ECOCHINA2_BRIEF_METHOD", "mf_var")
    if (!method %in% c("mf_var", "rf_var")) stop("Choose mf_var or rf_var.")
    scenarios <- c("normal", "2071-2100SSP245", "2071-2100SSP585")
    zones <- c(1:7, 9:50, 52:55); ks <- c(1L, 3L, 5L); threshold <- 0.4
    out <- Sys.getenv("ECOCHINA2_BRIEF_OUT", file.path(root, "assessment/collaborator_update_fast_v1"))
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    cache <- file.path(out, "cache"); dir.create(cache, showWarnings = FALSE)
    old_options <- terraOptions(print = FALSE)
    terraOptions(tempdir = cache, memfrac = 0.10)
    on.exit({
      td <- old_options$tempdir; if (is.null(td) || !dir.exists(td)) td <- tempdir()
      terraOptions(tempdir = td, memfrac = old_options$memfrac)
    }, add = TRUE)
    ref_file <- file.path(root, "raster/ecosys_ori.tif"); ref <- rast(ref_file)
    if (!is.lonlat(ref) || !grepl("WGS.?84|World Geodetic System 1984", crs(ref))) {
      stop("This focused script expects the verified WGS84 reference grid.")
    }
    pop_file <- file.path(root, "future tree niche var/tables/population_projection_lookup_var.csv")
    if (!file.exists(pop_file)) pop_file <- file.path(root, "species_zone_population_long.csv")
    pop <- fread(pop_file)
    if (!all(c("PopulationID", "Species") %in% names(pop))) stop("Invalid population table.")
    if (!"source_zone" %in% names(pop)) {
      if (!"Zone" %in% names(pop)) stop("Population table lacks source_zone or Zone.")
      pop[, source_zone := as.integer(sub("^Zone", "", Zone))]
    }
    pop[, `:=`(Species = as.character(Species), PopulationID = as.character(PopulationID), source_zone = as.integer(source_zone))]
    if ("projected" %in% names(pop)) pop <- pop[projected %in% c(TRUE, "TRUE", "true", 1, "1")]
    pop <- pop[source_zone %in% zones]
    if (!nrow(pop) || anyNA(pop[, .(PopulationID, Species, source_zone)]) ||
        anyDuplicated(pop[, .(Species, source_zone)])) stop("Invalid/duplicate projected populations.")
    setorder(pop, Species, source_zone)
    species <- sort(unique(pop$Species)); sources <- sort(unique(pop$source_zone))
    species_columns <- lapply(species, function(s) match(pop[Species == s, source_zone], sources))
    names(species_columns) <- species
    palette_file <- file.path(root, "color_palette_China.csv"); palette <- fread(palette_file)
    if (!all(c("zoneID", "zone", "COLOR") %in% names(palette))) stop("Invalid ecotype colour palette.")
    palette[, zoneID := as.integer(zoneID)]
    if (!99L %in% palette$zoneID) palette <- rbindlist(list(palette,
      data.table(zoneID = 99L, zone = "Novel ecosystems", COLOR = "#555555")), fill = TRUE)
    if (any(!sources %in% palette$zoneID)) stop("Population source ecotypes lack palette entries.")
    label_file <- file.path(root, "species_population_summary.csv")
    species_info <- if (file.exists(label_file)) fread(label_file) else data.table(
      SpeciesCode = c("querVar", "pinuYun", "lariGme", "phylPub", "lariOlg", "robiPse", "chamFor", "cyclLon", "pinuSyl", "saliMat"),
      LatinName = c("Quercus variabilis", "Pinus yunnanensis", "Larix gmelinii", "Phyllostachys pubescens",
        "Larix olgensis var. changpaiensis", "Robinia pseudoacacia", "Chamaecyparis formosensis",
        "Cyclobalanopsis longinux", "Pinus sylvestris var. mongolica", "Salix matsudana"),
      CommonName = c("Chinese Cork Oak", "Yunnan Pine", "Dahurian Larch", "Moso Bamboo", "Hinggan Larch",
        "Black Locust", "Taiwan Cypress", "Longinux Oak", "Mongolian Pine", "Chinese Willow"))
    if (!all(c("SpeciesCode", "LatinName", "CommonName") %in% names(species_info))) {
      species_info <- data.table(SpeciesCode = species, LatinName = species, CommonName = species)
    }
    unknown <- setdiff(species, species_info$SpeciesCode)
    if (length(unknown)) species_info <- rbindlist(list(species_info,
      data.table(SpeciesCode = unknown, LatinName = unknown, CommonName = unknown)), fill = TRUE)
    map_file <- function(s) file.path(root, "result maps", method, paste0("assigned_zone_", s,
      "_threshold0.4_tol1e-04_novel99_maskNA8_noNovelNormal.tif"))
    rank_file <- function(s) file.path(root, "dual suit ranking var", method, s, "ranked_zone.tif")
    dual_files <- function(s) file.path(root, "dual suit", method, s, paste0("dual_suitability_zone", sources, ".tif"))
    inputs <- unique(c(ref_file, pop_file, palette_file, if (file.exists(label_file)) label_file,
      unlist(lapply(scenarios, function(s) c(map_file(s), rank_file(s), dual_files(s))), use.names = FALSE)))
    if (!all(file.exists(inputs))) stop("Missing repaired input(s):\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
    for (file in setdiff(inputs, c(pop_file, palette_file, label_file))) {
      if (!compareGeom(rast(file)[[1]], ref, stopOnError = FALSE)) stop("Geometry mismatch: ", file)
    }
    if (any(vapply(scenarios, function(s) nlyr(rast(rank_file(s))) < 5L, logical(1)))) stop("Five saved ranks are required.")
    input_signature <- signature(inputs)
    plot_factor <- as.integer(Sys.getenv("ECOCHINA2_BRIEF_PLOT_FACTOR", "12"))
    if (is.na(plot_factor) || plot_factor < 1L) stop("Invalid plot factor.")
    config <- list(policy = policy, method = method, scenarios = scenarios, sources = sources,
      population_keys = paste(pop$Species, pop$source_zone), threshold = threshold, ks = ks, plot_factor = plot_factor)
    run_state_file <- file.path(cache, "run_state.rds")
    if (file.exists(run_state_file)) {
      prior <- readRDS(run_state_file)
      if (!identical(prior$config, config) || !identical(prior$input_signature, input_signature)) {
        stop("Report inputs changed. Set ECOCHINA2_BRIEF_OUT to a new report directory.")
      }
    } else saveRDS(list(config = config, input_signature = input_signature), run_state_file)
    unlink(file.path(out, "RUN_COMPLETE.txt"))
    # Exact WGS84 cell area depends on row latitude and cell width, not longitude.
    area_column <- rast(nrows = nrow(ref), ncols = 1L, xmin = xmin(ref), xmax = xmin(ref) + res(ref)[1],
      ymin = ymin(ref), ymax = ymax(ref), crs = crs(ref))
    row_area <- as.numeric(values(cellSize(area_column, unit = "km")))
    nr_block <- max(1L, min(64L, floor(6e6 / (ncol(ref) * (length(sources) + 12L)))))
    plot_template <- rast(nrows = ceiling(nrow(ref) / plot_factor), ncols = ceiling(ncol(ref) / plot_factor),
      xmin = xmin(ref), xmax = xmin(ref) + ceiling(ncol(ref) / plot_factor) * res(ref)[1] * plot_factor,
      ymax = ymax(ref), ymin = ymax(ref) - ceiling(nrow(ref) / plot_factor) * res(ref)[2] * plot_factor, crs = crs(ref))
    plot_cells <- cellFromXY(ref, xyFromCell(plot_template, seq_len(ncell(plot_template))))
    plot_index_by_block <- split(which(is.finite(plot_cells)),
      ceiling(rowFromCell(ref, plot_cells[is.finite(plot_cells)]) / nr_block))
    normal_cache <- file.path(cache, "normal_topk.tif")
    footprint_zones <- c(zones, 99L)
    all_stats <- list(); display <- list()
    tag <- function(s) gsub("[^A-Za-z0-9_-]", "_", s)
    cat(sprintf("[BRIEF] %s; %d species, %d populations, %d unique dual surfaces, 3 periods.\n",
                method, length(species), nrow(pop), length(sources)))
    for (scenario in scenarios) {
      checkpoint <- file.path(cache, paste0(tag(scenario), "_done.rds"))
      if (file.exists(checkpoint)) {
        old <- readRDS(checkpoint)
        if (identical(old$input_signature, input_signature) &&
            (scenario != "normal" || file.exists(normal_cache))) {
          all_stats[[scenario]] <- old$stats; display[[scenario]] <- old$display
          cat("[BRIEF REUSE]", scenario, "\n"); next
        }
      }
      started <- proc.time()[[3]]
      cat("[BRIEF COMPUTE]", scenario, "\n"); flush.console()
      normal_mode <- identical(scenario, "normal")
      input <- c(ref, rast(map_file(scenario)), rast(rank_file(scenario))[[1:5]], rast(dual_files(scenario)))
      if (!normal_mode) input <- c(input, rast(normal_cache))
      readStart(input)
      map_temp <- NULL; out_normal <- NULL
      if (normal_mode) {
        map_temp <- file.path(cache, "normal_topk_partial.tif")
        out_normal <- rast(ref, nlyrs = 3L); names(out_normal) <- paste0("top", ks)
        writeStart(out_normal, map_temp, overwrite = TRUE, datatype = "INT1U", NAflag = 255,
          gdal = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER"))
      }
      e_cells <- e_area <- matrix(0, nrow = length(footprint_zones), ncol = 3L)
      t_cells <- t_area <- matrix(0, nrow = 3L, ncol = 3L)
      agreement_cells <- agreement_area <- numeric(3L)
      compared_cells <- compared_area <- eligible_cells <- eligible_area <- 0
      stat_names <- c("valid_cells", "valid_area_km2", "suitable_cells", "suitable_area_km2",
        "suitability_weighted_area_km2", "sum_dual_suitability", "max_dual_suitability")
      source_stats <- matrix(0, nrow = length(sources), ncol = length(stat_names), dimnames = list(NULL, stat_names))
      species_stats <- matrix(0, nrow = length(species), ncol = length(stat_names), dimnames = list(NULL, stat_names))
      source_stats[, "max_dual_suitability"] <- species_stats[, "max_dual_suitability"] <- -Inf
      view <- list(ecotype_top3 = rep(NA_real_, ncell(plot_template)),
        change_top3 = rep(NA_real_, ncell(plot_template)),
        species_suitability = matrix(NA_real_, ncell(plot_template), length(species), dimnames = list(NULL, species)),
        population_best = matrix(NA_real_, ncell(plot_template), length(species), dimnames = list(NULL, species)))
      starts <- seq.int(1L, nrow(ref), by = nr_block)
      next_progress <- 0L
      tryCatch({
        for (bi in seq_along(starts)) {
          start <- starts[bi]; nr <- min(nr_block, nrow(ref) - start + 1L)
          v <- readValues(input, row = start, nrows = nr, mat = TRUE)
          observed <- v[, 1]; assigned <- v[, 2]; ranks <- v[, 3:7, drop = FALSE]
          ar <- rep(row_area[start:(start + nr - 1L)], each = ncol(ref))
          labels <- matrix(assigned, nrow(v), 3L)
          for (ki in 2:3) {
            hit <- rowSums(ranks[, seq_len(ks[ki]), drop = FALSE] == observed, na.rm = TRUE) > 0L
            restore <- is.finite(assigned) & assigned != 99 & observed %in% zones & hit
            labels[restore, ki] <- observed[restore]
          }
          labels[!is.finite(labels)] <- NA_real_
          if (normal_mode) writeValues(out_normal, labels, start = start, nrows = nr)
          plot_rows <- plot_index_by_block[[as.character(bi)]]
          plot_local <- if (length(plot_rows)) plot_cells[plot_rows] - (start - 1L) * ncol(ref) else integer()
          if (length(plot_rows)) view$ecotype_top3[plot_rows] <- labels[plot_local, 2L]
          if (normal_mode) {
            eligible <- observed %in% zones & is.finite(ar)
            common <- eligible & assigned %in% zones
            eligible_cells <- eligible_cells + sum(eligible); eligible_area <- eligible_area + sum(ar[eligible])
            compared_cells <- compared_cells + sum(common); compared_area <- compared_area + sum(ar[common])
            for (ki in 1:3) {
              hit_match <- common & labels[, ki] == observed
              agreement_cells[ki] <- agreement_cells[ki] + sum(hit_match, na.rm = TRUE)
              agreement_area[ki] <- agreement_area[ki] + sum(ar[which(hit_match)])
            }
          }
          for (ki in 1:3) {
            zindex <- match(labels[, ki], footprint_zones)
            ok <- is.finite(zindex) & is.finite(ar)
            e_cells[, ki] <- e_cells[, ki] + tabulate(zindex[ok], nbins = length(footprint_zones))
            if (any(ok)) {
              a <- rowsum(ar[ok], zindex[ok], reorder = FALSE)
              e_area[as.integer(rownames(a)), ki] <- e_area[as.integer(rownames(a)), ki] + as.numeric(a)
            }
          }
          if (!normal_mode) {
            baseline <- v[, 7L + length(sources) + 1:3, drop = FALSE]
            for (ki in 1:3) {
              ok <- baseline[, ki] %in% zones & labels[, ki] %in% footprint_zones & is.finite(ar)
              code <- ifelse(labels[, ki] == 99, 3L, ifelse(labels[, ki] == baseline[, ki], 1L, 2L))
              t_cells[, ki] <- t_cells[, ki] + tabulate(code[ok], nbins = 3L)
              if (any(ok)) {
                a <- rowsum(ar[ok], code[ok], reorder = FALSE)
                t_area[as.integer(rownames(a)), ki] <- t_area[as.integer(rownames(a)), ki] + as.numeric(a)
              }
              if (ki == 2L && length(plot_rows)) {
                shown <- code[plot_local]; shown[!ok[plot_local]] <- NA_real_
                view$change_top3[plot_rows] <- shown
              }
            }
          }
          dual <- v[, 7L + seq_along(sources), drop = FALSE]
          dual[!is.finite(observed), ] <- NA_real_
          surface_stats <- function(x) {
            good <- is.finite(x) & is.finite(ar); suitable <- good & x >= threshold
            c(sum(good), sum(ar[good]), sum(suitable), sum(ar[suitable]),
              sum(x[good] * ar[good]), sum(x[good]), if (any(good)) max(x[good]) else -Inf)
          }
          for (zi in seq_along(sources)) {
            s <- surface_stats(dual[, zi]); source_stats[zi, 1:6] <- source_stats[zi, 1:6] + s[1:6]
            source_stats[zi, 7] <- max(source_stats[zi, 7], s[7])
          }
          for (si in seq_along(species)) {
            cols <- species_columns[[si]]
            best <- dual[, cols[1]]
            best_zone <- rep(sources[cols[1]], length(best))
            if (length(cols) > 1L) for (zi in cols[-1]) {
              better <- is.finite(dual[, zi]) & (!is.finite(best) | dual[, zi] > best)
              best[better] <- dual[better, zi]; best_zone[better] <- sources[zi]
            }
            s <- surface_stats(best); species_stats[si, 1:6] <- species_stats[si, 1:6] + s[1:6]
            species_stats[si, 7] <- max(species_stats[si, 7], s[7])
            if (length(plot_rows)) {
              shown <- best[plot_local]; source_shown <- best_zone[plot_local]
              source_shown[!is.finite(shown) | shown < threshold] <- NA_real_
              shown[!is.finite(shown)] <- NA_real_
              view$species_suitability[plot_rows, si] <- shown
              view$population_best[plot_rows, si] <- source_shown
            }
          }
          progress <- floor(100 * bi / length(starts))
          if (progress >= next_progress || bi == length(starts)) {
            elapsed <- proc.time()[[3]] - started
            remain <- if (bi > 0) elapsed * (length(starts) - bi) / bi else NA_real_
            cat(sprintf("  %s %3d%%; elapsed %.1f min; this-period remainder ~%.1f min\n", scenario,
                        progress, elapsed / 60, remain / 60)); flush.console()
            next_progress <- progress + 10L
          }
        }
      }, finally = {
        readStop(input)
        if (normal_mode && !is.null(out_normal)) writeStop(out_normal)
      })
      if (normal_mode) {
        if (file.exists(normal_cache)) unlink(normal_cache)
        if (!file.rename(map_temp, normal_cache)) stop("Cannot finalize local normal-map cache.")
      }
      area <- rbindlist(lapply(1:3, function(ki) data.table(method = method, scenario = scenario,
        k = ks[ki], zone = footprint_zones, cells = e_cells[, ki], area_km2 = e_area[, ki],
        domain = "finite corrected assigned map; unmodeled observed classes retain Top-1")))
      transitions <- if (!normal_mode) rbindlist(lapply(1:3, function(ki) data.table(method = method,
        scenario = scenario, k = ks[ki], change_class = c("retained", "different_current", "novel"),
        cells = t_cells[, ki], area_km2 = t_area[, ki], common_cells = sum(t_cells[, ki]),
        common_area_km2 = sum(t_area[, ki]), pixel_share = t_cells[, ki] / sum(t_cells[, ki]),
        area_share = t_area[, ki] / sum(t_area[, ki])))) else data.table()
      agreement <- if (normal_mode) data.table(method = method, scenario = scenario, k = ks,
        eligible_reference_cells = eligible_cells, eligible_reference_area_km2 = eligible_area,
        compared_cells = compared_cells, compared_area_km2 = compared_area,
        matching_cells = agreement_cells, matching_area_km2 = agreement_area,
        pixel_agreement = agreement_cells / compared_cells, area_agreement = agreement_area / compared_area) else data.table()
      summarize_surfaces <- function(x) {
        a <- as.data.table(x)
        a[, mean_dual_suitability := fifelse(valid_cells > 0, sum_dual_suitability / valid_cells, NA_real_)]
        a[!is.finite(max_dual_suitability), max_dual_suitability := NA_real_]
        a[, `:=`(method = method, scenario = scenario, suitability_threshold = threshold,
          domain = "reference non-NA; finite source surface / maximum of available species sources")]
        a[]
      }
      by_source <- summarize_surfaces(source_stats)
      by_source[, source_zone := sources]
      populations <- merge(pop, by_source, by = "source_zone", all.x = TRUE, sort = FALSE)
      species_table <- summarize_surfaces(species_stats)
      species_table[, `:=`(Species = species, populations_projected = lengths(species_columns))]
      stats <- list(agreement = agreement, area = area, transition = transitions,
                    species = species_table, population = populations)
      if (!identical(input_signature, signature(inputs))) stop("Inputs changed during the report; do not use this incomplete report.")
      saveRDS(list(input_signature = input_signature, stats = stats, display = view), checkpoint)
      all_stats[[scenario]] <- stats; display[[scenario]] <- view
      cat(sprintf("[BRIEF PERIOD COMPLETE] %s in %.1f min\n", scenario, (proc.time()[[3]] - started) / 60))
      rm(input, v, dual); gc()
    }
    normal_agreement <- rbindlist(lapply(all_stats, `[[`, "agreement"), fill = TRUE)
    ecosystem_area <- rbindlist(lapply(all_stats, `[[`, "area"), fill = TRUE)
    transition <- rbindlist(lapply(all_stats, `[[`, "transition"), fill = TRUE)
    species_niche <- rbindlist(lapply(all_stats, `[[`, "species"), fill = TRUE)
    population_niche <- rbindlist(lapply(all_stats, `[[`, "population"), fill = TRUE)
    stopifnot(all(diff(normal_agreement$pixel_agreement) >= -1e-12),
      all(diff(normal_agreement$area_agreement) >= -1e-12),
      all(is.finite(transition$area_share)), all(transition$common_cells > 0))
    fwrite(normal_agreement, file.path(out, "brief_reference_agreement.csv"))
    fwrite(ecosystem_area, file.path(out, "brief_ecotype_area.csv"))
    fwrite(transition, file.path(out, "brief_matched_future_change.csv"))
    fwrite(species_niche, file.path(out, "brief_continuous_species_niche.csv"))
    fwrite(population_niche, file.path(out, "brief_population_niche.csv"))
    fwrite(pop, file.path(out, "brief_population_sources.csv"))
    fwrite(species_info, file.path(out, "brief_species_labels.csv"))
    info <- file.info(inputs)
    fwrite(data.table(path = normalizePath(inputs, winslash = "/"), size_bytes = as.numeric(info$size),
      modified = as.character(info$mtime)), file.path(out, "input_files.csv"))
    context <- file.path(out, "completed_results"); dir.create(context, showWarnings = FALSE)
    completed_files <- c(file.path(root, "assessment", c("soil_gate_sensitivity.csv", "soil_gate_sensitivity_common_domain.csv")),
      list.files(file.path(root, "assessment/missing_climate_repair_v1"), "[.]csv$", full.names = TRUE),
      list.files(file.path(root, "assessment_var"), "^normal_map_.*[.]csv$", full.names = TRUE),
      list.files(file.path(root, "future tree niche var/tables"), "[.]csv$", full.names = TRUE))
    for (f in completed_files[file.exists(completed_files)]) file.copy(f, context, overwrite = TRUE)
    cat("[BRIEF FIGURES] Full-resolution statistics saved; rendering sampled map displays.\n")
    if (exists("eco_brief_plots", mode = "function")) {
      eco_brief_plots(out, normal_agreement, ecosystem_area, transition, species_niche,
        population_niche, display, plot_template, palette, species_info)
    } else stop("Plotting function missing from the delivered script.")
    fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits, big.mark = ",")
    summary <- c("COLLABORATOR UPDATE FROM REPAIRED OUTPUTS", "",
      paste("Model:", if (method == "mf_var") "Plain RF Multi-Forest (main workflow)" else "Plain RF"),
      "Scope: reference period and 2071-2100 SSP245/SSP585. Statistics use full-resolution cells.",
      "Displays are sampled on a coarser grid; display samples are never used to estimate areas.", "",
      "REFERENCE-MAP AGREEMENT:",
      sprintf("Compared cells: %s; comparison area: %s km2.", fmt(normal_agreement$compared_cells[1], 0), fmt(normal_agreement$compared_area_km2[1])),
      vapply(seq_len(nrow(normal_agreement)), function(i) sprintf("Top-%d: pixel agreement %s%%; area agreement %s%%.",
        normal_agreement$k[i], fmt(100 * normal_agreement$pixel_agreement[i]), fmt(100 * normal_agreement$area_agreement[i])), character(1)),
      "Top-3/5 restore the observed ecotype when it occurs among the leading suitability ranks.",
      "They are reference-assisted map diagnostics, not independent validation accuracies.", "",
      "LATE-CENTURY MATCHED TOP-3 CHANGE:")
    for (s in scenarios[-1]) {
      a <- transition[scenario == s & k == 3L]
      summary <- c(summary, sprintf("%s: %s compared cells, %s km2.", s, fmt(a$common_cells[1], 0), fmt(a$common_area_km2[1])),
        sprintf("  %s: %s%% (%s km2).", a$change_class, fmt(100 * a$area_share), fmt(a$area_km2)))
    }
    summary <- c(summary, "", "POPULATION AND SPECIES NICHES:",
      sprintf("%d projected species x ecotype populations, %d species, %d distinct source ecotypes.", nrow(pop), length(species), length(sources)),
      "Population suitability equals the corrected dual surface of its source ecotype.",
      "Species suitability is the maximum available source-population suitability; all sources missing gives NA.",
      "Suitable area uses dual suitability >=0.40. These continuous envelopes are distinct from assigned Top-k niche maps.",
      "Each source/species uses its own finite-data domain; valid cells and valid area are reported for every period.",
      "Area changes across native period domains should be interpreted alongside those valid-domain columns.",
      "Population-colour maps show the source ecotype with the greatest modeled suitability, not abundance or ecological dominance.", "",
      "REPAIR EXPLANATION:",
      "The old soil-gate branch could turn missing climate suitability into zero, incorrectly admitting cells into maps/novel-area statistics.",
      "The corrected workflow preserves missing required inputs as NA. Missing data are not novel ecosystems.",
      "Fitted models, raw climate/soil inputs, and held-out performance were not changed by this repair.",
      "Independent diagnostics verified all30 soil-layer grids and found raw climate NoData at the sampled problem points.",
      "The upstream cause of the raw climate mask is still unresolved; a specific water/coastal cause has not been established.", "",
      "STATUS:",
      "The full12-stage refresh may still be running. This focused report is complete only for the explicitly listed scope.",
      "The source completed_results folder contains existing completed repair/normal-map/assigned-niche tables.")
    writeLines(summary, file.path(out, "COLLABORATOR_UPDATE.txt"))
    if (!identical(input_signature, signature(inputs))) stop("Inputs changed while plotting; report not marked complete.")
    writeLines(c("Focused report complete.", paste("UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      "No fitted model or existing project output was overwritten."), file.path(out, "RUN_COMPLETE.txt"))
    cat("\n[BRIEF COMPLETE]", normalizePath(out, winslash = "/"), "\n")
    cat("Share top-level CSV/TXT, figures/, and completed_results/. Keep cache/ local.\n")
  }
  run_brief()
})
