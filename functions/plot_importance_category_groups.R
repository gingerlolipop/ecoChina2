# Main manuscript figure: preserve horizontal stacked bars and a bottom legend.
# Requires the shared importance_variable_groups.R helper in the same environment.
# Plot consensus_mean as exported by script 5.3/5.31; never renormalize by group.

importance_category_group_plots <- function(grouped) {
  grouped <- data.table::copy(data.table::as.data.table(grouped))
  required <- c("category2", "niche", "variable_group", "consensus_mean")
  if (!all(required %in% names(grouped))) {
    stop("Category importance CSV needs: ", paste(required, collapse = ", "))
  }
  if (anyDuplicated(grouped[, .(category2, niche, variable_group)])) {
    stop("Duplicated category/niche/variable-group rows in importance CSV.")
  }
  if (any(!is.finite(grouped$consensus_mean)) || any(grouped$consensus_mean < 0)) {
    stop("Category importance shares must be finite and non-negative.")
  }
  if (any(as.character(grouped$variable_group) == "Temperature")) {
    stop("Old Temperature grouping detected. Run revised script 5.31 first.")
  }
  category_order <- c(
    "alpine_vegetation", "cropland", "desert", "forest",
    "grassland_meadow_steppe", "scrub", "wetland"
  )
  category_labels <- c(
    alpine_vegetation = "Alpine vegetation", cropland = "Cropland",
    desert = "Desert", forest = "Forest",
    grassland_meadow_steppe = "Grassland, meadow\nand steppe",
    scrub = "Scrub", wetland = "Wetland"
  )
  plots <- list()
  for (niche_value in c("climate", "soil")) {
    plot_data <- data.table::copy(grouped[niche == niche_value])
    if (!nrow(plot_data)) stop("Missing category importance for ", niche_value)
    category_codes <- as.character(plot_data$category2)
    if (length(setdiff(category_codes, names(category_labels)))) {
      stop("Unexpected vegetation category: ",
           paste(setdiff(category_codes, names(category_labels)), collapse = ", "))
    }
    expected_groups <- if (niche_value == "climate") {
      importance_climate_groups
    } else importance_soil_groups
    group_values <- as.character(plot_data$variable_group)
    unexpected <- setdiff(group_values, expected_groups)
    if (length(unexpected)) stop("Unknown variable group: ", paste(unexpected, collapse = ", "))
    group_order <- intersect(expected_groups, unique(group_values))
    colours <- importance_semantic_colours[group_order]
    if (anyNA(colours)) stop("Missing colour for a variable group.")
    totals <- plot_data[, .(total = sum(consensus_mean)), by = category2]
    if (any(abs(totals$total - 1) > 1e-6)) {
      stop("Category bars do not sum to one. Check aligned CSV completeness; do not renormalize the plot.")
    }
    present_categories <- intersect(category_order, unique(category_codes))
    plot_data[, category_label := factor(
      category_labels[as.character(category2)],
      levels = rev(unname(category_labels[present_categories]))
    )]
    plot_data[, variable_group := factor(variable_group, levels = group_order)]
    plots[[niche_value]] <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = consensus_mean, y = category_label, fill = variable_group)
    ) +
      ggplot2::geom_col(
        position = ggplot2::position_stack(reverse = TRUE), width = 0.70,
        colour = "white", linewidth = 0.22, orientation = "y"
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq(0, 1, by = 0.2), labels = function(x) sprintf("%.0f", 100 * x),
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.65)) +
      ggplot2::scale_fill_manual(
        values = colours, breaks = group_order, drop = FALSE, name = NULL
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2, byrow = FALSE)) +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::labs(
        title = if (niche_value == "climate") "Climate groups" else "Topsoil groups",
        x = "Share of normalized permutation importance (%)", y = NULL
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_line(colour = "#EEEEF2", linewidth = 0.45),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "#333333", linewidth = 0.4),
        axis.text = ggplot2::element_text(colour = "#222222"),
        axis.text.y = ggplot2::element_text(size = 13),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 7)),
        legend.position = "bottom", legend.justification = "left",
        legend.text = ggplot2::element_text(size = 11.5),
        legend.key.width = grid::unit(0.45, "cm"),
        legend.key.height = grid::unit(0.40, "cm"),
        legend.box.margin = ggplot2::margin(t = 12),
        plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0),
        plot.title.position = "panel",
        plot.margin = ggplot2::margin(t = 12, r = 14, b = 10, l = 12)
      )
  }
  plots
}

write_importance_category_group_plots <- function(plots, project_root, dpi = 320) {
  output_dir <- file.path(project_root, "figures", "importance")
  # Existing LaTeX graphicspath searches this directory before root figures/.
  manuscript_dir <- file.path(
    project_root, "visualization var threshold0.4", "figures", "importance"
  )
  for (directory in c(output_dir, manuscript_dir)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  }
  outputs <- character()
  for (niche_value in c("climate", "soil")) {
    filename <- paste0("category_variable_groups_", niche_value, ".png")
    output <- file.path(output_dir, filename)
    ggplot2::ggsave(
      output, plot = plots[[niche_value]], width = 12,
      height = if (niche_value == "climate") 7.8 else 7.2,
      units = "in", dpi = dpi, bg = "white"
    )
    info <- file.info(output)
    if (is.na(info$size) || info$size <= 0) stop("Empty figure: ", output)
    manuscript_output <- file.path(manuscript_dir, filename)
    if (!file.copy(output, manuscript_output, overwrite = TRUE)) {
      stop("Could not update manuscript figure: ", manuscript_output)
    }
    if (!identical(unname(tools::md5sum(output)), unname(tools::md5sum(manuscript_output)))) {
      stop("Manuscript figure differs from the newly generated figure.")
    }
    outputs <- c(outputs, output, manuscript_output)
    cat("[SAVED]", normalizePath(output, winslash = "/"), "\n")
    cat("[UPDATED MANUSCRIPT COPY]", normalizePath(manuscript_output, winslash = "/"), "\n")
  }
  invisible(outputs)
}
