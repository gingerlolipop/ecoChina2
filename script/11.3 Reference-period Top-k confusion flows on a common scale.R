# 11.3 Reference-period Top-k agreement decomposition
# ==============================================================================
# This standalone visualization script reads the reference-period Top-k
# confusion table produced by script 11.2. It does not rerun niche models and
# does not overwrite any 11.2 output.
#
# Diagnostic definition inherited from 11.2:
#   Top-1: the existing assigned map, including the 1e-4 tie rule.
#   Top-3/Top-5: the observed ecotype is counted as recovered when it occurs
#                within the first k dual-suitability ranks.
#
# The figure emphasizes cumulative agreement with the reference ecotype. Each
# 100% bar partitions evaluated pixels by the first rank interval at which the
# observed ecotype is recovered. The small remainder outside Top-5 is shown in
# very light grey so that residual disagreement is not visually dominant.
# ==============================================================================

library(data.table)
library(ggplot2)


# 0. Paths and settings =========================================================

base_dir <- "H:/Jing/ecoChina2"

visualization_dir <- file.path(
  base_dir,
  "visualization var threshold0.4"
)

figure_dir <- file.path(
  visualization_dir,
  "figures"
)

table_dir <- file.path(
  visualization_dir,
  "tables"
)

confusion_file <- file.path(
  table_dir,
  "Figure_var_4_reference_topk_confusion_long.csv"
)

output_figure <- file.path(
  figure_dir,
  "Figure_var_4_reference_topk_agreement_decomposition.png"
)

output_plot_data <- file.path(
  table_dir,
  "Figure_var_4_reference_topk_agreement_decomposition_data.csv"
)

output_summary <- file.path(
  table_dir,
  "Figure_var_4_reference_topk_agreement_summary.csv"
)

method_order <- c(
  "rf_var",
  "mf_var"
)

default_method_labels <- c(
  rf_var = "Plain RF",
  mf_var = "Plain MF RF"
)

rank_cutoffs <- c(
  1L,
  3L,
  5L
)

segment_order <- c(
  "Top-1 exact agreement",
  "Additional agreement at ranks 2--3",
  "Additional agreement at ranks 4--5",
  "Not recovered within Top-5"
)

# Restrained grayscale consistent with the manuscript's model-comparison plots.
segment_fills <- c(
  "Top-1 exact agreement" = "#202020",
  "Additional agreement at ranks 2--3" = "#737373",
  "Additional agreement at ranks 4--5" = "#BDBDBD",
  "Not recovered within Top-5" = "#F2F2F2"
)

segment_text_colors <- c(
  "Top-1 exact agreement" = "white",
  "Additional agreement at ranks 2--3" = "white",
  "Additional agreement at ranks 4--5" = "black",
  "Not recovered within Top-5" = "black"
)


# 1. Helpers ====================================================================

require_file <- function(file) {
  if (!file.exists(file)) {
    stop(
      "Missing required file: ",
      file
    )
  }
  
  file
}


format_integer <- function(x) {
  format(
    round(x),
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


# 2. Read and validate the 11.2 source table ===================================

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

confusion_all <- fread(
  require_file(
    confusion_file
  )
)

required_columns <- c(
  "original_zone",
  "predicted_zone",
  "n",
  "method",
  "rank_cutoff"
)

missing_columns <- setdiff(
  required_columns,
  names(confusion_all)
)

if (length(missing_columns) > 0L) {
  stop(
    "The confusion table is missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

confusion_all[
  ,
  `:=`(
    original_zone = as.integer(original_zone),
    predicted_zone = as.integer(predicted_zone),
    n = as.numeric(n),
    method = as.character(method),
    rank_cutoff = as.integer(rank_cutoff)
  )
]

if (any(
  !is.finite(confusion_all$n) |
  confusion_all$n < 0
)) {
  stop(
    "Column 'n' must contain finite, non-negative pixel counts."
  )
}

missing_methods <- setdiff(
  method_order,
  unique(confusion_all$method)
)

if (length(missing_methods) > 0L) {
  stop(
    "The confusion table is missing required methods: ",
    paste(
      missing_methods,
      collapse = ", "
    )
  )
}

confusion <- confusion_all[
  method %in% method_order &
    rank_cutoff %in% rank_cutoffs,
  .(
    pixel_count = sum(
      n,
      na.rm = TRUE
    )
  ),
  by = .(
    method,
    rank_cutoff,
    original_zone,
    predicted_zone
  )
]

rank_grid <- CJ(
  method = method_order,
  rank_cutoff = rank_cutoffs,
  unique = TRUE
)

observed_rank_grid <- unique(
  confusion[
    ,
    .(
      method,
      rank_cutoff
    )
  ]
)

missing_method_ranks <- rank_grid[
  !observed_rank_grid,
  on = .(
    method,
    rank_cutoff
  )
]

if (nrow(missing_method_ranks) > 0L) {
  stop(
    "Missing method/rank combinations: ",
    paste(
      paste0(
        missing_method_ranks$method,
        " Top-",
        missing_method_ranks$rank_cutoff
      ),
      collapse = ", "
    )
  )
}


# 3. Calculate cumulative Top-k agreement ======================================

agreement_summary <- confusion[
  ,
  .(
    compared_pixels = sum(
      pixel_count
    ),
    agreement_pixels = sum(
      pixel_count[
        original_zone == predicted_zone
      ]
    )
  ),
  by = .(
    method,
    rank_cutoff
  )
]

agreement_summary <- agreement_summary[
  order(
    match(
      method,
      method_order
    ),
    rank_cutoff
  )
]

agreement_summary[
  ,
  agreement_share := agreement_pixels /
    compared_pixels
]

agreement_summary[
  ,
  method_label := default_method_labels[method]
]

# Use method labels stored by 11.2 when each method has one unambiguous label.
if ("method_label" %in% names(confusion_all)) {
  labels_from_file <- unique(
    confusion_all[
      method %in% method_order &
        !is.na(method_label) &
        nzchar(
          trimws(
            as.character(method_label)
          )
        ),
      .(
        method,
        method_label = as.character(method_label)
      )
    ]
  )
  
  unambiguous_labels <- labels_from_file[
    ,
    if (.N == 1L) {
      .(
        method_label = method_label[[1L]]
      )
    },
    by = method
  ]
  
  agreement_summary[
    unambiguous_labels,
    method_label := i.method_label,
    on = "method"
  ]
}

# Every rank must evaluate the same pixels within a method.
denominator_check <- agreement_summary[
  ,
  .(
    n_denominators = uniqueN(
      compared_pixels
    )
  ),
  by = method
]

if (any(denominator_check$n_denominators != 1L)) {
  stop(
    "The evaluated-pixel denominator changes across Top-k ranks."
  )
}

# Agreement must be monotonic because increasing k only relaxes recovery.
monotonic_check <- agreement_summary[
  order(rank_cutoff),
  .(
    monotonic = all(
      diff(
        agreement_pixels
      ) >= 0
    )
  ),
  by = method
]

if (any(!monotonic_check$monotonic)) {
  stop(
    "Top-k agreement is not monotonic for at least one method."
  )
}

# A specific off-diagonal flow must also remain unchanged or decrease as k
# increases. This catches malformed Top-k tables that could pass the overall
# monotonicity check through compensating changes among zone pairs.
error_pair_wide <- dcast(
  confusion[
    original_zone != predicted_zone
  ],
  method + original_zone + predicted_zone ~ rank_cutoff,
  value.var = "pixel_count",
  fun.aggregate = sum,
  fill = 0
)

if (any(
  error_pair_wide[["3"]] >
  error_pair_wide[["1"]] |
  error_pair_wide[["5"]] >
  error_pair_wide[["3"]]
)) {
  stop(
    "At least one off-diagonal zone pair increases as the Top-k criterion ",
    "is relaxed. Check the 11.2 diagnostic maps."
  )
}


# 4. Partition pixels by first recovery interval ===============================

agreement_wide <- dcast(
  agreement_summary,
  method + method_label + compared_pixels ~ rank_cutoff,
  value.var = c(
    "agreement_pixels",
    "agreement_share"
  )
)

required_wide_columns <- c(
  "agreement_pixels_1",
  "agreement_pixels_3",
  "agreement_pixels_5",
  "agreement_share_1",
  "agreement_share_3",
  "agreement_share_5"
)

missing_wide_columns <- setdiff(
  required_wide_columns,
  names(agreement_wide)
)

if (length(missing_wide_columns) > 0L) {
  stop(
    "Failed to construct Top-k agreement columns: ",
    paste(
      missing_wide_columns,
      collapse = ", "
    )
  )
}

segment_data <- rbindlist(
  lapply(
    seq_len(
      nrow(agreement_wide)
    ),
    function(i) {
      row_i <- agreement_wide[i]
      
      segment_pixels <- c(
        row_i$agreement_pixels_1,
        row_i$agreement_pixels_3 -
          row_i$agreement_pixels_1,
        row_i$agreement_pixels_5 -
          row_i$agreement_pixels_3,
        row_i$compared_pixels -
          row_i$agreement_pixels_5
      )
      
      data.table(
        method = row_i$method,
        method_label = row_i$method_label,
        compared_pixels = row_i$compared_pixels,
        segment = segment_order,
        segment_pixels = as.numeric(
          segment_pixels
        ),
        segment_share = as.numeric(
          segment_pixels /
            row_i$compared_pixels
        )
      )
    }
  )
)

segment_data[
  ,
  `:=`(
    method_order_index = match(
      method,
      method_order
    ),
    segment_order_index = match(
      segment,
      segment_order
    )
  )
]

setorder(
  segment_data,
  method_order_index,
  segment_order_index
)

if (any(segment_data$segment_pixels < 0)) {
  stop(
    "At least one agreement segment has a negative pixel count."
  )
}

segment_check <- segment_data[
  ,
  .(
    pixel_sum = sum(
      segment_pixels
    ),
    share_sum = sum(
      segment_share
    ),
    compared_pixels = unique(
      compared_pixels
    )
  ),
  by = method
]

if (any(
  abs(
    segment_check$pixel_sum -
    segment_check$compared_pixels
  ) > 0.5
)) {
  stop(
    "Agreement segments do not sum to the evaluated-pixel total."
  )
}

if (any(
  abs(
    segment_check$share_sum -
    1
  ) > 1e-12
)) {
  stop(
    "Agreement segment shares do not sum to one."
  )
}

# Draw rectangles directly so segment order and annotation positions are exact.
segment_data[
  ,
  `:=`(
    x_min = shift(
      cumsum(
        segment_share
      ),
      fill = 0
    ),
    x_max = cumsum(
      segment_share
    )
  ),
  by = method
]

segment_data[
  ,
  `:=`(
    x_mid = (
      x_min +
        x_max
    ) / 2,
    y_position = length(method_order) -
      method_order_index +
      1,
    segment = factor(
      segment,
      levels = segment_order
    )
  )
]

segment_data[
  ,
  label := sprintf(
    "%.1f%%",
    100 * segment_share
  )
]

segment_data[
  ,
  text_color := unname(
    segment_text_colors[
      as.character(segment)
    ]
  )
]

top5_annotations <- agreement_summary[
  rank_cutoff == 5L,
  .(
    method,
    method_label,
    top5_share = agreement_share,
    y_position = length(method_order) -
      match(
        method,
        method_order
      ) +
      1
  )
]

top5_annotations[
  ,
  annotation := sprintf(
    "%.1f%% within Top-5",
    100 * top5_share
  )
]


# 5. Plot cumulative reference-period agreement =================================

method_axis <- unique(
  segment_data[
    ,
    .(
      y_position,
      method,
      method_label
    )
  ]
)

setorder(
  method_axis,
  y_position
)

figure_topk_agreement <- ggplot() +
  geom_rect(
    data = segment_data,
    aes(
      xmin = 100 * x_min,
      xmax = 100 * x_max,
      ymin = y_position - 0.24,
      ymax = y_position + 0.24,
      fill = segment
    ),
    color = "black",
    linewidth = 0.28
  ) +
  geom_text(
    data = segment_data,
    aes(
      x = 100 * x_mid,
      y = y_position,
      label = label,
      color = text_color
    ),
    size = 3.45,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_segment(
    data = top5_annotations,
    aes(
      x = 100 * top5_share,
      xend = 100 * top5_share,
      y = y_position + 0.25,
      yend = y_position + 0.36
    ),
    linewidth = 0.35,
    color = "black"
  ) +
  geom_text(
    data = top5_annotations,
    aes(
      x = 100 * top5_share,
      y = y_position + 0.43,
      label = annotation
    ),
    size = 3.15,
    fontface = "bold",
    vjust = 0,
    color = "black"
  ) +
  scale_fill_manual(
    values = segment_fills,
    breaks = segment_order,
    labels = c(
      "Top-1 exact",
      "Added at ranks 2--3",
      "Added at ranks 4--5",
      "Outside Top-5"
    ),
    drop = FALSE
  ) +
  scale_color_identity() +
  scale_x_continuous(
    name = "Share of evaluated pixels (%)",
    limits = c(
      0,
      100
    ),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    labels = function(x) {
      paste0(
        x,
        "%"
      )
    },
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = method_axis$y_position,
    labels = method_axis$method_label,
    limits = c(
      0.55,
      length(method_order) +
        0.75
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  labs(
    title = "Reference-period ecotype agreement by rank",
    subtitle = paste0(
      "Bars partition pixels by the first rank interval in which the ",
      "reference ecotype is recovered."
    ),
    fill = NULL,
    caption = paste0(
      "Top-1 uses the assigned map and the 1e-4 tie rule. At Top-3 and ",
      "Top-5, a pixel agrees when its reference ecotype occurs within the ",
      "first k dual-suitability ranks. Segment labels are percentages of all ",
      "evaluated pixels."
    )
  ) +
  theme_bw(
    base_size = 10.5
  ) +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(
      color = "#E3E3E3",
      linewidth = 0.35
    ),
    axis.line.x = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      size = 9.3,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 10.2,
      face = "bold",
      color = "black",
      margin = margin(
        r = 8
      )
    ),
    axis.title.x = element_text(
      size = 10.2,
      margin = margin(
        t = 8
      )
    ),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.key.width = grid::unit(
      1.15,
      "cm"
    ),
    legend.key.height = grid::unit(
      0.38,
      "cm"
    ),
    legend.text = element_text(
      size = 9.2
    ),
    plot.title = element_text(
      face = "bold",
      size = 13.2,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 9.8,
      hjust = 0,
      margin = margin(
        b = 6
      )
    ),
    plot.caption = element_text(
      size = 8.4,
      hjust = 0,
      lineheight = 1.06,
      margin = margin(
        t = 8
      )
    ),
    plot.margin = margin(
      10,
      14,
      10,
      10
    )
  )

ggsave(
  filename = output_figure,
  plot = figure_topk_agreement,
  width = 9.2,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# 6. Save figure-source tables ==================================================

plot_data_out <- copy(
  segment_data
)

plot_data_out[
  ,
  `:=`(
    segment = as.character(segment),
    segment_percent = 100 * segment_share,
    cumulative_percent = 100 * x_max
  )
]

setcolorder(
  plot_data_out,
  c(
    "method",
    "method_label",
    "segment",
    "segment_pixels",
    "segment_share",
    "segment_percent",
    "x_min",
    "x_max",
    "cumulative_percent",
    "compared_pixels"
  )
)

fwrite(
  plot_data_out,
  output_plot_data
)

agreement_summary_out <- copy(
  agreement_summary
)

agreement_summary_out[
  ,
  `:=`(
    agreement_percent = 100 * agreement_share,
    rank_label = paste0(
      "Top-",
      rank_cutoff
    )
  )
]

setcolorder(
  agreement_summary_out,
  c(
    "method",
    "method_label",
    "rank_cutoff",
    "rank_label",
    "compared_pixels",
    "agreement_pixels",
    "agreement_share",
    "agreement_percent"
  )
)

fwrite(
  agreement_summary_out,
  output_summary
)


# 7. Console audit ==============================================================

cat(
  "\nCOMPLETE\n",
  "Figure: ",
  output_figure,
  "\nPlot data: ",
  output_plot_data,
  "\nAgreement summary: ",
  output_summary,
  "\n\n",
  sep = ""
)

print(
  agreement_summary_out[
    ,
    .(
      method_label,
      rank_label,
      evaluated_pixels = format_integer(
        compared_pixels
      ),
      agreement_pixels = format_integer(
        agreement_pixels
      ),
      agreement_percent = round(
        agreement_percent,
        3
      )
    )
  ]
)
