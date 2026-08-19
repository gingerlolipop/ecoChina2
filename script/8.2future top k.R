# 8.2 Fast Top-k analysis using existing ranking rasters
# ==============================================================================
# Uses ONLY existing outputs:
#   - dual suit ranking var/{method}/{scenario}/ranked_zone.tif
#   - result maps/{method}/assigned_zone_*.tif
#   - raster/ecosys_ori.tif
#
# It does NOT:
#   - refit RF models
#   - rerun climate/soil prediction
#   - recalculate dual suitability
#   - rank the 53 zones again
#   - recalculate suitable area >= 0.4
#   - recalculate assigned-zone area change already produced elsewhere
#
# Outputs:
#   1. Reference-period Top-1/2/3/5 agreement.
#   2. Future retention of the normal-period assigned zone:
#        stable Top-1
#        changed, but former zone remains in Top-3
#        changed, former zone ranks 4-5
#        changed, former zone falls below Top-5
#        novel
#   3. The same future categories by normal-period assigned zone.
#
# Top-1 uses the existing assigned map, so the original tie-retention rule in
# script 4.1 is preserved. Top-2/3/5 use the saved ranking order from script 8.1.
# ==============================================================================

library(terra)
library(data.table)

rm(list = ls())
gc()


# 0. Paths and settings ==========================================================

base_dir <- "H:/Jing/ecoChina2"

reference_file <- file.path(
  base_dir,
  "raster/ecosys_ori.tif"
)

ranking_root <- file.path(
  base_dir,
  "dual suit ranking var"
)

result_map_root <- file.path(
  base_dir,
  "result maps"
)

out_dir <- file.path(
  base_dir,
  "assessment_var",
  "future_topk_analysis"
)

cache_dir <- file.path(
  out_dir,
  "cache_fast_v3"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  cache_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

zoneID <- c(
  1:7,
  9:50,
  52:55
)

method_order <- c(
  "rf_var",
  "mf_var"
)

method_labels <- c(
  rf_var = "Plain RF",
  mf_var = "Plain MF RF"
)

future_order <- c(
  "2011-2040SSP245",
  "2041-2070SSP245",
  "2071-2100SSP245",
  "2011-2040SSP585",
  "2041-2070SSP585",
  "2071-2100SSP585"
)

scenario_order <- c(
  "normal",
  future_order
)

dual_threshold <- 0.4
tie_tol <- 1e-4
novel_value <- 99L

# Read only the first five layers of ranked_zone.tif.
top_k <- 5L

# Small row blocks keep memory use low.
chunk_rows <- 30L


# 1. Helpers ====================================================================

div <- function(a, b) {
  
  ifelse(
    is.finite(b) & b > 0,
    a / b,
    NA_real_
  )
}


scenario_fields <- function(scenario) {
  
  data.table(
    period = sub(
      "SSP.*$",
      "",
      scenario
    ),
    ssp = sub(
      "^.*(SSP[0-9]+)$",
      "\\1",
      scenario
    )
  )
}


ranked_zone_file <- function(
    method,
    scenario) {
  
  file.path(
    ranking_root,
    method,
    scenario,
    "ranked_zone.tif"
  )
}


assigned_map_file <- function(
    method,
    scenario) {
  
  file.path(
    result_map_root,
    method,
    paste0(
      "assigned_zone_",
      scenario,
      "_threshold",
      dual_threshold,
      "_tol",
      tie_tol,
      "_novel",
      novel_value,
      "_maskNA8_noNovelNormal.tif"
    )
  )
}


cache_is_current <- function(
    cache_files,
    input_files) {
  
  cache_files <- unlist(
    cache_files
  )
  
  input_files <- unlist(
    input_files
  )
  
  if (!all(file.exists(cache_files))) {
    return(FALSE)
  }
  
  if (!all(file.exists(input_files))) {
    return(FALSE)
  }
  
  min(
    file.info(cache_files)$mtime
  ) >=
    max(
      file.info(input_files)$mtime
    )
}


check_raster <- function(
    file,
    template,
    min_layers = 1L) {
  
  if (!file.exists(file)) {
    stop(
      "Missing file: ",
      file
    )
  }
  
  x <- rast(
    file
  )
  
  if (
    nlyr(x) < min_layers
  ) {
    stop(
      "Too few layers in: ",
      file
    )
  }
  
  if (
    !compareGeom(
      x,
      template,
      stopOnError = FALSE
    )
  ) {
    stop(
      "Geometry mismatch: ",
      file
    )
  }
  
  invisible(TRUE)
}


add_area <- function(
    area_vector,
    category,
    keep,
    cell_area,
    n_category) {
  
  if (!any(keep)) {
    return(area_vector)
  }
  
  x <- rowsum(
    cell_area[keep],
    category[keep],
    reorder = FALSE
  )
  
  idx <- as.integer(
    rownames(x)
  )
  
  area_vector[idx] <-
    area_vector[idx] +
    as.numeric(
      x[, 1]
    )
  
  area_vector
}


# 2. Check existing inputs =======================================================

reference_map <- rast(
  reference_file
)

names(reference_map) <-
  "reference_zone"

cell_area <- cellSize(
  reference_map,
  unit = "km"
)

names(cell_area) <-
  "cell_area_km2"

for (method in method_order) {
  
  for (scenario in scenario_order) {
    
    check_raster(
      ranked_zone_file(
        method,
        scenario
      ),
      reference_map,
      min_layers = top_k
    )
    
    check_raster(
      assigned_map_file(
        method,
        scenario
      ),
      reference_map,
      min_layers = 1L
    )
  }
}


# 3. Reference-period Top-k agreement ===========================================

reference_results <- list()

for (method in method_order) {
  
  cache_file <- file.path(
    cache_dir,
    paste0(
      method,
      "_reference_topk.csv"
    )
  )
  
  rank_file <- ranked_zone_file(
    method,
    "normal"
  )
  
  assigned_file <- assigned_map_file(
    method,
    "normal"
  )
  
  if (
    cache_is_current(
      cache_file,
      c(
        reference_file,
        rank_file,
        assigned_file
      )
    )
  ) {
    
    cat(
      "[REUSE] reference Top-k | ",
      method,
      "\n",
      sep = ""
    )
    
    result <- fread(
      cache_file
    )
    
  } else {
    
    cat(
      "[RUN] reference Top-k | ",
      method,
      "\n",
      sep = ""
    )
    
    ranked_zone <- rast(
      rank_file
    )[[1:top_k]]
    
    assigned_map <- rast(
      assigned_file
    )
    
    compared_pixels <- 0
    compared_area <- 0
    
    top1_pixels <- 0
    top2_pixels <- 0
    top3_pixels <- 0
    top5_pixels <- 0
    
    top1_area <- 0
    top2_area <- 0
    top3_area <- 0
    top5_area <- 0
    
    readStart(reference_map)
    readStart(ranked_zone)
    readStart(assigned_map)
    readStart(cell_area)
    
    for (
      row_start in seq(
        1L,
        nrow(reference_map),
        by = chunk_rows
      )
    ) {
      
      nrows_now <- min(
        chunk_rows,
        nrow(reference_map) -
          row_start + 1L
      )
      
      original <- as.integer(
        readValues(
          reference_map,
          row = row_start,
          nrows = nrows_now,
          mat = FALSE
        )
      )
      
      assigned <- as.integer(
        readValues(
          assigned_map,
          row = row_start,
          nrows = nrows_now,
          mat = FALSE
        )
      )
      
      ranks <- readValues(
        ranked_zone,
        row = row_start,
        nrows = nrows_now,
        mat = TRUE
      )
      
      area <- readValues(
        cell_area,
        row = row_start,
        nrows = nrows_now,
        mat = FALSE
      )
      
      if (is.null(dim(ranks))) {
        ranks <- matrix(
          ranks,
          ncol = top_k
        )
      }
      
      valid <- (
        original %in% zoneID &
          assigned %in% zoneID
      )
      
      if (!any(valid)) {
        next
      }
      
      # Top-1 follows script 4.1 exactly, including original-zone tie retention.
      top1 <- (
        valid &
          assigned == original
      )
      
      # Higher-k results use the saved rank order from script 8.1.
      top2 <- (
        valid &
          rowSums(
            ranks[, 1:2, drop = FALSE] ==
              original,
            na.rm = TRUE
          ) > 0
      )
      
      top3 <- (
        valid &
          rowSums(
            ranks[, 1:3, drop = FALSE] ==
              original,
            na.rm = TRUE
          ) > 0
      )
      
      top5 <- (
        valid &
          rowSums(
            ranks[, 1:5, drop = FALSE] ==
              original,
            na.rm = TRUE
          ) > 0
      )
      
      compared_pixels <-
        compared_pixels +
        sum(valid)
      
      compared_area <-
        compared_area +
        sum(
          area[valid],
          na.rm = TRUE
        )
      
      top1_pixels <-
        top1_pixels +
        sum(top1)
      
      top2_pixels <-
        top2_pixels +
        sum(top2)
      
      top3_pixels <-
        top3_pixels +
        sum(top3)
      
      top5_pixels <-
        top5_pixels +
        sum(top5)
      
      top1_area <-
        top1_area +
        sum(
          area[top1],
          na.rm = TRUE
        )
      
      top2_area <-
        top2_area +
        sum(
          area[top2],
          na.rm = TRUE
        )
      
      top3_area <-
        top3_area +
        sum(
          area[top3],
          na.rm = TRUE
        )
      
      top5_area <-
        top5_area +
        sum(
          area[top5],
          na.rm = TRUE
        )
    }
    
    readStop(reference_map)
    readStop(ranked_zone)
    readStop(assigned_map)
    readStop(cell_area)
    
    result <- data.table(
      method = method,
      method_label =
        method_labels[[method]],
      rank_cutoff = c(
        1L,
        2L,
        3L,
        5L
      ),
      compared_pixels =
        compared_pixels,
      matched_pixels = c(
        top1_pixels,
        top2_pixels,
        top3_pixels,
        top5_pixels
      ),
      pixel_share = c(
        div(
          top1_pixels,
          compared_pixels
        ),
        div(
          top2_pixels,
          compared_pixels
        ),
        div(
          top3_pixels,
          compared_pixels
        ),
        div(
          top5_pixels,
          compared_pixels
        )
      ),
      compared_area_km2 =
        compared_area,
      matched_area_km2 = c(
        top1_area,
        top2_area,
        top3_area,
        top5_area
      ),
      area_share = c(
        div(
          top1_area,
          compared_area
        ),
        div(
          top2_area,
          compared_area
        ),
        div(
          top3_area,
          compared_area
        ),
        div(
          top5_area,
          compared_area
        )
      )
    )
    
    fwrite(
      result,
      cache_file
    )
    
    rm(
      ranked_zone,
      assigned_map
    )
    
    gc()
  }
  
  reference_results[[
    length(reference_results) + 1L
  ]] <- result
}


# 4. Future retention of normal-period assigned zone ============================

category_levels <- c(
  "stable_top1",
  "changed_former_in_top3",
  "changed_former_rank4_5",
  "changed_former_below_top5",
  "novel"
)

n_category <- length(
  category_levels
)

future_results <- list()
future_zone_results <- list()

for (method in method_order) {
  
  normal_file <- assigned_map_file(
    method,
    "normal"
  )
  
  normal_map <- rast(
    normal_file
  )
  
  for (scenario in future_order) {
    
    cache_overall <- file.path(
      cache_dir,
      paste0(
        method,
        "_",
        scenario,
        "_future_topk_overall.csv"
      )
    )
    
    cache_by_zone <- file.path(
      cache_dir,
      paste0(
        method,
        "_",
        scenario,
        "_future_topk_by_zone.csv"
      )
    )
    
    rank_file <- ranked_zone_file(
      method,
      scenario
    )
    
    future_file <- assigned_map_file(
      method,
      scenario
    )
    
    if (
      cache_is_current(
        c(
          cache_overall,
          cache_by_zone
        ),
        c(
          normal_file,
          rank_file,
          future_file
        )
      )
    ) {
      
      cat(
        "[REUSE] future Top-k | ",
        method,
        " | ",
        scenario,
        "\n",
        sep = ""
      )
      
      overall_result <- fread(
        cache_overall
      )
      
      zone_result <- fread(
        cache_by_zone
      )
      
    } else {
      
      cat(
        "[RUN] future Top-k | ",
        method,
        " | ",
        scenario,
        "\n",
        sep = ""
      )
      
      future_map <- rast(
        future_file
      )
      
      ranked_zone <- rast(
        rank_file
      )[[1:top_k]]
      
      category_pixels <- numeric(
        n_category
      )
      
      category_area <- numeric(
        n_category
      )
      
      zone_category_pixels <- numeric(
        length(zoneID) *
          n_category
      )
      
      zone_category_area <- numeric(
        length(zoneID) *
          n_category
      )
      
      readStart(normal_map)
      readStart(future_map)
      readStart(ranked_zone)
      readStart(cell_area)
      
      for (
        row_start in seq(
          1L,
          nrow(reference_map),
          by = chunk_rows
        )
      ) {
        
        nrows_now <- min(
          chunk_rows,
          nrow(reference_map) -
            row_start + 1L
        )
        
        normal_zone <- as.integer(
          readValues(
            normal_map,
            row = row_start,
            nrows = nrows_now,
            mat = FALSE
          )
        )
        
        future_zone <- as.integer(
          readValues(
            future_map,
            row = row_start,
            nrows = nrows_now,
            mat = FALSE
          )
        )
        
        ranks <- readValues(
          ranked_zone,
          row = row_start,
          nrows = nrows_now,
          mat = TRUE
        )
        
        area <- readValues(
          cell_area,
          row = row_start,
          nrows = nrows_now,
          mat = FALSE
        )
        
        if (is.null(dim(ranks))) {
          ranks <- matrix(
            ranks,
            ncol = top_k
          )
        }
        
        valid <- (
          normal_zone %in% zoneID &
            future_zone %in%
            c(
              zoneID,
              novel_value
            )
        )
        
        if (!any(valid)) {
          next
        }
        
        stable <- (
          valid &
            future_zone ==
            normal_zone
        )
        
        novel <- (
          valid &
            future_zone ==
            novel_value
        )
        
        changed <- (
          valid &
            !stable &
            !novel
        )
        
        former_top3 <- rep(
          FALSE,
          length(normal_zone)
        )
        
        former_top5 <- rep(
          FALSE,
          length(normal_zone)
        )
        
        changed_rows <- which(
          changed
        )
        
        if (length(changed_rows) > 0) {
          
          former_top3[
            changed_rows
          ] <- rowSums(
            ranks[
              changed_rows,
              1:3,
              drop = FALSE
            ] ==
              normal_zone[
                changed_rows
              ],
            na.rm = TRUE
          ) > 0
          
          former_top5[
            changed_rows
          ] <- rowSums(
            ranks[
              changed_rows,
              1:5,
              drop = FALSE
            ] ==
              normal_zone[
                changed_rows
              ],
            na.rm = TRUE
          ) > 0
        }
        
        category <- rep(
          NA_integer_,
          length(normal_zone)
        )
        
        category[stable] <- 1L
        
        category[
          changed &
            former_top3
        ] <- 2L
        
        category[
          changed &
            !former_top3 &
            former_top5
        ] <- 3L
        
        category[
          changed &
            !former_top5
        ] <- 4L
        
        category[novel] <- 5L
        
        keep <- (
          valid &
            !is.na(category)
        )
        
        category_pixels <-
          category_pixels +
          tabulate(
            category[keep],
            nbins = n_category
          )
        
        category_area <-
          add_area(
            category_area,
            category,
            keep,
            area,
            n_category
          )
        
        
        # By normal-period assigned zone -----------------------------------------
        
        normal_index <- match(
          normal_zone[keep],
          zoneID
        )
        
        code <- (
          normal_index - 1L
        ) *
          n_category +
          category[keep]
        
        zone_category_pixels <-
          zone_category_pixels +
          tabulate(
            code,
            nbins =
              length(zoneID) *
              n_category
          )
        
        zone_area_sum <- rowsum(
          area[keep],
          code,
          reorder = FALSE
        )
        
        idx <- as.integer(
          rownames(
            zone_area_sum
          )
        )
        
        zone_category_area[idx] <-
          zone_category_area[idx] +
          as.numeric(
            zone_area_sum[, 1]
          )
      }
      
      readStop(normal_map)
      readStop(future_map)
      readStop(ranked_zone)
      readStop(cell_area)
      
      fields <- scenario_fields(
        scenario
      )
      
      total_pixels <- sum(
        category_pixels
      )
      
      total_area <- sum(
        category_area
      )
      
      changed_pixels <- sum(
        category_pixels[2:4]
      )
      
      changed_area <- sum(
        category_area[2:4]
      )
      
      overall_result <- data.table(
        method = method,
        method_label =
          method_labels[[method]],
        scenario = scenario,
        period = fields$period,
        ssp = fields$ssp,
        
        common_valid_pixels =
          total_pixels,
        common_valid_area_km2 =
          total_area,
        
        stable_top1_pixels =
          category_pixels[1],
        stable_top1_pixel_share =
          div(
            category_pixels[1],
            total_pixels
          ),
        stable_top1_area_km2 =
          category_area[1],
        stable_top1_area_share =
          div(
            category_area[1],
            total_area
          ),
        
        changed_existing_pixels =
          changed_pixels,
        changed_existing_pixel_share =
          div(
            changed_pixels,
            total_pixels
          ),
        changed_existing_area_km2 =
          changed_area,
        changed_existing_area_share =
          div(
            changed_area,
            total_area
          ),
        
        novel_pixels =
          category_pixels[5],
        novel_pixel_share =
          div(
            category_pixels[5],
            total_pixels
          ),
        novel_area_km2 =
          category_area[5],
        novel_area_share =
          div(
            category_area[5],
            total_area
          ),
        
        former_top3_among_changed_pixel =
          div(
            category_pixels[2],
            changed_pixels
          ),
        former_top3_among_changed_area =
          div(
            category_area[2],
            changed_area
          ),
        
        former_top5_among_changed_pixel =
          div(
            category_pixels[2] +
              category_pixels[3],
            changed_pixels
          ),
        former_top5_among_changed_area =
          div(
            category_area[2] +
              category_area[3],
            changed_area
          ),
        
        former_below_top5_among_changed_pixel =
          div(
            category_pixels[4],
            changed_pixels
          ),
        former_below_top5_among_changed_area =
          div(
            category_area[4],
            changed_area
          )
      )
      
      
      zone_grid <- CJ(
        normal_assigned_zone =
          zoneID,
        category_index =
          seq_len(
            n_category
          ),
        sorted = FALSE
      )
      
      zone_grid[
        ,
        category :=
          category_levels[
            category_index
          ]
      ]
      
      zone_grid[
        ,
        `:=`(
          pixel_count =
            zone_category_pixels,
          area_km2 =
            zone_category_area
        )
      ]
      
      zone_grid[
        ,
        `:=`(
          pixel_share =
            pixel_count /
            sum(pixel_count),
          area_share =
            area_km2 /
            sum(area_km2)
        ),
        by = normal_assigned_zone
      ]
      
      zone_grid[
        ,
        `:=`(
          method = method,
          method_label =
            method_labels[[method]],
          scenario = scenario,
          period = fields$period,
          ssp = fields$ssp
        )
      ]
      
      zone_result <- zone_grid[
        ,
        .(
          method,
          method_label,
          scenario,
          period,
          ssp,
          normal_assigned_zone,
          category,
          pixel_count,
          pixel_share,
          area_km2,
          area_share
        )
      ]
      
      fwrite(
        overall_result,
        cache_overall
      )
      
      fwrite(
        zone_result,
        cache_by_zone
      )
      
      rm(
        future_map,
        ranked_zone
      )
      
      gc()
    }
    
    future_results[[
      length(future_results) + 1L
    ]] <- overall_result
    
    future_zone_results[[
      length(future_zone_results) + 1L
    ]] <- zone_result
  }
  
  rm(
    normal_map
  )
  
  gc()
}


# 5. Save final tables ===========================================================

reference_topk <- rbindlist(
  reference_results,
  fill = TRUE
)

future_topk <- rbindlist(
  future_results,
  fill = TRUE
)

future_topk_by_zone <- rbindlist(
  future_zone_results,
  fill = TRUE
)

setorder(
  reference_topk,
  method,
  rank_cutoff
)

setorder(
  future_topk,
  method,
  scenario
)

setorder(
  future_topk_by_zone,
  method,
  scenario,
  normal_assigned_zone
)

fwrite(
  reference_topk,
  file.path(
    out_dir,
    "reference_map_topk_agreement.csv"
  )
)

fwrite(
  future_topk,
  file.path(
    out_dir,
    "future_topk_retention_overall.csv"
  )
)

fwrite(
  future_topk_by_zone,
  file.path(
    out_dir,
    "future_topk_retention_by_normal_assigned_zone.csv"
  )
)


# 6. Console summary =============================================================

cat(
  "\n[REFERENCE TOP-K]\n"
)

print(
  reference_topk[
    ,
    .(
      method_label,
      rank_cutoff,
      pixel_pct =
        round(
          100 * pixel_share,
          2
        ),
      area_pct =
        round(
          100 * area_share,
          2
        )
    )
  ]
)

cat(
  "\n[FUTURE TOP-K RETENTION]\n"
)

print(
  future_topk[
    ,
    .(
      method_label,
      scenario,
      stable_top1_pct =
        round(
          100 *
            stable_top1_area_share,
          2
        ),
      changed_existing_pct =
        round(
          100 *
            changed_existing_area_share,
          2
        ),
      novel_pct =
        round(
          100 *
            novel_area_share,
          2
        ),
      former_top3_among_changed_pct =
        round(
          100 *
            former_top3_among_changed_area,
          2
        ),
      former_top5_among_changed_pct =
        round(
          100 *
            former_top5_among_changed_area,
          2
        )
    )
  ]
)

cat(
  "\nCOMPLETE\n",
  "Output: ",
  out_dir,
  "\n",
  sep = ""
)