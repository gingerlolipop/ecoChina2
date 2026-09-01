# ecoChina2

Reproducible ecological-niche workflow for China's 53 modelled vegetation zones. The main analysis uses only the selected-variable **Multi-Forest** one-vs-rest models; the former plain-RF and other superseded workflows are not part of this branch. A separate class-balanced multiclass random forest is retained as a robustness analysis.

## Analysis scripts

| Script | Purpose |
|---|---|
| `1. data.R` | Prepare the reference vegetation raster, zone lookup, cell coordinates and reference-period climate table. |
| `1.2 species to pop2.R` | Define tree populations as retained species-by-reference-zone combinations (at least 10 occupied cells). |
| `2.15 climate onehot rf.R` | Split climate data, run backward purging, select exactly nine climate predictors per zone, then fit the climate Multi-Forests. |
| `2.4 soil onehot rf.R` | Split soil data, run backward purging, select exactly nine soil predictors per zone, then fit the soil Multi-Forests. |
| `3.2 one hot clim and soil prediction and dual suitability.R` | Project climate and static-soil suitability for the reference period and six future scenarios; calculate primary dual suitability and the soil-gate sensitivity analysis. |
| `4. assign zone.R` | Assign the best-supported zone per pixel, including the future-only novel-zone rule and reproducible tie handling. |
| `5. assessment.R` | Calculate held-out binary-model metrics and reference-map reconstruction metrics. |
| `5.3 feature importance.R` | Compare binary Multi-Forest and multiclass importance, aggregate normalized importance by predictor group and vegetation category, and write the core feature-importance figures. |
| `6. future pop & species niche.R` | Summarize projected zone, population and species niche areas from dual suitability. |
| `7. multiclassrf.R` | Fit and assess the class-balanced multiclass RF and project it to the reference period and all six future scenarios. |
| `8. future top k.R` | Rank dual suitability for `k = 1, 2, 3, 4, 5`; calculate agreement, retention, matched-change, analogue, population and species diagnostics. |
| `11. visualization.R` | Rebuild manuscript-ready figures from the preceding outputs, including feature-importance, sensitivity, multiclass and Top-k figures. |
| `color_palette.R` | Create the canonical zone palette used by assessments and figures. |

Shared functions are in `functions/`.

## Inputs

Large input rasters and occurrence data are not stored in Git. From the repository root, place inputs in the following default layout or set the environment variables shown below.

```text
raster/veg_3                         # or raster/veg_3.tif
data raw/1. coord.csv
data raw/new_soil_raster.csv
species data/<species_code>_coord.csv
data/rasters/soil/<soil_variable>.tif
data/rasters/climate/
  Normal_1961_1990/<climate_variable>.tif
  8GCMs_ensemble_ssp245_2011-2040/<climate_variable>.tif
  8GCMs_ensemble_ssp245_2041-2070/<climate_variable>.tif
  8GCMs_ensemble_ssp245_2071-2100/<climate_variable>.tif
  8GCMs_ensemble_ssp585_2011-2040/<climate_variable>.tif
  8GCMs_ensemble_ssp585_2041-2070/<climate_variable>.tif
  8GCMs_ensemble_ssp585_2071-2100/<climate_variable>.tif
```

An existing `data raw/1. zoneID_Clm_800m_Normal_1961_1990SY.csv` is reused directly; `data/climate_reference.csv` is also accepted as a fallback instead of extracting it again. The temporary `data/raw/` layout used by an earlier clean draft remains accepted where the scripts list it as a fallback. The soil table must contain `zoneID`, `x`, `y` and the 15 topsoil predictors named in script 2.4. Species files must use the codes listed in script 1.2 and contain longitude, latitude and presence/status fields.

| Environment variable | Use |
|---|---|
| `ECOCHINA2_DIR` | Repository root when scripts are not launched from it. |
| `ECOCHINA2_CLIMATE_DIR` | Climate-raster root. |
| `ECOCHINA2_SOIL_TABLE` | Reference soil CSV. |
| `ECOCHINA2_SOIL_RASTER_DIR` | Soil-raster root used for projection. |
| `ECOCHINA2_SPECIES_DIR` | Species-occurrence directory. |
| `ECOCHINA2_FORCE_DATA` | Set to `true` only to rebuild prepared reference data. |
| `ECOCHINA2_FORCE_POPULATIONS` | Set to `true` only to rebuild the species-by-zone population lookup. |
| `ECOCHINA2_FORCE_SPLIT` | Set to `true` only to replace the saved 70:30 splits; dependent selection and final training are then rebuilt. |
| `ECOCHINA2_FORCE_SELECTION` | Set to `true` only to rerun backward purging; dependent final training is then rebuilt. |
| `ECOCHINA2_FORCE_CLIMATE_TRAINING` | Set to `true` only to refit the 53 climate Multi-Forests. |
| `ECOCHINA2_FORCE_SOIL_TRAINING` | Set to `true` only to refit the 53 soil Multi-Forests. |
| `ECOCHINA2_REUSE_SUITABILITY` | Reuse valid climate and soil suitability rasters; defaults to `true`. |
| `ECOCHINA2_FORCE_SUITABILITY` | Set to `true` only to rerun climate/soil prediction, primary dual suitability and the dependent gate-sensitivity table. |
| `ECOCHINA2_REUSE_DUAL` | Reuse valid primary dual-suitability rasters; defaults to `true`. |
| `ECOCHINA2_FORCE_DUAL` | Set to `true` only to rebuild primary dual-suitability rasters. |
| `ECOCHINA2_FORCE_SENSITIVITY` | Set to `true` only to recalculate the soil-gate sensitivity table without rebuilding suitability rasters. |
| `ECOCHINA2_WRITE_GATE_RASTERS` | Set to `true` to materialize every alternative-gate raster; complete stacks newer than their suitability inputs are reused. Defaults to `false`. |
| `ECOCHINA2_REUSE_MAPS` | Reuse valid assigned-zone rasters; defaults to `true`. |
| `ECOCHINA2_FORCE_MAPS` | Set to `true` only to rebuild assigned-zone rasters. |
| `ECOCHINA2_FORCE_ASSESSMENT` | Set to `true` only to rebuild assessment tables; otherwise input modification times are checked. |
| `ECOCHINA2_REUSE_POPULATION` | Reuse valid population/species raster outputs; defaults to `true`. |
| `ECOCHINA2_FORCE_MULTICLASS` | Set to `true` only to refit and reproject the multiclass RF. |
| `ECOCHINA2_FORCE_RANKING` | Set to `true` only after the 53 dual-suitability surfaces themselves have changed; defaults to reusing valid rank stacks. |
| `ECOCHINA2_FIGURE_SECTIONS` | Comma-separated figure families: `assessment`, `importance`, `maps`, `niche`, `topk`, or `all`. |
| `ECOCHINA2_STACK_SCENARIOS` | Comma-separated scenarios for the detailed 5×2 and exploded Top-k figures; defaults to `normal,2071-2100SSP585`. |

All rasters used together must have compatible coordinate reference systems, extent and resolution. Run the workflow from the repository root even when `ECOCHINA2_DIR` is set.
Top-k rank reuse is automatic: script 8 accepts valid existing rank stacks with at least five layers and does not overwrite them.

## Software

The scripts use R (4.x recommended) and the packages `terra`, `data.table`, `randomForest`, `CEMT`, `doSNOW`, `foreach`, `ggplot2`, `patchwork`, `scales`, `circlize` and `jsonlite`. Install the packages available from your configured R repository before running; `CEMT` must also be installed from the same package source/version used for the analysis.

`CEMT` was installed from an author-held local source archive in the original workflow; that archive and its version metadata are not currently in this repository or on CRAN. Before reviewer release, add the exact distributable archive (or a permanent source/DOI plus version) to the software record. Without it, a completely cold-environment rerun of scripts 2.15 and 2.4 cannot be guaranteed, although existing fitted-model checkpoints remain reusable.

```r
install.packages(c(
  "terra", "data.table", "randomForest", "doSNOW", "foreach",
  "ggplot2", "patchwork", "scales", "circlize", "jsonlite"
))
```

## Run order

Run each command from the repository root. Climate and soil training (steps 3 and 4 below) are independent and may be run in parallel if resources allow.

```sh
Rscript "script/1. data.R"
Rscript "script/color_palette.R"
Rscript "script/1.2 species to pop2.R"
Rscript "script/2.15 climate onehot rf.R"
Rscript "script/2.4 soil onehot rf.R"
Rscript "script/3.2 one hot clim and soil prediction and dual suitability.R"
Rscript "script/4. assign zone.R"
Rscript "script/5. assessment.R"
Rscript "script/6. future pop & species niche.R"
Rscript "script/7. multiclassrf.R"
Rscript "script/5.3 feature importance.R"
Rscript "script/8. future top k.R"
Rscript "script/11. visualization.R"
```

The numerical exception in this order is intentional: script 7 precedes 5.3 because the consensus feature-importance analysis needs the multiclass model. Scripts use fixed seeds for splitting, sampling and tie resolution.

### Restart boundaries

Heavy artifacts are restart checkpoints. Existing valid models and rasters are reused by default at their established paths. Do not delete or move them merely to use this branch.

- To extend the existing analysis to `k = 1:5`, run only:

  ```sh
  Rscript "script/8. future top k.R"
  ECOCHINA2_FIGURE_SECTIONS=topk Rscript "script/11. visualization.R"
  ```

  On Windows Command Prompt, set the selector first with `set ECOCHINA2_FIGURE_SECTIONS=topk`. Script 8 reads the first five layers of the existing ten-layer rank rasters and calculates the missing `k = 2` and `k = 4` products. It does not refit RFs, predict suitability, rebuild dual suitability, reassign Top-1 maps, or rerank the 53-zone raster stack.

  Top-k area summaries are cached per scenario. After this first extension, a changed rank stack invalidates only its own scenario rather than forcing the other six scenarios through the raster-area calculations again.

- For figure-only revisions, run script 11. Use `ECOCHINA2_FIGURE_SECTIONS` to avoid reopening unrelated rasters.
- After an assignment-rule change, run script 4 with `ECOCHINA2_FORCE_MAPS=true`, followed by 5, 6, 7, 8 and 11. Script 7 reuses its fitted model and valid projection maps while refreshing lightweight comparisons; script 8 reuses the unchanged ranks.
- After changing a climate or soil model, start at the corresponding 2.x script, then run script 3.2 with `ECOCHINA2_FORCE_SUITABILITY=true`. Continue with scripts 4–7, and set `ECOCHINA2_FORCE_RANKING=true` for script 8 because the underlying dual-suitability surfaces changed. Run script 11 last.
- Force flags are deliberately opt-in. Set only the flag for the stage that must be rebuilt; valid outputs from all other stages remain untouched.

## Model and decision definitions

- Modelled zones are `1:7`, `9:50` and `52:55` (53 zones); other reference-raster values are excluded from binary fitting, assigned-map outputs and their assessment masks.
- Each binary climate or soil model is a one-vs-rest Multi-Forest. The data are split 70:30, backward purging is performed on the training portion, exactly nine predictors are retained per zone, and the final object combines 10 resampled forests of 100 trees each. Held-out binary classification uses a probability threshold of `0.50`.
- The primary dual suitability is climate suitability where soil suitability is greater than `0.20`, and zero otherwise. Sensitivity gates are `0`, `0.10`, `0.15`, `0.20`, `0.25`, `0.30`, `0.40` and `0.50`; gate `0` is climate-only.
- A future pixel is assigned Zone 99 when its maximum dual suitability is below `0.40`. The novel-zone rule is not applied to the reference map. Values within `1e-4` of the maximum are tied; the observed reference zone is retained when possible, followed by seeded random resolution.
- Population and species suitable-area summaries use dual suitability greater than or equal to `0.40`.
- The multiclass robustness model uses 30 backward-selected climate/soil predictors, class-balanced training and 500 trees. It is projected to the reference period plus both SSPs for all three future periods. Multiclass projections do not use the binary `0.40` rule and never assign Zone 99.
- Top-k outputs retain ranks 1 through 5 and their dual-suitability values. For exact compatibility with the established ten-layer rank stacks, a forced ranking rebuild uses the original domain (`reference != NA` and `reference != 8`) and ranks only strictly positive dual suitability. The novel-area definition remains the Top-1 `0.40` rule and is therefore invariant to `k`.

## Established checkpoints and generated outputs

The older names below are intentionally retained for computationally expensive artifacts. The earlier request for shorter folder names applies to new, inexpensive tables and figures, not to model or raster caches that may take hours or days to reproduce.

| Directory | Contents |
|---|---|
| `results/` | Saved train/test splits and climate/soil backward-purging paths. |
| `rf/` | Climate models `clm_mfVar_zone*.Rdata`. |
| `rf_soil/` | Soil models `soil_mf_zone*.Rdata`. |
| `clim suitability/mf_var/` | Reference and future climate-suitability rasters. |
| `soil suitability/plain_mf/normal/` | Static soil-suitability rasters used by `mf_var`. `plain_mf` is the established cache label, not a retained ordinary-RF workflow. |
| `dual suit/mf_var/` | Primary dual-suitability rasters for all seven scenarios. |
| `result maps/mf_var/` | Assigned binary Multi-Forest maps, retaining the established threshold/tolerance filenames. |
| `dual suit ranking var/mf_var/` | Existing ten-layer zone and suitability ranks. Scripts 8 and 11 read only ranks 1–5. |
| `future tree niche var/` | Established assigned-map population/species niche rasters and area tables. |
| `future tree niche dual suitability var/` | Reusable population/species niche and population-ranking rasters. |
| `rf_multiclass/`, `tmp_multiclass_rf/` | Multiclass model and reusable train/test or alignment caches. |
| `result maps/multiclass_rf/` | Reference and six future multiclass projection maps. Complete, current class-probability stacks are indexed when available; a valid legacy map is not recomputed solely to recreate an optional probability stack. |
| `assessment/` | Short, rebuildable assessment, feature-importance and multiclass tables. New Top-k tables are in `assessment/topk/`; legacy Top-k caches remain readable in `assessment_var/`. |
| `figures/` | Manuscript and supplementary figures generated by script 11. |

These generated directories are ignored by Git. A valid existing checkpoint wins over a fallback clean-branch path. Reproduction is complete when the required scripts finish without errors and script 11 recreates the requested figure sections from the saved tables and rasters.
