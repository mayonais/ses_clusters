# Los Angeles heat and ED analysis

Run the R scripts from the project root. This project contains SES clustering,
daily ED-count reconstruction, and Poisson heat-response models.

## Current inputs and results

The active folder is:

```text
Los_Angeles_discrete_none/
  Los Angeles_ACS_zcta_clustered_for_ED.rds
  Los Angeles_daily_ed_no_heat.rds
  Los Angeles_ACS_diss_mats.rds
  Los Angeles_ACS_diss_mats_MIN.rds
  create_map.R
  linear_cluster/
  climate_zone/
  climate_zone_added/
  model_cache/
```

The cluster assignments are the 20-covariate, 14-cluster solution used by all
current heat models. The prepared daily ED input is unchanged from the file
used for the saved fits.

| Results folder | Model |
| --- | --- |
| `linear_cluster/` | Original linear acute-heat effects with SES-cluster heat slopes; no climate-zone heat curves |
| `climate_zone/` | Shared RW2 curve and pooled climate-zone RW2 deviations; no explicit common or SES-cluster linear acute-heat slopes |
| `climate_zone_added/` | Shared and climate-zone RW2 curves plus SES-cluster heat-slope interactions |
| `model_cache/` | Compatible cached fits for these model variants |

Each model folder contains its saved fits, coefficient summaries, WAIC/DIC TXT
summary, AR summaries/posterior draws, and applicable acute-heat tables and plot
objects. Existing model filenames retain the legacy
`Los_Angeles_discrete_none_HEAT_by_cluster` prefix; this is a filename prefix,
not the location of an active folder.

The original model's saved Census geometry and forest plots are in
`linear_cluster/acute_heat_poisson/` and `linear_cluster/acute_heat_nb/`.
Other variants use that saved geometry.

## Scripts

- `hri_ed_analysis.R`: fitting and model summaries; currently selects the
  RW2-only acute-heat `climate_zone` Poisson model. Negative-binomial refitting is disabled.
- `acute_heat.R`: redraws acute-heat maps from the saved fits.
- `acute_heat_irr_maps.R`: computes complete one-SD heat contrasts from INLA
  linear-combination marginals and exports the matching IRR/probability maps.
- `shared_functions.R`: data preparation, fitting, contrasts, and map helpers.
- `cluster_ar.R`: signed cluster-average excess-rate inference.
- `create_cluster.R`: upstream 20-covariate SES clustering; reads/writes the
  canonical output folder. It retains 2,000 burn-in iterations and 10,000 sweeps.
- `suppressed_data.R`: daily-count reconstruction, now using the canonical
  cluster/input folder. Its allocation calculations are unchanged.
- `cluster_diagnostics.R`: preserves the existing historical MIN-chain
  diagnostics using the explicitly named `Los Angeles_ACS_diss_mats_MIN.rds`.
  Those 14-covariate chains produced a 9-cluster solution and are different
  from the 20-covariate, 14-cluster assignments used by the heat models.

Redraw maps without starting clustering or refitting:

```r
source("acute_heat.R")
source("Los_Angeles_discrete_none/create_map.R")
```

Choose the desired `model_variant` in the map scripts. Both default to
`climate_zone`; `model_family` defaults to Poisson.
Fitting, count reconstruction, and clustering are separate actions; the map
scripts do not automatically run those upstream steps.

The separate simulation workflow reads the canonical assignments and writes
its own outputs into `CA_simulated_ed/results/`.

## Maps and interpretation

Map PNGs and PDFs remain under matching folders in `maps/`:
`linear_cluster/`, `climate_zone/`, and `climate_zone_added/`.
Geographic maps use a 10-by-8-inch canvas; PNGs are 300 dpi. Gray areas without model results have diagonal hatching. Additional IRR point-and-interval plots show climate-zone effects, and cluster-by-zone effects for `climate_zone_added`. AR point-and-interval plots have been removed; AR probability maps and saved posterior draws remain.

For climate-zone models, acute heat is each day's deviation from the ZIP's
average, expressed in the county-wide heat SD. The averages and SD use retained
May–October 2018 ZIP-day observations, without population weighting.
Zero heat-index values are treated as unavailable; genuine zero ED counts are
retained.

- IRRs compare a ZIP's average heat with one SD above that average, including
  both RW2 curve differences and, for `climate_zone_added`, the SES-cluster linear heat slope.
  INLA linear-combination marginals provide medians and credible intervals.
- The acute-heat probability map shows `Pr(IRR_1SD > 1 | Y)`, labeled 0%–100%.
- AR compares observed heat with each SES cluster's month-specific 2018 mean heat, expressed on the ZIP-centered acute-heat scale.
  It averages signed rate differences equally across ZIP-days within each SES
  cluster. Its map shows the probability that the cluster-average excess rate
  is positive, estimated from 2,000 joint posterior draws, with periodic
  garbage collection.

Climate-zone RW2 deviations share a precision. Only sum-to-zero constraints are
retained; `extraconstr` is removed and RW2 linear trends are included in the
complete contrasts. ZIPs 90010, 90704, and 91046 have no zone assignment and use
the shared curve without a zone deviation. Population exposure remains in the
model.

Concise mathematics and methods are in
`methods/climate_zone_addition/main.tex`.
Other files in `methods/` are earlier methods/supplement drafts and may refer to
older map locations or model versions.

## Archive and saved clustering runs

The reorganization on 2026-10-27 moved unused results into
`archive/output_reorganization_2026-10-27/`, including:

- the alternative `Los_Angeles_discrete_none_MIN` assignments and results;
- older August heat results from the canonical folder;
- superseded climate-model caches, unused unpooled NB outputs, and the replaced
  IRR sampling cache;
- intermediate diagnostic logs and the remaining former HEAT folder contents.

No files were deleted. `move_manifest.csv` inside that archive records original
and destination paths and SHA-256 hashes. The small MIN dissimilarity/run-metadata
file is also retained in the active folder for the existing diagnostics.

Raw saved MIN MCMC traces remain unchanged in
`archive/cleanup_2026-10-15/`. Those available traces describe 14 covariates;
do not relabel them as the 20-covariate chains or silently combine the two.
The canonical 20-covariate saved dissimilarity objects are retained separately.

## Dependencies and Git

Heat scripts use `INLA`, `dplyr`, `ggplot2`, `lubridate`, `sf`,
`tigris`, and `scales`. Upstream clustering and diagnostics require their
additional listed packages. Scripts do not install packages automatically.

Large fitted models, caches, raw MCMC traces, raw datasets, and local archive
contents are excluded from Git. Active canonical clustering inputs and summary
files are no longer excluded merely because of the canonical folder name.
Fitted models must be regenerated locally; raw datasets and `datasets/cz_zip.csv` are local inputs. Diagnostic scripts, diagnostic outputs, logs, and caches are excluded from publication.
