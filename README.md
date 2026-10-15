# CalHeatScore

Los Angeles ZCTA analyses of socioeconomic clusters, daily emergency department
(ED) counts, and acute heat. The current heat workflow fits Poisson and
negative-binomial models and produces cluster-specific maps and summaries.

## Run the acute-heat figures

Open `calheatscore.Rproj` in RStudio. At the top of `acute_heat.R`, set:

```r
model_family <- "both"  # "poisson", "nbinomial", or "both"
```

Then run the whole script from the project:

```r
source("acute_heat.R")
```

This regenerates the probability map, incidence rate ratio (IRR) map, and IRR
forest plot for the selected models. Each figure is saved as PNG and PDF.
Compatible fitted models are reused automatically; set `force_refit <- TRUE`
in the script only when a fresh fit is wanted. Running the script overwrites
the selected models' figure and summary files.

To redraw only the IRR maps for both models from existing summaries:

```r
source("acute_heat_irr_maps.R")
```

### Required inputs and packages

The prepared inputs are in `Los_Angeles_discrete_none_HEAT_by_cluster/`:

- `Los Angeles_ACS_zcta_clustered_for_ED.rds`: cluster assignments and area characteristics.
- `Los Angeles_daily_ed_no_heat.rds`: daily ED counts, population, dates, and heat measurements.

The acute-heat scripts require R packages `INLA`, `dplyr`, `ggplot2`, `lubridate`,
`sf`, `tigris`, and `scales`. The scripts do not install packages. If saved
geometry is unavailable, `tigris` downloads Census boundaries and internet access
is needed. The saved geometry uses 2020 ZCTAs.

### Where results go

Within `Los_Angeles_discrete_none_HEAT_by_cluster/`:

| Folder | Contents |
| --- | --- |
| `acute_heat_poisson/` | Poisson figures, summary tables, and map geometry |
| `acute_heat_nb/` | Negative-binomial figures, summary tables, and map geometry |
| `model_cache/` | Shared fitted models and compatibility signatures |

Each model's figure folder contains:

- `acute_heat_probability_map_zcta.png` / `.pdf`
- `acute_heat_irr_map.png` / `.pdf`
- `acute_heat_relative_risk.png` / `.pdf` (the IRR forest plot)
- `acute_heat_cluster_summary.csv`
- `acute_heat_zcta_probabilities.csv`
- Saved plot objects and geometry in `.rds` files.

### Interpretation

Heat is standardized using the mean and SD across all available LA ZCTA-day
heat observations, pooled across clusters without population weighting.
Chronic heat is each area's mean standardized heat; acute heat is the daily
departure from that mean. A one-unit acute increase uses the same pooled heat
SD for every cluster, holding chronic heat fixed.

- The probability map shows `P(beta_acute,k > 0 | Y)` directly from the posterior
  slope marginal, on a fixed linear 0–1 scale.
- The IRR map shows the posterior median of `exp(beta_acute,k)` and equal-tail
  95% credible intervals. IRR = 1 indicates no acute-heat association.
- `P(beta_acute,k > 0 | Y)` equals `P(IRR_k > 1 | Y)`.

## Project files

| File or folder | Role |
| --- | --- |
| `acute_heat.R` | Main entry point for current acute-heat figures |
| `acute_heat_irr_maps.R` | IRR maps from saved model summaries; called by the main script |
| `analysis_functions.R` | Shared data preparation, model fitting/cache, and posterior summaries |
| `hri_ed_analysis.R` | Broader Poisson/negative-binomial ED analysis and cluster summaries |
| `cluster_ar.R` | Posterior excess-rate calculations used by `hri_ed_analysis.R` |
| `create_cluster.R` | Upstream ACS profile clustering |
| `suppressed_data.R` | Upstream reconstruction of suppressed daily ED counts |
| `create_cluster outputs/` | Current `Los_Angeles_discrete_none_MIN` results on GitHub; older `none` results kept locally |
| `datasets/` | Raw data used by upstream preprocessing |
| `maps/` | Retained cluster maps for the `none` and `none_MIN` configurations |
| `CA_simulated_ed/` | Separate simulation workflow; retained with its own README |
| `documents/` | Reference documents |
| `methods/acute_heat_supplement.tex` | Supplementary acute-heat math, methods, calculation code, and Poisson maps |
| `methods/acute_heat_supplement.pdf` | Previously compiled supplement; may lag behind source edits |
| `methods/methodology.tex` | Broader methods draft |
| `analysis_updates/` | Earlier analysis reports and reproducible outputs retained on GitHub |
| `archive/` | Earlier scripts and archived results, retained locally |

The upstream clustering and suppression scripts currently select
`Los_Angeles_discrete_none_MIN`; the current heat scripts read prepared inputs
from `Los_Angeles_discrete_none_HEAT_by_cluster`. They are not an automatic
end-to-end pipeline. Check the configuration and deliberately prepare the
appropriate heat inputs before replacing them. Upstream clustering also uses
ACS/Census services and additional R packages listed in that script.

Compile the supplement from `methods/`; it references both PNGs through
`../Los_Angeles_discrete_none_HEAT_by_cluster/acute_heat_poisson/`. Keep this
relative folder structure when compiling or sharing its source. The built-in
single-file LaTeX preview has previously been unable to load these external
images; a compiler with access to the project files is needed for the full PDF.

## Archive and restore

The cleanup on 2026-10-15 moved loose PReMiuM run traces, older statewide
outputs/profile exports, an unreferenced daily-count file, the automatic
`Rplots.pdf`, and alternative clustering results/maps into
`archive/cleanup_2026-10-15/`. These were outside the current heat workflow and
the selected upstream configuration. No files were deleted.

`archive/cleanup_2026-10-15/move_manifest.csv` records each original path and its
archived path, relative to the project root. Directory entries include all
their contents. To restore an item, move it from `ArchivedPath` back to
`OriginalPath`, creating its parent directory if needed and checking for an
existing destination first. Archived PReMiuM traces may be needed to inspect
old run objects. Older scripts can refer to their original paths and may need
restoration before reuse. A future clustering run may create new traces in the
project root.

Raw data, references, simulation work, active heat results/caches, upstream
`none`/`MIN` results, and RStudio/Git metadata were retained. Their absence from
the current heat script alone was not treated as evidence that they were unused.

## GitHub contents

Large fitted-model files, model caches, and raw PReMiuM run traces are kept locally and excluded from Git. A fresh clone refits models when needed. The obsolete top-level Los_Angeles_discrete_none_HEAT folder is excluded; the current heat workflow uses Los_Angeles_discrete_none_HEAT_by_cluster.

The local `archive/` and older `create_cluster outputs/Los_Angeles_discrete_none/`
results are also excluded from GitHub because current R configurations do not
use them. `analysis_updates/` remains on GitHub. Root-level methods documents
are collected under `methods/`; generated map PDFs remain beside their figures.

