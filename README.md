# Los Angeles heat and ED analysis

Run the R scripts from the project root. This project contains SES clustering,
daily ED-count reconstruction, and Poisson/negative-binomial heat-response models.

## Corrected negative-binomial climate-zone model

The current requested NB specification is fitted by `hri_ed_analysis.R`. It removes
the common acute-heat slope and all 13 SES-by-heat interactions from the saved
`climate_zone_added` design, leaving 26 fixed coefficients. Acute heat enters through
the shared and climate-zone RW2 curves. The model retains the negative-binomial
likelihood, iid ZIP random intercept and exchangeable climate-zone deviations
using `group = climate_zone_group`.
The group index is the 1--5 encoding of the original climate-zone codes.

The shared inputs `Los_Angeles_discrete_none/heat_model_layout.rds` and
`Los_Angeles_discrete_none/heat_grid.rds` hold the preserved row layout, factor
coding, population exposure, fixed-design information, heat grid and curve-node
mappings. Fitting and validation read these compact files directly, without
loading the original `climate_zone_added` fitted model. The layout contains no
visit counts; current counts always come from the reconstructed daily input.
Map colors and themes are preserved separately in
`Los_Angeles_discrete_none/heat_map_style.rds`. The old
`Los_Angeles_discrete_none/climate_zone_added/` results folder has been removed;
the new fit and maps do not read it.

`suppressed_data.R` now assigns one visit to one randomly selected day in each
all-zero ZIP-month with a suppressed monthly total, before applying the raw-zero
branch. With seed 1, exactly seven daily counts change from 0 to 1; all other
counts are unchanged. Run this reconstruction before fitting if rebuilding inputs.
The corrected daily input has been used to refit both NB and matched Poisson
models. The NB maps and forest plots have been regenerated using 2,000 joint
posterior draws for the excess-rate outputs. The preserved `main.pdf` and retained
coefficient-variance comparisons describe the previous fit; `numbered_methods.pdf`
describes the revised model and current precision-prior comparisons.
The precision-prior comparisons have been replaced with runs of the revised model.

```powershell
& 'C:\Program Files\R\R-4.4.3\bin\Rscript.exe' --vanilla hri_ed_analysis.R --preflight
& 'C:\Program Files\R\R-4.4.3\bin\Rscript.exe' --vanilla hri_ed_analysis.R
```

The second command fits the model, samples 2,000 joint posterior draws for excess
rates, exports the three geographic maps and accompanying forest plots, and runs
the validation checks. Saved fits and completed posterior draws are reused.
There is one entry point for this NB workflow:

| Command | Purpose |
| --- | --- |
| `Rscript hri_ed_analysis.R` | Full workflow |
| `Rscript hri_ed_analysis.R --preflight` | Check the model specification before fitting |
| `Rscript hri_ed_analysis.R --fit-only` | Fit and save model summaries |
| `Rscript hri_ed_analysis.R --outputs-only` | Generate outputs from the saved fit |
| `Rscript hri_ed_analysis.R --validate-only` | Check the saved model and outputs |
| `Rscript hri_ed_analysis.R --sensitivity` | Compare precision priors using the corrected model |
| `Rscript hri_ed_analysis.R --likelihood-comparison` | Fit a matched Poisson model and compare DIC, WAIC and zeros with the saved NB baseline; no draws |
| `Rscript hri_ed_analysis.R --coefficient-sensitivity` | Compare NB Normal coefficient variances 100 and 10 against baseline 1000; no draws |
| `Rscript hri_ed_analysis.R --rerun-variance-100` | Restart variance 100 from its fitted solution and compare numerical reproducibility; no draws |

Fitting, output generation, validation and precision-prior sensitivity are all
handled by `hri_ed_analysis.R`. There is no separate validation script.

Results are saved separately in
`Los_Angeles_discrete_none/climate_zone_added_nb_pc/`; new map exports go to
`maps/nb_pc_prior/`, reusing the previous locations. Current precision-prior
comparisons are in `prior_sensitivity_nb/`; its `coefficient_variance/` subfolder
retains the earlier coefficient-prior checks for the previous model.

The four precision priors are `pc.prec`, with `P(SD > 1) = 0.01`. The zone
correlation has a two-sided PC prior with base zero, equal mass on each sign, and
`P(rho > 0.5) = 0.1`, supported on `(-0.25, 1)` for the five-zone exchangeable
matrix. Its internal-scale density includes the transformation Jacobian and is
numerically checked before fitting. NB size uses `pc.mgamma` with parameter 7.
Fixed coefficients retain Normal(0, variance 1000) priors and the intercept is
flat. These are PC hyperpriors, not PC priors on individual fixed coefficients.

Original heat scaling, observation rows, population offsets, RW2 constraints and
heat grid are preserved; reconstructed counts are updated. Full one-SD heat
contrasts use the shared/zone curves without fixed linear heat slopes. Clusters
in the same climate zone share the heat-response curve. The excess-rate reference remains the SES
cluster's monthly mean heat, with equal ZIP-day weights. The individual fixed
baseline cluster effects remain in the model.

The IRR legend is fixed at 1.0--2.0 in 0.2 intervals, using the positive side of
the `climate_zone_added` palette: white at 1.0 and red at 2.0. Both probability
maps retain their corresponding blue--off-white--red `climate_zone_added` palette
with the midpoint at 50% and limits 0--100%.

`Rscript hri_ed_analysis.R --validate-only` checks the reduced fixed design, updated counts, active
fitted hyperpriors, all 42 contrast coefficient vectors, Gaussian sampling
variances, draw counts and model provenance. `validation.json` records the
checks. The correlation-prior check includes INLA's native two-zone PC prior as
an independent benchmark. The model-fit manifest, fitted prior specification,
sampling provenance and validation reports remain with the current results in
`Los_Angeles_discrete_none/climate_zone_added_nb_pc/` after refitting.

The custom exchangeable correlation prior uses the determinant-distance PC
construction described by [Ventrucci et al.](https://arxiv.org/abs/1902.08828),
extended here to both valid correlation branches with equal branch weights.
Its numerical implementation uses INLA's documented
[internal-scale prior table interface](https://sites.google.com/r-inla.org/r-inla/faq).

## Current inputs and results

The active folder is:

```text
Los_Angeles_discrete_none/
  Los Angeles_ACS_zcta_clustered_for_ED.rds
  Los Angeles_daily_ed_no_heat.rds
  heat_model_layout.rds
  heat_grid.rds
  heat_map_style.rds
  Los Angeles_ACS_diss_mats.rds
  Los Angeles_ACS_diss_mats_MIN.rds
  los_angeles_zcta_2020.rds
  create_map.R
  linear_cluster/
  climate_zone_added_nb_pc/
```

The cluster assignments are the 20-covariate, 14-cluster solution used by all
current heat models. The prepared daily ED input includes the seven one-visit
corrections; previously saved fits retain their original counts internally.

| Results folder | Model |
| --- | --- |
| `linear_cluster/` | Original linear acute-heat effects with SES-cluster heat slopes; no climate-zone heat curves |
| `climate_zone_added_nb_pc/` | Current NB baseline without fixed linear heat slopes, using corrected counts; includes summaries, checks and 2,000 excess-rate draws |
| `climate_zone_added_poisson_pc/` | Matched Poisson fit and NB/Poisson DIC, WAIC and zero-count diagnostics; no Poisson maps or joint draws |

Each model folder contains its saved fits, coefficient summaries, WAIC/DIC TXT
summary, AR summaries/posterior draws, and applicable acute-heat tables and plot
objects. Existing model filenames retain the legacy
`Los_Angeles_discrete_none_HEAT_by_cluster` prefix; this is a filename prefix,
not the location of an active folder.

All map scripts read the shared saved Census geometry from
`Los_Angeles_discrete_none/los_angeles_zcta_2020.rds`. Identical model-specific
copies were consolidated into this single input. Historical linear-model plots
remain in `linear_cluster/acute_heat_poisson/` and `linear_cluster/acute_heat_nb/`.

The `climate_zone/` and `climate_zone_zip_iid/` result folders and the
`maps/climate_zone/` exports were removed during cleanup. The older `model_cache/`
was also removed; the current baseline and map workflow read saved results directly.

## Current-model precision-prior sensitivity

The four precision-prior scenarios have been rerun for the revised model with
26 fixed coefficients, no linear heat slopes, and corrected suppressed counts.
They replace the previous precision-prior fits and reports in `prior_sensitivity_nb/`.
The separate coefficient-variance results remain historical and have not been rerun.

Run `Rscript hri_ed_analysis.R --sensitivity`. This is a fit-only comparison: it
does not generate excess-rate posterior draws. It reuses the exact saved current
baseline fit and any completed sensitivity fits, and fits the remaining
alternative settings. All four scenarios preserve 26 fixed coefficients, the same data, population
offsets, ZIP iid intercept, exchangeable zone grouping, heat grid, PC correlation
prior and PC NB-size prior. Only the four precision priors vary jointly:

| Scenario | Prior on the four random-effect precisions |
| --- | --- |
| `pc_baseline` | PC: P(SD > 1) = 0.01 |
| `pc_stronger` | PC: P(SD > 0.5) = 0.01 |
| `pc_weaker` | PC: P(SD > 2) = 0.01 |
| `gamma_comparison` | Gamma(shape = 1, rate = 0.00005), a comparison only |

Sensitivity results are in the top-level `prior_sensitivity_nb/` folder: `comparison.csv`
and `comparison.txt`, all 42 full heat contrasts per scenario in
`heat_contrasts.csv`, and conditional cluster IRRs. Each scenario retains its
saved fit, active priors, summaries and validation reports. The main baseline's
excess-rate draws are preserved; the comparison generates and uses no such draws.
The weaker-PC and Gamma fits were restarted from their fitted solutions; retained
fits passed their curvature checks without negative-eigenvalue corrections.
Stronger/weaker PC alternatives changed heat IRR medians by at most 1.86% and
probabilities by 3.12 percentage points. The Gamma comparison changed these by
up to 6.55% and 10.26 points. No heat interval-exclusion decisions changed.
These assess the specified joint changes
to precision priors; they do not separately vary the correlation or dispersion
prior. The shared and climate-zone RW2 curves provide the full acute-heat response.
Fit-only validation checks the fitted model and INLA contrasts. It skips the
separate excess-rate sampling backend; that check remains enabled when validating
the current model and its sampled excess-rate outputs.

The optional correlation and NB-dispersion sensitivity workflow and its outputs
were removed. The baseline remains the primary model; those two prior groups
were not fully assessed by the retained sensitivity analysis.

## Previous-model coefficient-prior variance sensitivity

These retained results concern the previous 39 non-intercept coefficients and
have not been repeated for the revised model's 25 non-intercept coefficients.

Run `Rscript hri_ed_analysis.R --coefficient-sensitivity`. Results are saved in
`prior_sensitivity_nb/coefficient_variance/`, with new fits in `variance_100/`
and `variance_10/`. The exact current NB baseline with variance 1000 is reused
without copying or refitting its model. The two alternatives are fitted
sequentially and do not save joint-sampling configurations or generate
excess-rate draws.

Only the zero-centered Normal variances for the 39 non-intercept fixed
coefficients change. The flat intercept and all PC precision, correlation and
NB-size priors stay at baseline. Comparison tables include DIC/WAIC, fitted-mean
zero checks, full heat IRRs with credible intervals and exceedance probabilities,
and fixed-effect summaries. Hash checks preserve baseline results, current maps
and methods files. Existing precision-prior sensitivity results are separate.

The coefficient-variance numerical check runs with
`Rscript hri_ed_analysis.R --rerun-variance-100` and saves a separate fit under
`prior_sensitivity_nb/coefficient_variance/variance_100_rerun/`. It initializes
both hyperparameter and latent modes from the original variance-100 solution,
then allows further optimization (`restart = TRUE`, `fixed = FALSE`). The
original fit and sensitivity tables are preserved. The rerun saves DIC/WAIC,
zero diagnostics and complete heat contrasts without generating excess-rate draws.
The warning-producing original variance-100 fit was removed after this check;
the retained clean rerun is now the canonical fit in `variance_100/`, and the
coefficient-sensitivity comparison tables use that retained fit.

## Matched Poisson comparison

Run `Rscript hri_ed_analysis.R --likelihood-comparison`. Results are saved in
`Los_Angeles_discrete_none/climate_zone_added_poisson_pc/`. This reuses the exact
saved NB baseline and fits Poisson counts with the same observations, population
exposures, fixed effects, ZIP iid intercept, shared and zone RW2 terms, and common
priors. It does not generate excess-rate posterior draws or export maps.

`comparison.csv` and `comparison.txt` contain DIC, WAIC and observed versus
expected zeros. Zero diagnostics use posterior mean fitted rates and, for NB,
posterior mean size: these are in-sample plug-in checks, not fully integrated
posterior-predictive tests. Additional tables give zeros by SES cluster, ZIP,
climate zone and month, and count frequencies for 0, 1, 2, 3–4 and 5+ visits.
Hashes verify that the baseline results, current maps and methods files remain
unchanged. The Poisson fit omits joint-sampling configurations and derived heat
contrasts because neither is needed for this comparison.

## Scripts

`hri_ed_analysis.R` is the single main entry point for the current corrected NB
model, output generation and retained precision-prior sensitivity comparisons.
It replaces the disabled historical workflow and the former separate NB script.
The zero diagnostic in its summary tables uses fitted posterior mean rates and
posterior mean NB size; it is a plug-in diagnostic. ZIP SD summaries transform
the posterior precision marginal.

- `hri_ed_analysis.R`: current fitting, output generation and prior sensitivity.
- `acute_heat_irr_maps.R`: computes complete one-SD heat contrasts from INLA
  linear-combination marginals and exports the matching IRR/probability maps.
- `shared_functions.R`: data preparation, fitting, contrasts, and map helpers.
- `cluster_ar.R`: active helper called by `hri_ed_analysis.R` for signed
  cluster-average excess rates and their posterior probabilities, comparing
  observed heat with the SES cluster's monthly mean heat.
- `create_cluster.R`: upstream 20-covariate SES clustering; reads/writes the
  canonical output folder. It retains 2,000 burn-in iterations and 10,000 sweeps.
- `suppressed_data.R`: daily-count reconstruction using the canonical cluster/input
  folder, with the all-zero suppressed-month one-visit override before raw zeros.
- `Los_Angeles_discrete_none/create_map.R`: redraws the current NB maps and forest
  plots; the retained `linear_cluster` plots can be selected explicitly.
- `Los_Angeles_discrete_none/heat_coverage_diagnostics.R`: descriptive heat
  coverage and missing climate-zone report; it does not change model eligibility.
Regenerate the current NB outputs from the saved fit:

```powershell
Rscript hri_ed_analysis.R --outputs-only
```

The map helpers retain support for explicitly selected older model variants.
Fitting, count reconstruction, and clustering are separate actions; the map
scripts do not automatically run those upstream steps.

The separate simulation workflow reads the canonical assignments and writes
its own outputs into `CA_simulated_ed/results/`.

## Maps and interpretation

Map PNGs and PDFs remain under matching folders in `maps/`:
`linear_cluster/`, `climate_zone_added/`, and `nb_pc_prior/` (revised NB model).
The NB maps and three forest plots have been regenerated in that location, each
as PNG and PDF. The fixed 1.0--2.0 IRR legend uses the nearest end color for
estimates outside its display limits; numeric results retain their fitted values.
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

## Current clustering inputs

The current 20-covariate, 14-cluster assignments and their saved dissimilarity
objects remain under `Los_Angeles_discrete_none/` and are used by the current
heat analysis.

The historical clustering-diagnostics script and its `cluster_diagnostics/`
reports were removed during cleanup.

The old archive was removed during cleanup. It contained the raw historical
MIN traces; rerunning those diagnostics requires obtaining those original
trace files again. They are separate from the current clustering inputs and
are not needed for the NB model, its validation, or the current maps.

## Dependencies and Git

Heat scripts use `INLA`, `dplyr`, `ggplot2`, `lubridate`, `sf`,
`tigris`, and `scales`. Upstream clustering and diagnostics require their
additional listed packages. Scripts do not install packages automatically.

Large fitted models, caches, raw MCMC traces and raw datasets are excluded from
Git. Active canonical clustering inputs and summary
files are no longer excluded merely because of the canonical folder name.
Fitted models must be regenerated locally; raw datasets and `datasets/cz_zip.csv` are local inputs. Diagnostic scripts, diagnostic outputs, logs, and caches are excluded from publication.

